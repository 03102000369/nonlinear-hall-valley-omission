"""Observable-specific revision calculations, separate from immutable legacy data.

Gamma0 denotes the rate scale. A0 is a read-only legacy API alias only.
Weighted coordinates g=sqrt(w)*ell make Euclidean projection physically weighted.
"""
from dataclasses import dataclass, replace, asdict
import numpy as np
from scipy.linalg import solve, eigh, eigvalsh, svdvals, cho_factor, cho_solve, norm
from . import core

@dataclass(frozen=True)
class RevisionParameters:
    n: int = 256
    mu: float = 1.8
    m_A: float = 1.
    E_A: float = 0.
    t_A: float = .15
    v_A: float = 1.
    v_P: float = 1.
    delta_P: float = .8
    eta: float = 1.
    xi: float = .05
    Gamma0: float = .01
    q_cutoff: float = 2.5
    zero_rtol: float = 1e-11
    m_P: float = 0.
    rotation: float = 0.
    internal_model: str = 'laboratory'
    kernel: str = 'gaussian'
    kernel_power: int = 1

    @property
    def A0(self):
        """Compatibility with preserved legacy band/contour utilities."""
        return self.Gamma0

    def band(self, r):
        return (self.E_A,self.m_A,self.t_A,self.v_A) if r=='A' else (self.mu-self.delta_P,self.m_P,0.,self.v_P)


def contour_states(p):
    """Positive quadratic radii and numerical eigenspinors, independently of brentq."""
    if p.n<8 or p.n%2 or p.Gamma0<=0 or p.xi<0 or p.eta<0:
        raise ValueError('Invalid grid or disorder setting')
    th=2*np.pi*np.arange(p.n)/p.n
    fields=['sector','s','theta','qx','qy','kx','ky','energy','vx','vy','Omega','weight','spinor','spinor_lab','rotation_matrix']
    out={k:[] for k in fields}
    for r in ['A','P']:
        E,m,t,v=p.band(r);delta=p.mu-E
        if v<=0 or abs(t)>=v or delta<=abs(m):
            raise ValueError('Origin-enclosing conduction contour requires delta>|m| and |t|<v')
        U=np.cos(p.rotation/2)*np.eye(2)-1j*np.sin(p.rotation/2)*core.SX if r=='P' else np.eye(2,dtype=complex)
        for s in [1,-1]:
            q=(-delta*s*t*np.cos(th)+np.sqrt(v*v*(delta*delta-m*m)+m*m*t*t*np.cos(th)**2))/(v*v-t*t*np.cos(th)**2)
            xy=q[:,None]*np.array([np.cos(th),np.sin(th)]).T
            eps,vx,vy,om,spin=core.quantities(xy,r,s,p)
            mats=np.array([core.hamiltonian(x,r,s,p) for x in xy])
            ev,vec=np.linalg.eigh(mats);u=vec[:,:,1];up=u@U.T
            d=np.sqrt(m*m+v*v*q*q)
            w=q/(np.abs(s*t*np.cos(th)+v*v*q/d)*2*np.pi*p.n)
            values=[np.repeat(r,p.n),np.repeat(s,p.n),th,xy[:,0],xy[:,1],xy[:,0]+s*(6 if r=='A' else 12),xy[:,1],eps,vx,vy,om,w,u,up,np.repeat(U[None,:,:],p.n,axis=0)]
            for k,x in zip(fields,values):out[k].append(x)
    st={k:np.concatenate(v) for k,v in out.items()}
    return st


def transition_parts(st,p):
    u=st['spinor_lab']
    if p.internal_model=='gauge':
        # <U_i u_i| U_i I U_j^dagger |U_j u_j>; transform ALL disorder blocks.
        u=np.einsum('nji,nj->ni',st['rotation_matrix'].conj(),u)
    elif p.internal_model!='laboratory':raise ValueError('Unknown internal/disorder scenario')
    overlaps=np.abs(u.conj()@u.T)**2
    xy=np.column_stack([st['kx'],st['ky']]);d2=np.sum((xy[:,None,:]-xy[None,:,:])**2,axis=2)
    form=np.exp(-p.xi**2*d2) if p.kernel=='gaussian' else (1+p.xi**2*d2)**(-p.kernel_power)
    rate=p.Gamma0*overlaps*form
    same=st['sector'][:,None]==st['sector'][None,:]
    return rate*same,rate*(~same)


def laplacian(rate,w):
    sw=np.sqrt(w)
    return np.diag(rate@w)-sw[:,None]*rate*sw[None,:]


def factor_solve(mat,rhs,omega=0.):
    if omega:return solve(mat,rhs,assume_a='sym',check_finite=False)
    return cho_solve(cho_factor(mat,check_finite=False),rhs,check_finite=False)


def certificate(E,H0):
    return E/(abs(H0)-E) if E<abs(H0) else np.inf


def relative(a,b):
    return abs(a-b)/abs(a) if abs(a)>1e-14 else np.nan


class PreparedSystem:
    """Fixed bands/range; exact affine collision matrices L0+eta K cached across eta."""
    def __init__(self,p):
        self.p=p;self.st=st=contour_states(p);n=p.n;w=st['weight'];sw=np.sqrt(w)
        r0,rx=transition_parts(st,p);full0=laplacian(r0,w);fullK=laplacian(rx,w)
        self.L0=core.odd_matrix(full0,st,n);self.K=core.odd_matrix(fullK,st,n)
        self.b=core.odd_vector(sw*st['vx'],st,n);self.by=core.odd_vector(sw*st['vy'],st,n)
        self.h=core.odd_vector(sw*st['Omega'],st,n)
        self.g0=factor_solve(self.L0,self.b);self.g1=-factor_solve(self.L0,self.K@self.g0)
        self.gy0=factor_solve(self.L0,self.by)
        self.Aref=self.L0[:n,:n];self.ha=self.h[:n];self.ba=self.b[:n]
        self.H0=self.ha@self.g0[:n]
        self.rate_envelope_unit=p.Gamma0*w.sum()/2
        dA=np.sqrt(p.m_A**2+p.v_A**2*(st['qx'][:2*n]**2+st['qy'][:2*n]**2))
        self.gap=min(2*dA.min(),2*p.delta_P)
        plus,minus=core.odd_indices(st,n);perm=np.empty(4*n,int);perm[plus]=minus;perm[minus]=plus
        def symmetry(m):return norm(m-m.T)/max(norm(m),1e-300)
        def charge(m):return norm(m@sw)/max(norm(m)*norm(sw),1e-300)
        self.base_checks=dict(collision_symmetry=max(symmetry(full0),symmetry(fullK)),charge_residual=max(charge(full0),charge(fullK)),
            parity_residual=max(norm(m-m[np.ix_(perm,perm)])/norm(m) for m in [full0,fullK]),
            eigenspinor_energy_error=float(max(abs(st['energy']-p.mu))),
            origin_enclosure_margin=min(p.mu-p.E_A-abs(p.m_A),p.delta_P-abs(p.m_P)),
            contour_min=float(np.hypot(st['qx'],st['qy']).min()),qF_over_cutoff=float(np.hypot(st['qx'],st['qy']).max()/p.q_cutoff),
            pocket_separation=min(abs(sa*(6 if ra=='A' else 12)-sb*(6 if rb=='A' else 12))-float(np.hypot(st['qx'],st['qy'])[(st['sector']==ra)&(st['s']==sa)].max())-float(np.hypot(st['qx'],st['qy'])[(st['sector']==rb)&(st['s']==sb)].max()) for i,(ra,sa) in enumerate([('A',1),('A',-1),('P',1),('P',-1)]) for rb,sb in [('A',1),('A',-1),('P',1),('P',-1)][i+1:]))
        assert self.base_checks['collision_symmetry']<1e-12 and self.base_checks['charge_residual']<1e-12
        assert self.base_checks['parity_residual']<1e-11
        self.assembly_bytes=sum(x.nbytes for x in [r0,rx,full0,fullK,self.L0,self.K])

    def evaluate(self,eta=None,omega=0.,details=False):
        p=self.p;n=p.n;eta=p.eta if eta is None else eta
        L=self.L0+eta*self.K;mat=L+1j*omega*np.eye(2*n) if omega else L
        g=factor_solve(mat,self.b,omega)
        gy=factor_solve(L+2j*omega*np.eye(2*n) if omega else L,self.by,2*omega)
        ref=self.Aref+1j*omega*np.eye(n) if omega else self.Aref
        g0=factor_solve(ref,self.ba,omega) if omega else self.g0[:n]
        H0=self.ha@g0;Ha=self.ha@g[:n];Hp=self.h[n:]@g[n:];H=Ha+Hp
        alpha=np.vdot(g0,g[:n])/np.vdot(g0,g0);resid=g[:n]-alpha*g0
        ds=(alpha-1)*H0;dh=self.ha@resid;dt=Ha-H0
        budget=abs(ds)+abs(dh);floor=1e-12*max(abs(H0),norm(self.ha)*norm(g0))
        fs=abs(ds)/budget if budget>floor else np.nan;fh=abs(dh)/budget if budget>floor else np.nan;cancel=abs(dt)/budget if budget>floor else np.nan
        A,B,F=mat[:n,:n],mat[:n,n:],mat[n:,n:]
        fx=factor_solve(F,np.column_stack([B.T,self.b[n:]]),omega)
        Sigma=B@fx[:,:n];zeta=B@fx[:,n];M=A-Sigma;Delta=A-ref
        f=Delta@g0-Sigma@g0+zeta
        mz=factor_solve(M,np.column_stack([self.ba-zeta,self.ha]),omega)
        ga,z=mz[:,0],mz[:,1]  # M is complex symmetric, so this also solves M^T z=h.
        exactadj=-z@f;outs=-z@(Delta@g0);src=-z@zeta;feed=z@(Sigma@g0)
        smin=svdvals(M)[-1] if omega else eigvalsh(M,subset_by_index=[0,0],check_finite=False)[0]
        E=norm(self.ha)*norm(f)/smin;Ez=norm(z)*norm(f)
        Sx=self.b@g;Sy=self.by@gy
        Sx0=self.ba@g0;Sy0=self.by[:n]@factor_solve(self.Aref+2j*omega*np.eye(n) if omega else self.Aref,self.by[:n],2*omega)
        chi=-H/2;chi0=-H0/2;sigE=-chi/Sy;sigE0=-chi0/Sy0;sigJ=sigE/Sx**2;sigJ0=sigE0/Sx0**2
        eigen= eigvalsh(L,check_finite=False)
        row=dict(eta=eta,delta_P=p.delta_P,xi=p.xi,n=n,Gamma0=p.Gamma0,m_P=p.m_P,rotation=p.rotation,internal_model=p.internal_model,kernel=p.kernel,omega=omega,
          H=H,H0=H0,H_active=Ha,H_passive_direct=Hp,H_readout_zero=Ha,passive_to_active_readout_norm=norm(self.h[n:])/norm(self.ha),passive_to_active_max_curvature=float(max(abs(self.st['Omega'][2*n:]))/max(abs(self.st['Omega'][:2*n]))),Delta_H_total=dt,Delta_H_scalar=ds,Delta_H_shape=dh,
          alpha=alpha,shape_residual=norm(resid)/norm(g[:n]),fraction_scalar=fs,fraction_shape=fh,fraction_absolute_shape=fh,cancellation_factor=cancel,
          component_budget=budget,component_budget_over_H0=budget/abs(H0),zero_correction=budget<=floor,
          actual_error=relative(H,H0),active_omission_error=relative(Ha,H0),scalar_shape_identity=abs(dt-ds-dh),
          adjoint_identity=abs(dt-exactadj),block_error=norm(ga-g[:n])/norm(g[:n]),collision_identity=abs(dt-outs-src-feed),
          correction_scattering_out=outs,correction_passive_source=src,correction_feedback=feed,
          bound_singular_abs=E,bound_adjoint_abs=Ez,certificate_singular=certificate(E,H0),certificate_adjoint=certificate(Ez,H0),
          bound_singular_violation=max(0.,abs(dt)-E),bound_adjoint_violation=max(0.,abs(dt)-Ez),
          sigma_xx=Sx,sigma_yy=Sy,sigma_xx_Aonly=Sx0,sigma_yy_Aonly=Sy0,chi_yxx=chi,chi_yxx_Aonly=chi0,
          fixed_E_signal=sigE,fixed_E_Aonly=sigE0,fixed_j_signal=sigJ,fixed_j_Aonly=sigJ0,
          H_omission=relative(H,H0),chi_omission=relative(chi,chi0),sigma_xx_omission=relative(Sx,Sx0),sigma_yy_omission=relative(Sy,Sy0),fixed_E_omission=relative(sigE,sigE0),fixed_j_omission=relative(sigJ,sigJ0),
          sigma_xy_residual=abs(self.by@g)/max(abs(Sx),1e-300),
          collision_dominance=abs(dt)/(abs(dt)+abs(Hp)) if abs(dt)+abs(Hp)>floor else np.nan,
          driven_condition=eigen[-1]/eigen[0],minimum_driven_eigenvalue=eigen[0],solve_residual=norm(mat@g-self.b)/norm(self.b),
          rate_interband_upper=max(1.,eta)*self.rate_envelope_unit/self.gap,**self.base_checks)
        if not omega:
            H1=self.ha@(self.g0[:n]+eta*self.g1[:n])
            row.update(H_first_order=H1,weak_omission_estimate=relative(H1,H0),weak_response_error=relative(Ha,H1))
        good=row['origin_enclosure_margin']>0 and row['contour_min']>0 and row['qF_over_cutoff']<1 and row['pocket_separation']>0 and row['rate_interband_upper']<=.05 and row['solve_residual']<1e-10 and row['driven_condition']<1e8 and eigen[0]>0
        row['valid_physical_numerical']=bool(good)
        scale=max(abs(H0),norm(self.ha)*norm(g0),1.)
        for key in ['scalar_shape_identity','adjoint_identity','collision_identity','bound_singular_violation','bound_adjoint_violation']:
            if row[key]>1e-11*scale:raise AssertionError((key,row[key],scale))
        if row['block_error']>1e-10:raise AssertionError('Schur mismatch')
        if details:return row,dict(g=g,g0=g0,residual=resid,M=M,f=f,z=z,Lo=L,Sigma=Sigma,zeta=zeta,Delta=Delta)
        return row


def classify(r,thresholds):
    if r['component_budget_over_H0']<=thresholds['negligible_budget'] and r['shape_residual']<=thresholds['negligible_shape']:return 'negligible-redistribution'
    if r['Delta_H_scalar']*r['Delta_H_shape']<0 and r['cancellation_factor']<=thresholds['cancellation_max'] and r['component_budget_over_H0']>thresholds['negligible_budget']:return 'strong-cancellation'
    if r['fraction_scalar']>=thresholds['dominant_fraction']:return 'scalar-dominated'
    if r['fraction_shape']>=thresholds['dominant_fraction']:return 'shape-dominated'
    return 'mixed'


def weak_series(system,etas,omega=0.):
    """Neumann and positive-operator remainder bounds; spectral radius is diagnostic."""
    from scipy.linalg import eigvals
    L0=system.L0;K=system.K;b=system.b;h=system.h
    mat=L0+1j*omega*np.eye(len(L0)) if omega else L0
    RK=factor_solve(mat,K,omega);g0=factor_solve(mat,b,omega);g1=-RK@g0
    norm_unit=svdvals(RK)[0]
    spectral_unit=max(abs(eigvals(RK)))
    drive_unit=norm(RK@g0)/norm(g0)
    h0=h@g0
    if not omega:
        eig,u=eigh(L0);ih=(u.T@h)/np.sqrt(eig);ib=(u.T@b)/np.sqrt(eig)
        energy_prefactor=norm(ih)*norm(ib)
    rows=[]
    for eta in etas:
        g=factor_solve(mat+eta*K,b,omega);hfull=h@g;h1=h@(g0+eta*g1)
        rho=eta*norm_unit;rho_s=eta*spectral_unit
        err=abs(hfull-h1)
        E=norm(h)*eta**2*norm(RK@RK@g0)/(1-rho) if rho<1 else np.inf
        energy=energy_prefactor*rho_s**2/(1+rho_s) if not omega else np.nan
        if np.isfinite(E) and err>E+1e-11*max(abs(h0),1):raise AssertionError('Neumann remainder violation')
        if not omega and err>energy+1e-11*max(abs(h0),1):raise AssertionError('Energy remainder violation')
        rows.append(dict(eta=eta,omega=omega,H=hfull,H0=h0,H_first_order=h1,first_order_error=relative(hfull,h1),absolute_remainder=err,
          rho_norm=rho,rho_spectral=rho_s,rho_drive=eta*drive_unit,remainder_neumann=E,remainder_energy=energy,
          remainder_certificate=certificate(E,h1),neumann_violation=max(0.,err-E),energy_violation=max(0.,err-energy) if not omega else np.nan))
    return rows
