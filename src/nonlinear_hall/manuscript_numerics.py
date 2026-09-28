"""Independent audit and manuscript diagnostics; legacy model remains intact."""
from dataclasses import dataclass, replace
import hashlib
import numpy as np
from scipy.linalg import eigh, solve, norm
from . import core

@dataclass(frozen=True)
class ManuscriptParameters(core.Parameters):
    kernel: str = 'gaussian'
    kernel_power: int = 1


def checksum(*arrays):
    h=hashlib.sha256()
    for a in arrays:
        a=np.ascontiguousarray(a)
        h.update(str(a.shape).encode());h.update(str(a.dtype).encode());h.update(a.tobytes())
    return h.hexdigest()


def independent_system(p, active_only=False):
    """No production contour, band, overlap or collision calls; no cache.

    Positive quadratic roots, 2x2 eigenspinors and |u_i^dagger u_j|^2.
    The analytic derivative of the eigenenergy supplies the contour Jacobian.
    """
    n=p.n; th=2*np.pi*np.arange(n)/n
    vals={k:[] for k in ['sector','s','theta','qx','qy','kx','ky','energy','vx','vy','Omega','weight','spinor']}
    for sector in (['A'] if active_only else ['A','P']):
        E,m,t,v=(p.E_A,p.m_A,p.t_A,p.v_A) if sector=='A' else (p.mu-p.delta_P,0.,0.,p.v_P)
        delta=p.mu-E
        for s in [1,-1]:
            a=v*v-t*t*np.cos(th)**2
            q=(-delta*s*t*np.cos(th)+np.sqrt(v*v*(delta*delta-m*m)+m*m*t*t*np.cos(th)**2))/a
            x,y=q*np.cos(th),q*np.sin(th)
            ds=np.sqrt(m*m+v*v*q*q)
            spinors=[];energies=[];vx=[];vy=[];om=[]
            dx=np.array([[s*t,-1j*v],[1j*v,s*t]],complex)
            dy=np.array([[0,-v],[-v,0]],complex)
            for xx,yy in zip(x,y):
                mat=np.array([[E+s*t*xx+s*m,-v*yy-1j*v*xx],[-v*yy+1j*v*xx,E+s*t*xx-s*m]])
                ev,u=eigh(mat); up,lo=u[:,1],u[:,0]
                spinors.append(up);energies.append(ev[1])
                vx.append(np.vdot(up,dx@up).real);vy.append(np.vdot(up,dy@up).real)
                om.append(-2*np.imag(np.vdot(up,dx@lo)*np.vdot(lo,dy@up))/(ev[1]-ev[0])**2)
            w=q/(abs(s*t*np.cos(th)+v*v*q/ds)*2*np.pi*n)
            for key,value in zip(vals,[np.repeat(sector,n),np.repeat(s,n),th,x,y,x+s*(6 if sector=='A' else 12),y,energies,vx,vy,om,w,spinors]):vals[key].append(np.asarray(value))
    st={k:np.concatenate(a) for k,a in vals.items()}
    k=np.column_stack([st['kx'],st['ky']]); d2=np.sum((k[:,None]-k[None,:])**2,axis=2)
    overlap=abs(st['spinor'].conj()@st['spinor'].T)**2
    f=np.exp(-p.xi*p.xi*d2) if p.kernel=='gaussian' else (1+p.xi*p.xi*d2)**(-p.kernel_power)
    W=p.A0*np.where(st['sector'][:,None]==st['sector'][None,:],1.,p.eta)*f*overlap
    w=st['weight'];sw=np.sqrt(w)
    C=np.diag(W@w)-W*w[None,:]
    L=sw[:,None]*C/sw[None,:]
    b=sw*st['vx'];h=sw*st['Omega']
    ev,u=eigh(L);keep=ev>p.zero_rtol*ev[-1]
    g=u[:,keep]@((u[:,keep].T@b)/ev[keep])
    return st,W,C,L,ev,u,b,h,g


def independent_observation(p):
    st,W,C,L,ev,u,b,h,g=independent_system(p)
    act=st['sector']=='A'; S=b@g;H=h@g;SA=b[act]@g[act];D=h@b
    _,_,_,_,_,_,ba,ha,ga=independent_system(p,True)
    H0=ha@ga
    row=dict(D_x=D,sigma_xx=S,H=H,Q=H*(b@b)/(S*D),Q_A=H*(b[act]@b[act])/(SA*D),delta_H=abs(H-H0)/abs(H),
             H_Aonly=H0,total_states=len(b),matrix_dimension=len(L),
             contour_checksum=checksum(st['qx'],st['qy']),weight_checksum=checksum(st['weight']),
             matrix_checksum=checksum(L),spinor_checksum=checksum(st['spinor']),
             eigenvector_checksum=checksum(u),contour_energy_residual=max(abs(st['energy']-p.mu)),
             solve_residual=norm(L@g-b)/norm(b),minimum_eigenvalue=ev[0],lambda_max=ev[-1])
    raw=dict(**st,W=W,C=C,L=L,eigenvalues=ev,eigenvectors=u,b=b,h=h,g=g)
    return row,raw


def perturbation_system(p,omega=0.):
    """K assembled solely from AP edge rates, independently of L1-L0."""
    p0=replace(p,eta=0.); st=core.contours(p0);b,h=core.vectors(st)
    _,_,L0=core.collision(st,p0)
    k=np.column_stack([st['kx'],st['ky']]); d2=np.sum((k[:,None]-k[None,:])**2,axis=2)
    f=np.exp(-p.xi**2*d2) if p.kernel=='gaussian' else (1+p.xi**2*d2)**(-p.kernel_power)
    cross=st['sector'][:,None]!=st['sector'][None,:]
    W1=p.A0*cross*f*(1+st['spin']@st['spin'].T)/2
    sw=np.sqrt(st['weight']); K=np.diag(W1@st['weight'])-sw[:,None]*W1*sw[None,:]
    Lo=core.odd_matrix(L0,st,p.n);Ko=core.odd_matrix(K,st,p.n)
    bo=core.odd_vector(b,st,p.n);ho=core.odd_vector(h,st,p.n)
    Rmat=Lo+1j*omega*np.eye(len(Lo)) if omega else Lo
    g0=solve(Rmat,bo);g1=-solve(Rmat,Ko@g0)
    n=p.n; a=Rmat[:n,:n]; f0=Rmat[n:,n:];crossblock=Ko[:n,n:]
    out=-ho[:n]@solve(a,Ko[:n,:n]@g0[:n])
    src=-ho[:n]@solve(a,crossblock@g0[n:])
    feedback2=ho[:n]@solve(a,crossblock@solve(f0,crossblock.T@g0[:n]))
    return dict(st=st,Lo=Lo,Ko=Ko,bo=bo,ho=ho,g0=g0,g1=g1,H0=ho@g0,H1=ho@g1,
                active_out=out,passive_source=src,feedback_first=0.,feedback_second=feedback2)


def shape_diagnostic(raw,n):
    g=raw['g'][:2*n]; g0=core.lift(raw['g0'],{k:v[:2*n] for k,v in raw.items() if k in ['weight']},n)
    alpha=np.vdot(g0,g).real/np.vdot(g0,g0).real
    diff=g-alpha*g0
    return float(norm(diff)/norm(g)),float(alpha),g0,diff


def validity(p,row,raw,converged=False):
    q=np.hypot(raw['qx'],raw['qy']);maxq=float(q.max())
    disks=[]
    for r in ['A','P']:
        for s in [1,-1]:
            m=(raw['sector']==r)&(raw['s']==s)
            disks.append((s*(6 if r=='A' else 12),float(q[m].max())))
    separation=min(abs(a[0]-b[0])-a[1]-b[1] for i,a in enumerate(disks) for b in disks[i+1:])
    gap=row['interband_threshold'];r=row['lambda_max']/gap
    checks=dict(contour_exists=p.mu-p.E_A>p.m_A and p.delta_P>0,
                positive_q=bool(q.min()>0),continuum_valid=maxq<p.q_cutoff,
                pockets_separated=separation>0,type_I=abs(p.t_A)<p.v_A and p.v_P>0,
                interband_positive=gap>0,weak_scattering=r<=.05,
                condition_valid=row['driven_condition']<=1e8,solver_valid=row['solve_residual']<=1e-10,
                response_converged=bool(converged))
    return dict(**checks,valid_for_main_text=all(checks.values()),r_scatt=r,Gamma_relevant=row['lambda_max'],
                Delta_interband_min=gap,qF_min=float(q.min()),qF_max=maxq,qF_over_cutoff=maxq/p.q_cutoff,
                pocket_separation_min=separation,invalid_reasons=';'.join(k for k,v in checks.items() if not v))


def common_ratio(groups):
    sig=[m for m in groups if m['significant']]
    l=np.array([m['R_L'] for m in sig]);h=np.array([m['R_H'] for m in sig])
    r=(l@h)/(l@l)
    return float(r),float(norm(h-r*l)/norm(h))


def compare_rows(lo,hi,keys=('D_x','sigma_xx','H','Q','Q_A','delta_H','delta_Q_A','shape_residual')):
    out=[]
    for k in keys:
        a,b=lo[k],hi[k];small=max(abs(a),abs(b))<1e-8
        error=abs(a-b) if small else abs(a-b)/max(abs(b),1e-300)
        tol=1e-10 if small else .01
        out.append(dict(observable=k,error=float(error),metric='absolute' if small else 'relative',threshold=tol,status='PASS' if error<=tol else 'FAIL',value_from=a,value_to=b))
    return out
