"""Physical mode projection and tolerance controls; the validated solver is unchanged."""
import numpy as np
from scipy.linalg import solve, eigh, eigvalsh, norm
from . import core
from .referee_numerics import factor_solve


def cancellation(H0, S, T, tolerance, floor=1e-12):
    H = H0 + S + T
    a = abs(S) + abs(T)
    scale = max(abs(H0), a, 1e-300)
    conditioned = min(abs(H0), abs(H)) > floor * scale
    actual = abs(S + T) / abs(H) if abs(H) else np.inf
    envelope = a / (abs(H0) - a) if a < abs(H0) else np.inf
    discrete = max(abs(d) / abs(H0 + d) if H0 + d else np.inf
                   for d in [a, -a, abs(S)-abs(T), -abs(S)+abs(T)])
    decision = ('UNDEFINED' if not conditioned else 'UNSAFE' if actual > tolerance
                else 'ROBUSTLY_SAFE' if envelope <= tolerance else 'CANCELLATION_DEPENDENT_SAFE')
    return dict(H0=H0, H=H, S=S, T=T, DeltaH=S+T, a=a,
                B=a/abs(H0) if H0 else np.inf, C=abs(S+T)/a if a else np.nan,
                actual_omission=actual, no_cancellation_envelope=envelope,
                discrete_sign_worst=discrete, opposite_signs=bool(S*T < 0),
                tolerance=tolerance, decision=decision)


def harmonic(system, sector, k, kind="cos"):
    """A weighted angular displacement with imposed time-reversal odd parity."""
    st, n = system.st, system.p.n
    # s*cos[k(theta-pi*1_{s=-})] is time-reversal odd and mirror even.
    trig = np.cos if kind == 'cos' else np.sin
    values = st['s'] * trig(k*(st['theta']-np.pi*(st['s'] < 0)))
    return core.odd_vector(np.sqrt(st['weight'])*values*(st['sector']==sector), st, n)


def physical_basis(system, dimension=4):
    n=system.p.n
    current_A=np.r_[system.g0[:n], np.zeros(n)]
    current_P=np.r_[np.zeros(n), system.b[n:]]
    raw=[current_A, harmonic(system,'A',0), current_P, harmonic(system,'P',0)]
    labels=['active_reference_current','active_valley_distortion','passive_current','passive_valley_population']
    for k in [2,3]:
        for sector in ['A','P']:
            raw.append(harmonic(system,sector,k));labels.append(f'{sector}_cos{k}')
    basis=[]
    for v in raw[:dimension]:
        v=v.copy()
        for _ in range(2):
            for q in basis:v-=q*(q@v)
        if norm(v)<1e-13:raise ValueError('Degenerate physical basis')
        basis.append(v/norm(v))
    return np.column_stack(basis),labels[:dimension]


def full_response(system, eta):
    n=system.p.n;L=system.L0+eta*system.K
    g=factor_solve(L,system.b);gA0=system.g0[:n]
    alpha=(gA0@g[:n])/(gA0@gA0)
    S=(alpha-1)*system.H0;T=system.ha@(g[:n]-alpha*gA0)
    a=abs(S)+abs(T)
    return dict(H=system.h@g,DeltaH=system.h@g-system.H0,S=S,T=T,
                shape_fraction=abs(T)/a if a>1e-12*abs(system.H0) else np.nan,
                scalar_fraction=abs(S)/a if a>1e-12*abs(system.H0) else np.nan,
                actual_omission=abs(system.h@g-system.H0)/abs(system.h@g)),g,L


def reduced_response(system, eta, dimension=4):
    Q,labels=physical_basis(system,dimension)
    L=system.L0+eta*system.K;K=Q.T@L@Q;d=Q.T@system.b;j=Q.T@system.h
    x=solve(K,d,assume_a='pos');G=norm(system.g0[:system.p.n])
    S=j[0]*(x[0]-G)
    # All non-reference active modes are orthogonal to the active reference.
    T=sum(j[i]*x[i] for i in range(1,dimension))
    H=j@x;a=abs(S)+abs(T)
    out=dict(H=H,DeltaH=H-system.H0,S=S,T=T,
             shape_fraction=abs(T)/a if a>1e-12*abs(system.H0) else np.nan,
             scalar_fraction=abs(S)/a if a>1e-12*abs(system.H0) else np.nan,
             actual_omission=abs(H-system.H0)/abs(H),dimension=dimension,
             closure_residual=norm(L@Q-Q@K)/norm(L@Q),
             drive_projection_residual=norm(system.b-Q@d)/norm(system.b))
    # Eliminate every mode except the two physically distinguished active ones.
    B=K[:2,2:];F=K[2:,2:]
    eff=K[:2,:2]-B@solve(F,B.T,assume_a='pos')
    source=B@solve(F,d[2:],assume_a='pos');p=d[:2]-source
    aa,c,dd=eff[0,0],eff[0,1],eff[1,1]
    stiff=dd-c*c/aa;forcing=p[1]-c*p[0]/aa
    out.update(a_rate=aa,c_rate=c,d_rate=dd,shape_stiffness=stiff,
               shape_forcing=forcing,p0=p[0],p1=p[1],passive_source0=source[0],
               passive_source1=source[1],h0=j[0],h1=j[1],G=G,x0=x[0],x1=x[1],
               scalar_amplitude=x[0]-G,shape_amplitude=x[1],
               shape_mode_H=j[1]*x[1],higher_shape_H=T-j[1]*x[1],
               R=abs(T)/abs(S) if S else np.inf)
    return out,dict(Q=Q,K=K,b=d,h=j,x=x,labels=labels)


def characterize_modes(system,eta,label):
    Q,labels=physical_basis(system,8);values,V=eigh(system.L0+eta*system.K)
    drive=V.T@system.b;hall=V.T@system.h;n=system.p.n
    harmonics=[]
    for sector in ['A','P']:
        for kind,orders in [('cos',range(n//2+1)),('sin',range(1,n//2))]:
            for k in orders:
                v=harmonic(system,sector,k,kind);harmonics.append((sector,kind,k,v/norm(v)))
    harmonic_weights=np.abs(np.array([item[3] for item in harmonics])@V)**2
    populations=[harmonic(system,r,0)/norm(harmonic(system,r,0)) for r in ['A','P']]
    rows=[]
    for i,(rate,v) in enumerate(zip(values,V.T)):
        overlaps=harmonic_weights[:,i];z=int(np.argmax(overlaps))
        sector,kind,k,_=harmonics[z]
        # Equal-local-angle valley exchange, distinct from the imposed TR odd parity.
        full=core.lift(v,system.st,n);exchanged=full.reshape(2,2,n)[:,::-1,:].reshape(-1)
        row=dict(label=label,eta=eta,delta_P=system.p.delta_P,xi=system.p.xi,
                 eigenmode=i,eigenvalue=rate,drive_overlap=drive[i],Hall_overlap=hall[i],
                 Hall_residue=drive[i]*hall[i]/rate,active_weight=float(v[:n]@v[:n]),
                 passive_weight=float(v[n:]@v[n:]),time_reversal_parity=-1.,
                 valley_exchange_expectation=float(full@exchanged),
                 dominant_sector=sector,dominant_harmonic_kind=kind,dominant_angular_harmonic=k,
                 mirror_local_y_expectation=float(v@v.reshape(2,n)[:,(-np.arange(n))%n].reshape(-1)),
                 dominant_harmonic_overlap=overlaps[z],
                 valley_population_weight=sum(abs(u@v)**2 for u in populations))
        row.update({name+'_weight':float(abs(Q[:,j]@v)**2) for j,name in enumerate(labels)})
        rows.append(row)
    return rows


def passive_population(system,g):
    ell=core.lift(g,system.st,system.p.n)/np.sqrt(system.st['weight'])
    st=system.st
    N={s:float(np.sum(st['weight'][(st['sector']=='P')&(st['s']==s)]*ell[(st['sector']=='P')&(st['s']==s)])) for s in [1,-1]}
    return N[1],N[-1]


def mirror_even_basis(n):
    q=[]
    for sector in [0,1]:
        for i in range(n//2+1):
            v=np.zeros(2*n);j=(-i)%n
            v[sector*n+i]=1.;v[sector*n+j]=1.;v/=norm(v);q.append(v)
    return np.column_stack(q)
