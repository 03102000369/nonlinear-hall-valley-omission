"""Analytical bands, contour quadrature, reversible collisions and transport.

No finite differences, eigenvalue clipping, or artificial relaxation are used.
All transposes in response functions are bilinear, not Hermitian conjugates.
"""
from dataclasses import dataclass, asdict
import numpy as np
from scipy.linalg import eigh, eigvalsh, solve, norm, svdvals
from scipy.optimize import brentq

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.diag([1., -1.]).astype(complex)


@dataclass(frozen=True)
class Parameters:
    n: int = 64
    mu: float = 1.8
    m_A: float = 1.0
    E_A: float = 0.0
    t_A: float = 0.15
    v_A: float = 1.0
    v_P: float = 1.0
    delta_P: float = 0.8
    eta: float = 1.0
    xi: float = 0.05
    A0: float = 1.0
    q_cutoff: float = 2.5
    zero_rtol: float = 1e-11

    def band(self, r):
        return (self.E_A, self.m_A, self.t_A, self.v_A) if r == "A" else (self.mu-self.delta_P, 0., 0., self.v_P)


def quantities(q, r, s, p):
    q = np.asarray(q)
    E, m, t, v = p.band(r)
    d = np.sqrt(m*m+v*v*np.sum(q*q, axis=-1))
    if np.any(d == 0):
        raise ValueError("Degenerate Dirac node is outside the band gauge domain")
    eps = E+s*t*q[..., 0]+d
    vx = s*t+v*v*q[..., 0]/d
    vy = v*v*q[..., 1]/d
    berry = -s*m*v*v/(2*d**3)
    spin = np.stack((-v*q[..., 1]/d, v*q[..., 0]/d, s*m*np.ones_like(d)/d), axis=-1)
    return eps, vx, vy, berry, spin


def hamiltonian(q, r, s, p):
    E, m, t, v = p.band(r)
    x, y = q
    return (E+s*t*x)*np.eye(2)-v*y*SX+v*x*SY+s*m*SZ


def berry_kubo(q, r, s, p):
    """Independent eigenstate perturbation formula, A=i<u|grad u>."""
    _, _, t, v = p.band(r)
    ev, u = eigh(hamiltonian(q, r, s, p))
    dx = s*t*np.eye(2)+v*SY
    dy = -v*SX
    product = np.vdot(u[:, 1], dx@u[:, 0])*np.vdot(u[:, 0], dy@u[:, 1])
    return -2*np.imag(product)/(ev[1]-ev[0])**2


def contours(p, active_only=False):
    if p.n % 2 or p.n < 8:
        raise ValueError("An even angular grid >=8 is required")
    if p.eta < 0 or p.xi < 0 or p.A0 <= 0:
        raise ValueError("Invalid scattering parameter")
    arrays = {k: [] for k in ("sector", "s", "theta", "qx", "qy", "kx", "ky", "energy", "vx", "vy", "Omega", "weight", "spin")}
    th = 2*np.pi*np.arange(p.n)/p.n
    for r in (("A",) if active_only else ("A", "P")):
        E, m, t, v = p.band(r)
        if v <= 0 or abs(t) >= v or p.mu-E <= m:
            raise ValueError("Required single positive-radius type-I contour does not exist")
        for s in (1, -1):
            radii = []
            for theta in th:
                fun = lambda z: E+s*t*z*np.cos(theta)+np.sqrt(m*m+v*v*z*z)-p.mu
                if fun(p.q_cutoff) <= 0:
                    raise ValueError("Fermi contour reaches prescribed continuum cutoff")
                radii.append(brentq(fun, 0., p.q_cutoff, xtol=5e-15, rtol=1e-14))
            qf = np.asarray(radii)
            q = qf[:, None]*np.stack((np.cos(th), np.sin(th)), axis=1)
            eps, vx, vy, om, spin = quantities(q, r, s, p)
            radial = s*t*np.cos(th)+v*v*qf/np.sqrt(m*m+v*v*qf*qf)
            w = (2*np.pi/p.n)/(2*np.pi)**2*qf/abs(radial)
            if np.min(w) <= 0 or np.max(abs(eps-p.mu)) > 1e-12:
                raise AssertionError("Invalid contour quadrature")
            vals = (np.repeat(r, p.n), np.repeat(s, p.n), th, q[:, 0], q[:, 1], q[:, 0]+(6 if r == "A" else 12)*s, q[:, 1], eps, vx, vy, om, w, spin)
            for key, val in zip(arrays, vals):
                arrays[key].append(val)
    return {k: np.concatenate(v) for k, v in arrays.items()}


def vectors(st):
    sw = np.sqrt(st["weight"])
    return sw*st["vx"], sw*st["Omega"]


def collision(st, p):
    w = st["weight"]
    k = np.stack((st["kx"], st["ky"]), axis=1)
    d2 = np.sum((k[:, None, :]-k[None, :, :])**2, axis=2)
    overlap = (1+st["spin"]@st["spin"].T)/2
    # Preserve tiny signed roundoff in exact antipodal spin overlaps; audit it.
    kernel = getattr(p, 'kernel', 'gaussian')
    if kernel == 'gaussian':
        form = np.exp(-p.xi**2*d2)
    elif kernel == 'lorentzian':
        power = getattr(p, 'kernel_power', 1)
        if power <= 0: raise ValueError('Lorentzian power must be positive')
        form = (1+p.xi**2*d2)**(-power)
    else:
        raise ValueError(f'Unknown scattering kernel: {kernel}')
    # eta multiplies AP transition rates (squared matrix-element weight).
    W = p.A0*np.where(st["sector"][:, None] == st["sector"][None, :], 1., p.eta)*form*overlap
    C = np.diag(W@w)-W*w[None, :]
    L = np.diag(W@w)-np.sqrt(w)[:, None]*W*np.sqrt(w)[None, :]
    return W, C, L


def odd_indices(st, n):
    plus = np.concatenate([np.arange(j, j+n) for j in range(0, len(st["weight"]), 2*n)])
    minus = np.concatenate([j+n+(np.arange(n)+n//2) % n for j in range(0, len(st["weight"]), 2*n)])
    return plus, minus


def odd_matrix(L, st, n):
    a, z = odd_indices(st, n)
    return (L[np.ix_(a,a)]-L[np.ix_(a,z)]-L[np.ix_(z,a)]+L[np.ix_(z,z)])/2


def odd_vector(v, st, n):
    a, z = odd_indices(st, n)
    return (v[a]-v[z])/np.sqrt(2)


def lift(v, st, n):
    a, z = odd_indices(st, n)
    out = np.zeros(len(st["weight"]), dtype=v.dtype)
    out[a], out[z] = v/np.sqrt(2), -v/np.sqrt(2)
    return out


def audit_collision(st, p, W, C, L, spectral=True):
    w, (b,h) = st["weight"], vectors(st)
    sw = np.sqrt(w)
    scale = norm(L)
    ev = eigvalsh(L, check_finite=False) if spectral else None
    vmax = ev[-1] if spectral else norm(L, 2)
    threshold = p.zero_rtol*vmax
    a, z = odd_indices(st, p.n)
    permutation = np.empty(len(w), int)
    permutation[a], permutation[z] = z, a
    odd = odd_matrix(L, st, p.n)
    oe = eigvalsh(odd, check_finite=False)
    vals = {
        "W_min": float(W.min()),
        "W_symmetry": float(norm(W-W.T)/max(norm(W), 1e-300)),
        "collision_symmetry": float(norm(L-L.T)/scale),
        "weighted_transform_error": float(norm(L-sw[:,None]*C/sw[None,:])/scale),
        "detailed_balance": float(norm(w[:,None]*C-(w[:,None]*C).T)/max(norm(w[:,None]*C),1e-300)),
        "charge_residual": float(norm(L@sw)/(scale*norm(sw))),
        "drive_charge_overlap": float(abs(sw@b)/(norm(sw)*norm(b))),
        "drive_sector_overlap": float(max(abs((sw*(st['sector']==r))@b)/(norm(sw*(st['sector']==r))*norm(b)) for r in np.unique(st['sector']))),
        "TR_collision_residual": float(norm(L-L[np.ix_(permutation,permutation)])/scale),
        "drive_even_residual": float(norm(b+b[permutation])/norm(b)),
        "minimum_eigenvalue": float(ev[0]) if spectral else np.nan,
        "lambda_max": float(vmax),
        "zero_threshold": float(threshold),
        "zero_mode_count": int(np.sum(abs(ev)<=threshold)) if spectral else -1,
        "driven_condition": float(oe[-1]/oe[0]),
        "driven_lambda_min": float(oe[0]),
    }
    expected = len(np.unique(st['sector'])) if p.eta == 0 else 1
    for key in ("W_symmetry", "collision_symmetry", "weighted_transform_error", "detailed_balance", "charge_residual", "drive_charge_overlap", "drive_sector_overlap", "TR_collision_residual", "drive_even_residual"):
        if vals[key] > 1e-11:
            raise AssertionError(f"Collision audit failed: {key}={vals[key]}")
    if W.min() < -1e-12*max(W.max(),1) or (spectral and (ev[0] < -1e-11*vmax or vals['zero_mode_count'] != expected)) or oe[0] <= threshold:
        raise AssertionError(f"Collision positivity/zero modes failed: {vals}")
    return vals, ev


def full_solution(L, b, omega=0., zero_rtol=1e-11):
    if omega:
        return solve(L.astype(complex)+1j*omega*np.eye(len(L)), b, assume_a='sym')
    ev, u = eigh(L, check_finite=False)
    keep = ev > zero_rtol*ev[-1]
    if norm(u[:,~keep].T@b) > 1e-10*norm(b):
        raise AssertionError("Drive overlaps a conserved mode")
    return u[:,keep]@((u[:,keep].T@b)/ev[keep])


def observations(st, g):
    b, h = vectors(st)
    act = st['sector']=='A'
    S, H = b@g, h@g
    WA, W = b[act]@b[act], b@b
    SA, D = b[act]@g[act], h@b
    conditioned = abs(D) > 1e-10*norm(h)*norm(b)
    Q = H*W/(S*D) if conditioned and abs(S)>1e-14*norm(b)*norm(g) else np.nan
    QA = H*WA/(SA*D) if conditioned and abs(SA)>1e-14*norm(b[act])*norm(g[act]) else np.nan
    return dict(D_x=D, D_A=h[act]@b[act], D_P=h[~act]@b[~act], sigma_xx=S, sigma_A=SA, H=H, chi_yxx=-H/2,
                W_xx=W, W_A=WA, tau_fit=S/W, Q=Q, Q_A=QA, R_abs=H-S*D/W)


def compute(p, omega=0., independent=True, bounds=True):
    st = contours(p)
    b, h = vectors(st)
    W, C, L = collision(st,p)
    audit, ev = audit_collision(st,p,W,C,L)
    Lo = odd_matrix(L,st,p.n)
    bo, ho = odd_vector(b,st,p.n), odd_vector(h,st,p.n)
    mat = Lo.astype(complex)+1j*omega*np.eye(len(Lo)) if omega else Lo
    go = solve(mat, bo, assume_a='sym' if omega else 'pos')
    g = lift(go, st, p.n)
    residual = norm((L@ g+1j*omega*g if omega else L@g)-b)/norm(b)
    audit['solve_residual'] = float(residual)
    # Independent full-space spectral pseudoinverse/direct solve.
    if independent:
        gf = full_solution(L,b,omega,p.zero_rtol)
        audit['full_projected_error'] = float(norm(gf-g)/norm(g))
        if audit['full_projected_error'] > 1e-9:
            raise AssertionError(f"Full/projected solve disagreement: {audit}")
    n=p.n
    A, B, F = mat[:n,:n], mat[:n,n:], mat[n:,n:]
    # Exact Schur calculation in physical odd subspace; F is invertible even eta=0.
    sol = solve(F, np.column_stack((B.T,bo[n:])), assume_a='sym' if omega else 'pos')
    Sigma, zeta = B@sol[:,:n], B@sol[:,n]
    M = A-Sigma
    ga = solve(M,bo[:n]-zeta,assume_a='sym' if omega else 'pos')
    audit['block_error'] = float(norm(ga-go[:n])/norm(go[:n]))
    if audit['block_error'] > 1e-9 or residual > 1e-10:
        raise AssertionError(f"Block/linear solve failed: {audit}")
    out = observations(st,g)
    out.update(audit)
    out['H_reduced'] = ho[:n]@ga
    out['dos_ratio'] = st['weight'][2*n:].sum()/st['weight'][:2*n].sum()
    d = np.sqrt(p.m_A**2+p.v_A**2*(st['qx'][:2*n]**2+st['qy'][:2*n]**2))
    gap = min(2*d.min(),2*p.delta_P)
    out['interband_threshold'] = gap
    out['rate_interband_ratio'] = audit['lambda_max']/gap
    raw = dict(**st,g=g,b=b,h=h,eigenvalues=ev,odd_eigenvalues=eigvalsh(Lo),Sigma_P=Sigma,zeta_P=zeta)
    if bounds:
        sta = {k:v[:2*n] for k,v in st.items()}
        _,_,L0 = collision(sta,p)
        L0o = odd_matrix(L0,sta,n)
        L0m = L0o.astype(complex)+1j*omega*np.eye(n) if omega else L0o
        g0 = solve(L0m,bo[:n],assume_a='sym' if omega else 'pos')
        H0 = ho[:n]@g0
        Delta = A-L0m
        V = Delta-Sigma
        f = V@g0+zeta
        smin = svdvals(M)[-1] if omega else eigvalsh(M,subset_by_index=[0,0])[0]
        absbound = norm(ho[:n])*norm(f)/smin
        K = solve(L0m,V)
        kappa = svdvals(K)[0]
        neum = norm(ho[:n])*norm(solve(L0m,f))/(1-kappa) if kappa < 1 else np.inf
        diff=abs(out['H']-H0)
        fullscale = norm(ho[:n])*norm(go[:n])
        ratio = diff/abs(out['H']) if abs(out['H'])>1e-10*fullscale else np.nan
        boundrel=absbound/(abs(H0)-absbound) if absbound<abs(H0) else np.inf
        out.update(H_Aonly=H0,sigma_Aonly=bo[:n]@g0,Q_Aonly=H0*(bo[:n]@bo[:n])/((bo[:n]@g0)*(ho[:n]@bo[:n])),Delta_H_abs=diff,delta_H=ratio,
                   bound_abs=absbound,bound_relative_certified=boundrel,bound_neumann_abs=neum,bound_kappa=kappa,
                   Delta_A_norm=norm(Delta,2),Sigma_P_norm=norm(Sigma,2),zeta_P_norm=norm(zeta),schur_min_singular=smin,
                   bound_violation=max(0.,diff-absbound),H0=H0)
        if diff>absbound+1e-10*max(abs(H0),1):
            raise AssertionError("Sufficient Schur bound violated")
        raw.update(g0=g0,Delta_A=Delta)
    return out,raw


def shared_tau(st,tau,omega):
    b,h=vectors(st)
    # Independent numerical matrix solve benchmark rather than evaluating g analytically.
    g=solve((1/tau+1j*omega)*np.eye(len(b)),b)
    return observations(st,g),g


def mode_groups(L,b,h,zero_rtol=1e-11,group_rtol=1e-9):
    ev,u=eigh(L)
    c,d=u.T@b,u.T@h
    groups=[]
    start=0
    while start<len(ev):
        stop=start+1
        while stop<len(ev) and abs(ev[stop]-ev[start])<=group_rtol*ev[-1]:
            stop+=1
        rl=float(c[start:stop]@c[start:stop]); rh=float(d[start:stop]@c[start:stop])
        groups.append(dict(lambda_rate=float(np.mean(ev[start:stop])),multiplicity=stop-start,R_L=rl,R_H=rh,
                           residue_ratio=rh/rl if rl>1e-14*(b@b) else np.nan,zero_mode=abs(ev[start])<=zero_rtol*ev[-1]))
        start=stop
    ah=sum(abs(x['R_H']) for x in groups)
    for x in groups:
        x['longitudinal_fraction']=x['R_L']/(b@b)
        x['hall_abs_fraction']=abs(x['R_H'])/ah if ah else 0.
        x['significant']=not x['zero_mode'] and (x['longitudinal_fraction']>=1e-6 or x['hall_abs_fraction']>=1e-6)
    return groups
