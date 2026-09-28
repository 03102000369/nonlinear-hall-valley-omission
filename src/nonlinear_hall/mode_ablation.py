"""Fixed physical hierarchy and observable-weighted leakage diagnostics.

All vectors use Euclidean coordinates after the sqrt(contour-weight) transform
and the orthonormal time-reversal-odd projection. No coupled-response fitting.
"""
import numpy as np
from scipy.linalg import norm, solve, eigh, svdvals
from . import core
from .final_physics import harmonic
from .referee_numerics import factor_solve

NAMES = ('active_reference_current', 'active_valley_distortion',
         'passive_current', 'passive_valley_population')


def basis(system):
    n = system.p.n
    raw = np.column_stack((np.r_[system.g0[:n], np.zeros(n)],
                           harmonic(system, 'A', 0),
                           np.r_[np.zeros(n), system.b[n:]],
                           harmonic(system, 'P', 0)))
    U = np.zeros_like(raw); transform = np.eye(4); scales = []
    for j in range(4):
        v = raw[:, j].copy(); c = transform[:, j].copy()
        for _ in range(2):
            for i in range(j):
                a = U[:, i] @ v
                v -= a * U[:, i]; c -= a * transform[:, i]
        scale = norm(v)
        if scale < 1e-13: raise ValueError('Physical basis is degenerate')
        U[:, j] = v / scale; transform[:, j] = c / scale; scales.append(scale)
    return U, raw, transform, np.array(scales)


def observe(system, g):
    n = system.p.n; b = system.b; h = system.h
    H = h @ g; sig = b @ g; sigA = b[:n] @ g[:n]
    D = h @ b; W = b @ b; WA = b[:n] @ b[:n]
    return dict(H=H, sigma_xx=sig, chi_yxx=-H/2,
                Q=H*W/(sig*D), Q_A=H*WA/(sigA*D), sigma_A=sigA)


def hierarchy(system, eta, omega=0., prepared=None):
    U, raw, transform, scales = basis(system) if prepared is None else prepared
    L = system.L0 + eta * system.K
    A = L + 1j*omega*np.eye(len(L)) if omega else L
    sol = factor_solve(A, np.column_stack((system.b, system.h)), omega)
    g, z = sol[:, 0], sol[:, 1]; full = observe(system, g)
    LU = L @ U; K = U.T @ LU; d = U.T @ system.b
    hs = []; reconstructed = []
    previous = dict.fromkeys(full, 0.)
    for m in range(1, 5):
        km = K[:m, :m] + 1j*omega*np.eye(m) if omega else K[:m, :m]
        x = factor_solve(km, d[:m], omega); gm = U[:, :m] @ x
        row = dict(dimension=m, **observe(system, gm))
        for key in full:
            row['error_'+key] = abs(row[key]-full[key])/abs(full[key])
            row['increment_'+key] = row[key]-previous[key]
        row['distribution_error'] = norm(gm-g)/norm(g)
        row['active_distribution_error'] = norm(gm[:system.p.n]-g[:system.p.n])/norm(g[:system.p.n])
        row['H_error_over_H0'] = abs(row['H']-full['H'])/abs(system.H0)
        row['reduced_equation_residual'] = norm(km@x-d[:m])/norm(d[:m])
        hs.append(row); reconstructed.append(gm); previous = row
    R = LU-U@K; bperp = system.b-U@d
    x4 = U.T @ reconstructed[-1]; residual = system.b-A@reconstructed[-1]
    error_identity = z @ residual
    zperp = z-U@(U.T@z)
    dual_x = factor_solve(K+1j*omega*np.eye(4) if omega else K, U.T@system.h, omega)
    dual_residual = system.h-A@(U@dual_x)
    defect = dict(rho_abs=svdvals(R)[0], rho_rel=svdvals(R)[0]/svdvals(LU)[0],
                  rho_F_abs=norm(R), rho_F_rel=norm(R)/norm(LU),
                  drive_projection_relative=norm(bperp)/norm(system.b),
                  solution_residual_relative=norm(residual)/norm(system.b),
                  leakage_forcing_norm=norm(R@x4), drive_miss_norm=norm(bperp),
                  orthogonality_error=norm(U.T@U-np.eye(4)),
                  transform_error=norm(raw@transform-U)/norm(U),
                  H_error=hs[-1]['error_H'], distribution_error=hs[-1]['distribution_error'],
                  Q_A_error=hs[-1]['error_Q_A'],
                  signed_H_error=full['H']-hs[-1]['H'], adjoint_error=error_identity,
                  adjoint_identity_defect=abs(full['H']-hs[-1]['H']-error_identity),
                  adjoint_bound=norm(z)*norm(residual),
                  orthogonal_adjoint_bound=norm(zperp)*norm(residual),
                  dual_projection_relative=norm(zperp)/norm(z),
                  dual_residual_relative=norm(dual_residual)/norm(system.h),
                  residual_decomposition_defect=norm(residual-(bperp-R@x4)),
                  M4_H_error_over_H0=hs[-1]['H_error_over_H0'])
    return full, hs, defect, dict(U=U,L=L,g=g,g4=reconstructed[-1],R=R,
         residual=residual,dual_residual=dual_residual,bperp=bperp,x4=x4,raw=raw,transform=transform,scales=scales)


def characterize(system, meta, prepared=None):
    U, raw, T, scales = basis(system) if prepared is None else prepared
    n = system.p.n; rows=[]; transforms=[]
    for i, name in enumerate(NAMES):
        v=U[:,i];full=core.lift(v,system.st,n)
        exchange=full.reshape(2,2,n)[:,::-1,:].reshape(-1)
        # physical mirror is TR followed by local qy reflection.
        local_y=v.reshape(2,n)[:,(-np.arange(n))%n].reshape(-1)
        sector='A' if i<2 else 'P';pop=harmonic(system,sector,0)
        current=np.r_[system.b[:n],np.zeros(n)] if i<2 else np.r_[np.zeros(n),system.b[n:]]
        rows.append(dict(**meta,mode=i+1,name=name,sector=sector,
           raw_norm=norm(raw[:,i]),orthogonalization_norm=scales[i],
           time_reversal_parity=-1.,physical_mirror_expectation=-v@local_y,
           local_qy_reflection_expectation=v@local_y,
           equal_angle_valley_exchange_expectation=full@exchange,
           drive_overlap=system.b@v,Hall_overlap=system.h@v,
           normalized_velocity_overlap=current@v/norm(current),
           population_overlap=pop@v,normalized_population_overlap=pop@v/norm(pop)))
        for j in range(4):transforms.append(dict(**meta,raw_mode=j+1,orthonormal_mode=i+1,coefficient=T[j,i]))
    return rows, transforms


def spectrum(system, eta, meta):
    _, _, defect, raw=hierarchy(system,eta)
    rates,V=eigh(raw['L'],check_finite=False)
    b=V.T@system.b;h=V.T@system.h;r=V.T@raw['residual']
    leak=V.T@raw['R']; miss=V.T@raw['bperp'];forcing=leak@raw['x4']
    terms=h*r/rates;power=np.sum(leak**2,axis=1)
    rows=[]
    for j,rate in enumerate(rates):
        rows.append(dict(**meta,eigenmode=j,rate=rate,drive_overlap=b[j],Hall_overlap=h[j],
             leakage_power=power[j],leakage_fraction=power[j]/sum(power),
             drive_miss_overlap=miss[j],leakage_forcing_overlap=forcing[j],residual_overlap=r[j],
             signed_error_contribution=terms[j],absolute_error_contribution=abs(terms[j]),
             full_H_contribution=h[j]*b[j]/rate,
             complement_weight=1-norm(raw['U'].T@V[:,j])**2,
             **{f'leakage_column_{k+1}_overlap':leak[j,k] for k in range(4)}))
    fast=rates>np.median(rates)
    # Weighted rates avoid dependence on arbitrarily rotated degenerate eigenvectors.
    summary=dict(**meta,lambda_min=rates[0],lambda_max=rates[-1],rate_median=np.median(rates),
        leakage_weighted_rate=power@rates/sum(power),
        leakage_fraction_above_median=sum(power[fast])/sum(power),
        sum_signed_error=sum(terms),sum_absolute_error=sum(abs(terms)),
        sum_abs_drive_miss_error=sum(abs(h*miss/rates)),
        sum_abs_leakage_forcing_error=sum(abs(h*forcing/rates)),
        spectral_identity_defect=abs(sum(terms)-defect['signed_H_error']),
        resolvent_bound=norm(system.h)*norm(raw['residual'])/rates[0],
        primal_dual_bound=norm(raw['dual_residual'])*norm(raw['residual'])/rates[0],
        orthogonal_adjoint_bound=defect['orthogonal_adjoint_bound'],
        dual_projection_relative=defect['dual_projection_relative'],
        dual_residual_relative=defect['dual_residual_relative'],
        residual_norm=norm(raw['residual']),
        drive_miss_norm=norm(raw['bperp']),leakage_forcing_norm=norm(raw['R']@raw['x4']),
        residual_cancellation_ratio=norm(raw['residual'])/(norm(raw['bperp'])+norm(raw['R']@raw['x4'])),
        Hall_norm=norm(system.h),signed_H_error=defect['signed_H_error'])
    return rows,summary
