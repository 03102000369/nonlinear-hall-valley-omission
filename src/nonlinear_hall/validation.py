from dataclasses import replace
import numpy as np
from scipy.linalg import norm, eigvalsh
from .core import *


def check(rows,name,value,tolerance,notes='',lower=False):
    passed=value>=tolerance if lower else abs(value)<=tolerance
    rows.append(dict(test=name,measured_value=float(value),tolerance=float(tolerance),criterion='>=' if lower else 'abs <=',status='PASS' if passed else 'FAIL',notes=notes))
    if not passed:raise AssertionError(f'{name}: {value}, tolerance {tolerance}')


def band_checks(p,rows):
    T,M=1j*SY,1j*SX
    check(rows,'T_squared_minus_one',norm(T@T.conj()+np.eye(2)),1e-12)
    eigerr=tr=mir=omtr=velerr=0.
    berry_rows=[]
    for r in ['A','P']:
        for s in [1,-1]:
            for q in [np.array([.37,.21]),np.array([-.8,.41]),np.array([1.3,-.57]),np.array([-.23,-1.1])]:
                H=hamiltonian(q,r,s,p); E,m,t,v=p.band(r)
                rad=np.sqrt(m*m+v*v*(q@q))
                exact=E+s*t*q[0]+np.array([-rad,rad])
                eigerr=max(eigerr,np.max(abs(eigvalsh(H)-exact)))
                tr=max(tr,norm(T@H.conj()@T.conj().T-hamiltonian(-q,r,-s,p)))
                mir=max(mir,norm(M@H@M.conj().T-hamiltonian(q*np.array([-1,1]),r,-s,p)))
                om=quantities(q,r,s,p)[3]; om2=quantities(-q,r,-s,p)[3]
                omtr=max(omtr,abs(om+om2)); num=berry_kubo(q,r,s,p)
                _,u=np.linalg.eigh(H)
                vx=np.vdot(u[:,1],(s*t*np.eye(2)+v*SY)@u[:,1]).real
                vy=np.vdot(u[:,1],(-v*SX)@u[:,1]).real
                vxa,vya=quantities(q,r,s,p)[1:3]
                velerr=max(velerr,abs(vx-vxa),abs(vy-vya))
                berry_rows.append(dict(sector=r,s=s,qx=q[0],qy=q[1],analytical=om,kubo=num,absolute_error=abs(om-num),relative_error=abs((om-num)/om) if om else np.nan))
    for name,val in [('Hamiltonian_eigenvalues',eigerr),('time_reversal',tr),('mirror',mir),('Berry_time_reversal',omtr),('Hellmann_Feynman_velocity',velerr)]:check(rows,name,val,1e-12)
    check(rows,'Berry_Kubo_absolute',max(x['absolute_error'] for x in berry_rows),1e-12)
    check(rows,'Berry_Kubo_relative',max(x['relative_error'] for x in berry_rows if x['sector']=='A'),1e-10)
    st=contours(p); b,h=vectors(st); mask=st['sector']=='P'
    check(rows,'passive_Berry_curvature',np.max(abs(st['Omega'][mask])),1e-14)
    check(rows,'passive_dipole',abs(h[mask]@b[mask]),1e-14)
    check(rows,'equilibrium_Hall_contour_cancellation',abs(st['weight']@st['Omega']),1e-12,'TR also cancels the occupied-sea integral pairwise')
    check(rows,'forbidden_D_y',abs(st['weight']@(st['Omega']*st['vy'])),1e-12)
    check(rows,'positive_weights',np.min(st['weight']),0.,lower=True)
    circular=contours(replace(p,t_A=0.))
    expected=2*(p.mu-p.E_A)/(2*np.pi*p.v_A**2)
    check(rows,'circular_DOS',abs(circular['weight'][:2*p.n].sum()-expected)/expected,1e-12)
    # First-order-in-tilt analytic dipole of the active valley pair.
    # D_pair = -3 t m (mu^2-m^2)/(4 pi mu^4)+O(t^3), v=1,E=0.
    small=replace(p,t_A=1e-5)
    bs,hs=vectors(contours(small))
    d=small.mu-small.E_A
    predicted=-3*small.t_A*small.m_A*(d*d-small.m_A**2)/(4*np.pi*d**4)
    check(rows,'small_tilt_dipole_analytic',abs((hs@bs-predicted)/predicted),1e-7,'Independent expansion through first order; truncation O(t^2) relative')
    return st,berry_rows


def tau_checks(st,rows):
    result=[]; b,h=vectors(st)
    for tau in [.3,1.,2.7]:
        for w in [0.,.017,.31]:
            out,g=shared_tau(st,tau,w); exact=tau/(1+1j*w*tau)
            err=norm(g-exact*b)/norm(exact*b)
            herr=abs(out['H']-exact*(h@b))/abs(exact*(h@b))
            qerr=abs(out['Q']-1)
            result.append(dict(tau=tau,omega=w,H_real=out['H'].real,H_imag=out['H'].imag,H_exact_real=(exact*(h@b)).real,H_exact_imag=(exact*(h@b)).imag,g_relative_error=err,H_relative_error=herr,Q_error=qerr))
    check(rows,'constant_tau_solution',max(x['g_relative_error'] for x in result),1e-10)
    check(rows,'constant_tau_H',max(x['H_relative_error'] for x in result),1e-10)
    check(rows,'constant_tau_Q',max(x['Q_error'] for x in result),1e-10)
    return result


def collision_checks(p,rows):
    output=[]
    for eta in [0.,1e-3,1.,10.]:
        for xi in [0.,.05,.1]:
            pp=replace(p,eta=eta,xi=xi)
            st=contours(pp); W,C,L=collision(st,pp); audit,ev=audit_collision(st,pp,W,C,L)
            output.append(dict(eta=eta,xi=xi,**audit))
            for key in ['collision_symmetry','charge_residual','detailed_balance','drive_charge_overlap','weighted_transform_error']:
                check(rows,f'{key}_eta{eta}_xi{xi}',audit[key],1e-11)
            check(rows,f'minimum_eigenvalue_eta{eta}_xi{xi}',audit['minimum_eigenvalue'],-1e-11*audit['lambda_max'],lower=True,notes='Signed raw value; no clipping')
            check(rows,f'zero_mode_count_eta{eta}_xi{xi}',audit['zero_mode_count']-(2 if eta==0 else 1),0.)
    return output


def solver_checks(p,rows):
    data=[]
    dipoles=[]
    for eta in [0.,.1,10.]:
        for omega in [0.,.03]:
            pp=replace(p,eta=eta)
            out,raw=compute(pp,omega)
            check(rows,f'block_eta{eta}_omega{omega}',out['block_error'],1e-9)
            check(rows,f'full_projected_eta{eta}_omega{omega}',out['full_projected_error'],1e-9)
            check(rows,f'solve_eta{eta}_omega{omega}',out['solve_residual'],1e-10)
            st={k:raw[k] for k in contours(pp)}; _,_,L=collision(st,pp)
            gy=full_solution(L,np.sqrt(st['weight'])*st['vy'],omega,p.zero_rtol)
            forbidden=abs(raw['h']@gy)
            check(rows,f'forbidden_H_y_drive_eta{eta}_omega{omega}',forbidden,1e-12)
            if omega:
                # Direct full-space Schur test at a well-conditioned finite frequency.
                n=2*p.n; M=L.astype(complex)+1j*omega*np.eye(len(L))
                sol=solve(M[n:,n:],np.column_stack((M[n:,:n],raw['b'][n:])))
                ga=solve(M[:n,:n]-M[:n,n:]@sol[:,:n],raw['b'][:n]-M[:n,n:]@sol[:,n])
                check(rows,f'full_space_finite_frequency_Schur_eta{eta}',norm(ga-raw['g'][:n])/norm(ga),1e-10)
            dipoles.append(out['D_x'])
            data.append(dict(eta=eta,omega=omega,block_error=out['block_error'],full_projected_error=out['full_projected_error']))
    check(rows,'D_x_eta_invariance',max(dipoles)-min(dipoles),1e-14)
    return data
