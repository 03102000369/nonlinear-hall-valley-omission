from dataclasses import replace
import numpy as np
import pytest
from scipy.linalg import norm,solve
from nonlinear_hall import core
from nonlinear_hall.manuscript_numerics import ManuscriptParameters,independent_system,independent_observation,perturbation_system,shape_diagnostic,validity,common_ratio

@pytest.mark.parametrize('kernel',['gaussian','lorentzian'])
def test_independent_root_spinor_operator_and_response(kernel):
    p=ManuscriptParameters(n=24,kernel=kernel,A0=.01,eta=.3)
    st=core.contours(p);W,C,L=core.collision(st,p)
    si,wi,ci,li,*_=independent_system(p)
    assert np.max(abs(st['weight']-si['weight']))<1e-14
    assert norm(L-li)/norm(L)<1e-13
    obs,raw=core.compute(p);ind,_=independent_observation(p)
    for k in ['D_x','H','sigma_xx','Q','Q_A','delta_H']:
        assert abs(obs[k]-ind[k])/abs(obs[k])<1e-10

@pytest.mark.parametrize('kernel',['gaussian','lorentzian'])
@pytest.mark.parametrize('omega',[0.,.008])
def test_cross_edge_derivative_and_first_order_remainder(kernel,omega):
    p=ManuscriptParameters(n=24,kernel=kernel,A0=.01)
    s=perturbation_system(p,omega)
    _,_,L=core.collision(s['st'],replace(p,eta=.13))
    assert norm(core.odd_matrix(L,s['st'],p.n)-s['Lo']-.13*s['Ko'])<1e-14
    assert abs(s['H1']-s['active_out']-s['passive_source'])<1e-12
    errors=[]
    for eta in [1e-3,5e-4]:
        g=solve(s['Lo']+eta*s['Ko']+1j*omega*np.eye(2*p.n),s['bo'])
        errors.append(abs(s['ho']@g-s['H0']-eta*s['H1']))
    assert 3.8<errors[0]/errors[1]<4.2


def test_shape_and_mask_have_nontrivial_guards():
    p=ManuscriptParameters(n=24,A0=.01)
    o,r=core.compute(p);shape,alpha,g0,diff=shape_diagnostic(r,p.n)
    assert shape>1e-3
    assert abs(np.vdot(g0,diff))<1e-10*norm(g0)*norm(diff)
    assert validity(p,o,r,True)['valid_for_main_text']
    assert not validity(p,o,r,False)['valid_for_main_text']
    o['lambda_max']=o['interband_threshold']
    assert not validity(p,o,r,True)['weak_scattering']


def test_kernel_contact_identity_and_common_ratio():
    p=ManuscriptParameters(n=24,xi=0)
    st=core.contours(p)
    assert np.array_equal(core.collision(st,p)[2],core.collision(st,replace(p,kernel='lorentzian'))[2])
    groups=[dict(R_L=1.,R_H=2.,significant=True),dict(R_L=3.,R_H=6.,significant=True)]
    r,defect=common_ratio(groups)
    assert r==2 and defect==0
