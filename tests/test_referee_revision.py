import numpy as np
import pytest
from dataclasses import replace
from threadpoolctl import threadpool_limits
from nonlinear_hall.referee_numerics import RevisionParameters,PreparedSystem,contour_states,transition_parts,weak_series
from nonlinear_hall import core
from nonlinear_hall.manuscript_numerics import ManuscriptParameters

@pytest.fixture(autouse=True)
def one_thread():
    with threadpool_limits(limits=1):yield

@pytest.mark.parametrize('eta,delta,xi',[(10,1.2,0),(.372759372,.4,0),(.001,.4,.1),(0,.8,.05)])
def test_new_solver_matches_independent_legacy(eta,delta,xi):
    p=RevisionParameters(n=32,eta=eta,delta_P=delta,xi=xi)
    new,raw=PreparedSystem(p).evaluate(details=True)
    old,_=core.compute(ManuscriptParameters(n=32,A0=.01,eta=eta,delta_P=delta,xi=xi))
    assert new['H']==pytest.approx(old['H'],rel=1e-11,abs=1e-12)
    assert new['H0']==pytest.approx(old['H0'],rel=1e-11)
    assert new['Delta_H_total']==pytest.approx(new['Delta_H_scalar']+new['Delta_H_shape'],abs=1e-12)
    assert abs(np.vdot(raw['g0'],raw['residual']))<1e-11*np.linalg.norm(raw['g0'])*np.linalg.norm(raw['g'][:32])
    assert new['adjoint_identity']<1e-11
    assert new['collision_identity']<1e-11


def test_gauge_covariance_and_physical_rotation_are_distinct():
    p=RevisionParameters(n=32,m_P=.02)
    base=PreparedSystem(p);rb=base.evaluate()
    gauge=PreparedSystem(replace(p,rotation=np.pi/3,internal_model='gauge'));rg=gauge.evaluate()
    for a,b in zip(transition_parts(base.st,p),transition_parts(gauge.st,gauge.p)):
        assert np.allclose(a,b,rtol=1e-12,atol=1e-15)
    for key in ['H','sigma_xx','sigma_yy','fixed_E_signal','fixed_j_signal']:
        assert rb[key]==pytest.approx(rg[key],rel=1e-11)
    physical=PreparedSystem(replace(p,rotation=np.pi/3)).evaluate()
    assert abs(physical['H']-rb['H'])>.01*abs(rb['H'])
    assert physical['sigma_xy_residual']<1e-12


def test_origin_enclosure_stronger_than_type_I():
    with pytest.raises(ValueError,match='Origin-enclosing'):
        contour_states(RevisionParameters(n=16,m_P=.5,delta_P=.4))
    with pytest.raises(ValueError,match='Origin-enclosing'):
        contour_states(RevisionParameters(n=16,m_P=-.5,delta_P=.4))


def test_remainder_bounds_dc_and_complex():
    s=PreparedSystem(RevisionParameters(n=32,xi=.1,delta_P=.4))
    for omega in [0,.01*s.gap]:
        for r in weak_series(s,[0,1e-4,.01,.1,1],omega):
            assert r['neumann_violation']<1e-11
            if omega==0:assert r['energy_violation']<1e-11


def test_second_harmonic_sign_and_complex_amplitude_from_time_domain():
    s=PreparedSystem(RevisionParameters(n=32));omega=.003
    row,raw=s.evaluate(omega=omega,details=True)
    # Reconstruct physical distribution and anomalous velocity independently of chi.
    g=core.lift(raw['g'],s.st,32);ell=g/np.sqrt(s.st['weight'])
    phase=2*np.pi*np.arange(4096)/4096;E0=.7+.2j
    field=np.real(E0*np.exp(1j*phase))
    delta_f=-np.real(E0*ell[:,None]*np.exp(1j*phase)[None,:])
    dotk=np.column_stack([-field,np.zeros(len(phase)),np.zeros(len(phase))])
    om=np.column_stack([np.zeros(len(g)),np.zeros(len(g)),s.st['Omega']])
    anomalous=-np.cross(dotk[None,:,:],om[:,None,:])[:,:,1]
    current=-np.sum(s.st['weight'][:,None]*anomalous*delta_f,axis=0)
    harmonic=2*np.mean(current*np.exp(-2j*phase))
    assert harmonic==pytest.approx(row['chi_yxx']*E0**2,rel=1e-12,abs=1e-12)
    # Complex current-driven normalization uses sigma^2, not |sigma|^2.
    assert row['fixed_j_signal']==pytest.approx(-row['chi_yxx']/(row['sigma_yy']*row['sigma_xx']**2))


def test_decoupled_H_does_not_imply_decoupled_voltage():
    r=PreparedSystem(RevisionParameters(n=32,eta=0)).evaluate()
    assert r['H_omission']<1e-12
    assert r['fixed_E_omission']>.1
    assert r['fixed_j_omission']>r['fixed_E_omission']
