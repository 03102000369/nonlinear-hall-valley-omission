import numpy as np
from nonlinear_hall.referee_numerics import PreparedSystem,RevisionParameters
from nonlinear_hall.final_physics import cancellation,full_response,reduced_response,passive_population


def test_contact_physical_subspace_is_closed():
    system=PreparedSystem(RevisionParameters(n=48,xi=0,delta_P=.4))
    for eta in [.001,.372759372,10]:
        exact,_,_=full_response(system,eta);red,_=reduced_response(system,eta)
        assert red['closure_residual']<1e-12
        for key in ['H','S','T']:
            assert abs(red[key]-exact[key])<1e-11
        assert abs(red['x1']-red['shape_forcing']/red['shape_stiffness'])<1e-10


def test_passive_population_readout_and_mass_baseline():
    system=PreparedSystem(RevisionParameters(n=64,xi=.1,delta_P=.4,m_P=.01))
    HP0=system.h[64:]@system.g0[64:]
    for eta in [0,.001]:
        row,raw=system.evaluate(eta,details=True);Np,Nm=passive_population(system,raw['g'])
        HP=-.01*(Np-Nm)/(2*.4**3)
        assert abs(HP-row['H_passive_direct'])<1e-12
        assert abs(row['H']-row['H0']-(HP0+row['H_active']-row['H0']+HP-HP0))<1e-12
    assert abs(HP0)>2.9


def test_cancellation_envelope_and_relative_denominator():
    row=cancellation(-1,.02,-.019,.01)
    assert row['decision']=='CANCELLATION_DEPENDENT_SAFE'
    assert abs(row['no_cancellation_envelope']-.039/.961)<1e-15
    assert row['discrete_sign_worst']==row['no_cancellation_envelope']
    assert cancellation(0,.1,.1,.01)['decision']=='UNDEFINED'
    assert cancellation(-1,.002,-.001,.01)['decision']=='ROBUSTLY_SAFE'
    assert np.isinf(cancellation(-1,2,-1.999,.01)['no_cancellation_envelope'])
    assert np.isfinite(cancellation(-1,2,-1.999,.01)['discrete_sign_worst'])
