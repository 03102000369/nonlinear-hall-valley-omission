import numpy as np
from nonlinear_hall.referee_numerics import PreparedSystem, RevisionParameters
from nonlinear_hall.mode_ablation import basis, hierarchy, spectrum, characterize
from nonlinear_hall.final_physics import physical_basis


def test_weighted_physical_modes_and_mirror():
    s=PreparedSystem(RevisionParameters(n=48,xi=.1))
    U,raw,T,scales=basis(s)
    assert np.linalg.norm(raw@T-U)<1e-13
    assert np.linalg.norm(U.T@U-np.eye(4))<1e-13
    assert np.linalg.norm(U-physical_basis(s)[0])<1e-13
    rows,_=characterize(s,{})
    assert all(abs(r['physical_mirror_expectation']+1)<1e-12 for r in rows)
    assert all(abs(r['local_qy_reflection_expectation']-1)<1e-12 for r in rows)
    assert np.linalg.norm(U[48:,:2])==0 and np.linalg.norm(U[:48,2:])==0


def test_three_contact_modes_close_at_dc_and_frequency():
    s=PreparedSystem(RevisionParameters(n=48,xi=0,delta_P=.4))
    for omega in [0.,.003]:
        for eta in [.001,1.,10.]:
            full,hs,r,_=hierarchy(s,eta,omega)
            assert r['rho_rel']<1e-12
            assert hs[2]['distribution_error']<1e-11
            assert hs[3]['error_H']<1e-11
            assert r['adjoint_identity_defect']<1e-11
            assert abs(full['chi_yxx']+full['H']/2)<1e-15


def test_finite_range_residual_and_spectral_identity():
    s=PreparedSystem(RevisionParameters(n=64,xi=.1,delta_P=.4))
    full,hs,r,raw=hierarchy(s,.001)
    assert r['rho_rel']>1e-4 and r['drive_projection_relative']>1e-4
    assert r['residual_decomposition_defect']<1e-12
    assert r['adjoint_identity_defect']<1e-11
    rows,summary=spectrum(s,.001,{})
    assert summary['spectral_identity_defect']<1e-11
    assert abs(sum(x['signed_error_contribution'] for x in rows)-(full['H']-hs[3]['H']))<1e-11
    assert abs(r['signed_H_error'])<=summary['resolvent_bound']
