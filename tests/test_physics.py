from dataclasses import replace,asdict
import numpy as np
import pytest
from scipy.linalg import eigvalsh, norm
from nonlinear_hall.core import *
from nonlinear_hall.validation import band_checks,tau_checks,solver_checks
from nonlinear_hall.io import write_json,save_npz

@pytest.fixture
def p(): return Parameters(n=16)

def test_analytical_eigenvalues_and_all_band_symmetries(p):
    rows=[]; band_checks(replace(p,n=64),rows)
    assert all(r['status']=='PASS' for r in rows)

@pytest.mark.parametrize('r',['A','P'])
def test_time_reversal_and_berry_transform(p,r):
    q=np.array([.7,-.3]);T=1j*SY
    assert np.allclose(T@hamiltonian(q,r,1,p).conj()@T.conj().T,hamiltonian(-q,r,-1,p),atol=1e-14)
    assert quantities(q,r,1,p)[3]==-quantities(-q,r,-1,p)[3]

def test_passive_curvature_and_positive_weights(p):
    st=contours(p)
    assert np.min(st['weight'])>0
    assert np.all(st['Omega'][st['sector']=='P']==0)

def test_circular_quadrature(p):
    p=replace(p,t_A=0);st=contours(p)
    assert abs(sum(st['weight'][:2*p.n])-p.mu/np.pi)<1e-14

def test_constant_tau_and_Q(p):
    rows=[];tau_checks(contours(p),rows)

@pytest.mark.parametrize('eta,expected',[(0.,2),(.001,1),(1.,1),(10.,1)])
def test_scattering_collision_positivity_symmetry_conservation(p,eta,expected):
    p=replace(p,eta=eta);st=contours(p);W,C,L=collision(st,p)
    vals,ev=audit_collision(st,p,W,C,L)
    assert W.min()>-1e-14
    assert norm(W-W.T)<1e-12
    assert vals['zero_mode_count']==expected
    assert ev[0]>-1e-12
    assert norm(L@np.sqrt(st['weight']))<1e-12

def test_block_equals_full_finite_frequency_and_DC(p):
    solver_checks(p,[])

def test_dipole_eta_invariant(p):
    ds=[vectors(contours(replace(p,eta=x))) for x in [0.,.1,10.]]
    assert all(a@b==ds[0][0]@ds[0][1] for a,b in ds)

def test_metadata_and_raw_exact_reload(tmp_path,p):
    write_json(tmp_path/'metadata.json',asdict(p))
    save_npz(tmp_path/'raw.npz',**contours(p))
    with pytest.raises(FileExistsError):save_npz(tmp_path/'raw.npz',a=np.array([1]))

def test_reject_absent_or_out_of_domain_contour(p):
    for pp in [replace(p,delta_P=0.),replace(p,mu=.5),replace(p,t_A=1.),replace(p,q_cutoff=.1)]:
        with pytest.raises(ValueError):contours(pp)

def test_sufficient_error_bound(p):
    for eta in [.001,.1,10.]:
        out,raw=compute(replace(p,eta=eta))
        assert out['Delta_H_abs']<=out['bound_abs']+1e-12

def test_A0_scaling(p):
    a,_=compute(p);b,_=compute(replace(p,A0=.01))
    for key in ['Q','Q_A','delta_H']:
        assert abs(a[key]-b[key])<1e-10
    assert abs(a['H']-b['H']*.01)<1e-10
