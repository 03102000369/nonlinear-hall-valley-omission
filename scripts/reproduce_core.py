#!/usr/bin/env python3
"""Recompute the fixed physical hierarchy without modifying validated datasets."""
import json, time
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy.stats import spearmanr
from nonlinear_hall.referee_numerics import PreparedSystem, RevisionParameters
from nonlinear_hall.final_physics import cancellation
from nonlinear_hall.mode_ablation import basis, hierarchy, characterize, spectrum
from nonlinear_hall.io import read_csv, write_csv, write_json, sha
ROOT=Path(__file__).resolve().parents[1]
def main():
    start=time.perf_counter()
    from release import verify_data, compare_tables
    verify_data(ROOT)
    cfg=json.loads((ROOT/'config/study.json').read_text())
    def parameters(**kwargs):
        return RevisionParameters(**(cfg['model'] | {'n':cfg['n']} | kwargs))
    run=ROOT/'build/reproduced';proc=run/'processed';proc.mkdir(parents=True,exist_ok=True)
    save=lambda name,rows:write_csv(proc/name,rows)
    principal=[dict(kernel=cfg['principal']['kernel'],delta_P=d,xi=x,eta=e)
               for x in cfg['principal']['xi'] for d in cfg['principal']['delta_P'] for e in cfg['principal']['eta']]
    comparisons=[dict(kernel=k,delta_P=d,xi=x,eta=e)
                 for k in cfg['comparison']['kernels'] for d in cfg['comparison']['delta_P'] for x in cfg['comparison']['xi'] for e in cfg['comparison']['eta']]
    reps=[{k:v for k,v in r.items() if k!='kernel'} for r in cfg['representatives']]
    groups=defaultdict(list)
    for r in principal:groups[(r['kernel'],r['delta_P'],r['xi'])].append(dict(**r,sample='principal',label=''))
    for r in comparisons:groups[(r['kernel'],r['delta_P'],r['xi'])].append(dict(**r,sample='kernel_comparison',label=''))
    for r in reps:groups[('gaussian',r['delta_P'],r['xi'])].append(dict(**r,kernel='gaussian',sample='representative'))
    wide=[];long=[];residuals=[];modes=[];transforms=[];decisions=[];gates=[]
    for index,((kernel,delta,xi),points) in enumerate(groups.items()):
        system=PreparedSystem(parameters(delta_P=delta,xi=xi,kernel=kernel))
        prepared=basis(system);meta=dict(kernel=kernel,delta_P=delta,xi=xi)
        c,t=characterize(system,meta,prepared);modes.extend(c);transforms.extend(t)
        for source in points:
            eta=source['eta'];tag=dict(**meta,eta=eta,sample=source['sample'],label=source['label'])
            full,hs,defect,raw=hierarchy(system,eta,prepared=prepared)
            row=dict(**tag,H0=system.H0,**{'full_'+k:v for k,v in full.items()})
            for level in hs:
                m=int(level['dimension']);long.append(dict(**tag,**level))
                row.update({f'M{m}_'+k:v for k,v in level.items() if k!='dimension'})
            wide.append(row);residuals.append(dict(**tag,**defect))
            assert defect['orthogonality_error']<1e-12 and defect['transform_error']<1e-12
            assert defect['adjoint_identity_defect']<1e-10 and defect['residual_decomposition_defect']<1e-11
            assert abs(defect['signed_H_error'])<=defect['adjoint_bound']+1e-11
            if xi==0:
                assert defect['rho_rel']<1e-12 and hs[2]['distribution_error']<1e-10
            if source['sample']=='principal':
                # Active scalar/shape split evaluated independently from the solved response.
                g=raw['g'];g0=system.g0[:cfg['n']]
                alpha=(g0@g[:cfg['n']])/(g0@g0)
                S=(alpha-1)*system.H0;T=system.ha@(g[:cfg['n']]-alpha*g0)
                for tol in cfg['tolerances']:decisions.append(dict(**tag,**cancellation(system.H0,S,T,tol)))
            gates.append(dict(check='adjoint_projection_identity',**meta,eta=eta,value=defect['adjoint_identity_defect'],threshold=1e-10,passed=True))
            for key,threshold in [('orthogonality_error',1e-12),('transform_error',1e-12),('residual_decomposition_defect',1e-11)]:
                gates.append(dict(check=key,**meta,eta=eta,value=defect[key],threshold=threshold,passed=defect[key]<threshold))
            violation=max(0.,abs(defect['signed_H_error'])-defect['orthogonal_adjoint_bound'])
            assert violation<1e-11
            gates.append(dict(check='orthogonal_adjoint_bound_violation',**meta,eta=eta,value=violation,threshold=1e-11,passed=True))
            if xi==0:
                gates.append(dict(check='contact_invariant_residual',**meta,eta=eta,value=defect['rho_rel'],threshold=1e-12,passed=True))
                gates.append(dict(check='contact_M3_distribution_error',**meta,eta=eta,value=hs[2]['distribution_error'],threshold=1e-10,passed=True))
        if (index+1)%10==0: print('Completed physical configurations',index+1,'/',len(groups),flush=True)
    save('final_mode_ablation.csv',wide)
    save('final_invariant_subspace_residual.csv',residuals)
    save('final_mode_characterization.csv',modes);save('final_basis_transformation.csv',transforms)
    save('final_tolerance_cancellation.csv',decisions)
    statistics=[]
    for subset in ['positive_eta','all']:
        for tol in cfg['tolerances']:
            rr=[r for r in decisions if r['tolerance']==tol and (subset=='all' or r['eta']>0)]
            safe=sum(r['decision'] in ['ROBUSTLY_SAFE','CANCELLATION_DEPENDENT_SAFE'] for r in rr)
            for category in ['SAFE','ROBUSTLY_SAFE','CANCELLATION_DEPENDENT_SAFE','UNSAFE','UNDEFINED']:
                count=safe if category=='SAFE' else sum(r['decision']==category for r in rr)
                statistics.append(dict(subset=subset,tolerance=tol,category=category,numerator=count,denominator=len(rr),percentage=100*count/len(rr),
                   safe_denominator=safe,percentage_of_safe=100*count/safe if category in ['SAFE','ROBUSTLY_SAFE','CANCELLATION_DEPENDENT_SAFE'] else np.nan))
    save('final_cancellation_statistics.csv',statistics)
    # Deliberately small spectral diagnosis: two existing finite-range cases and matched kernel controls.
    spectral=[];summaries=[]
    for label,kernel,eta,delta,xi in cfg['spectral_cases']:
        system=PreparedSystem(parameters(delta_P=delta,xi=xi,kernel=kernel))
        meta=dict(label=label,kernel=kernel,eta=eta,delta_P=delta,xi=xi)
        rr,ss=spectrum(system,eta,meta);spectral.extend(rr);summaries.append(ss)
        assert ss['spectral_identity_defect']<1e-10
        gates.append(dict(check='spectral_error_identity',**meta,value=ss['spectral_identity_defect'],threshold=1e-10,passed=True))
    save('final_leakage_spectrum.csv',spectral);save('final_leakage_summary.csv',summaries)
    validation=[]
    for subset in ['principal_positive','kernel_unique_positive','representative']:
        rr=[r for r in long if (subset=='principal_positive' and r['sample']=='principal' and r['eta']>0) or
            (subset=='kernel_unique_positive' and r['sample']=='kernel_comparison' and r['eta']>0 and not(r['kernel']=='lorentzian' and r['xi']==0)) or
            (subset=='representative' and r['sample']=='representative')]
        for dim in range(1,5):
            rows=[r for r in rr if r['dimension']==dim]
            for metric in ['error_H','error_sigma_xx','error_Q','error_Q_A','distribution_error','H_error_over_H0']:
                values=[r[metric] for r in rows]
                validation.append(dict(subset=subset,dimension=dim,metric=metric,n=len(rows),mean=np.mean(values),p95=np.quantile(values,.95),maximum=max(values),
                      below_one_percent=sum(v<=.01 for v in values),below_five_percent=sum(v<=.05 for v in values)))
    save('final_reduced_model_validation.csv',validation)
    correlations=[]
    for xi in [None,.05,.1]:
        rr=[r for r in residuals if r['sample']=='principal' and r['eta']>0 and r['xi']>0 and (xi is None or r['xi']==xi)]
        for metric in ['H_error','distribution_error','Q_A_error']:
            correlations.append(dict(subset='principal_positive_finite_range',xi='all_finite' if xi is None else xi,metric=metric,n=len(rr),spearman_rho=spearmanr([r['rho_rel'] for r in rr],[r[metric] for r in rr]).statistic))
    save('final_leakage_error_correlations.csv',correlations)
    agreements=[]
    for m in range(1,5):
        for tol in cfg['tolerances']:
            rr=[r for r in wide if r['sample']=='principal' and r['eta']>0]
            fullsafe=np.array([abs(r['full_H']-r['H0'])/abs(r['full_H'])<=tol for r in rr])
            predsafe=np.array([abs(r[f'M{m}_H']-r['H0'])/abs(r[f'M{m}_H'])<=tol for r in rr])
            agreements.append(dict(dimension=m,tolerance=tol,n=len(rr),agree=int(sum(fullsafe==predsafe)),false_safe=int(sum(predsafe&~fullsafe)),false_unsafe=int(sum(~predsafe&fullsafe))))
    save('final_mode_decision_agreement.csv',agreements)
    # Recompute the retained mass and physical amplitude controls from explicit settings.
    from nonlinear_hall.final_physics import passive_population, full_response, reduced_response
    mass=[]
    for case in cfg['mass_cases']:
        for m in cfg['masses']:
            system=PreparedSystem(parameters(delta_P=case['delta_P'],xi=case['xi'],m_P=m))
            HP0=float(system.h[cfg['n']:]@system.g0[cfg['n']:])
            for eta in [0,case['eta']]:
                r,raw=system.evaluate(eta,details=True);Np,Nm=passive_population(system,raw['g']);HP=r['H_passive_direct']
                direct=HP0;ac=r['H_active']-r['H0'];pc=HP-HP0
                via=-m*system.p.v_P**2*(Np-Nm)/(2*system.p.delta_P**3)
                row=dict(label=case['label'],eta=eta,delta_P=case['delta_P'],xi=case['xi'],m_P=m,H_A=r['H_active'],H_P=HP,H_full=r['H'],H0=r['H0'],
                    DIRECT_PASSIVE=direct,ACTIVE_COUPLING=ac,PASSIVE_COUPLING=pc,DeltaH=r['H']-r['H0'],relative_omission=r['actual_error'],
                    H_P_over_m=HP/m if m else np.nan,N_plus=Np,N_minus=Nm,valley_polarization=Np-Nm,polarization_prediction=via,
                    decomposition_residual=abs(r['H']-r['H0']-direct-ac-pc),polarization_residual=abs(HP-via),
                    valid=r['valid_physical_numerical'],minimum_driven_eigenvalue=r['minimum_driven_eigenvalue'],solve_residual=r['solve_residual'])
                assert row['decomposition_residual']<1e-11 and row['polarization_residual']<1e-11 and row['valid']
                mass.append(row)
    save('final_mass_scaling.csv',mass)
    coefficients=[]
    for c in cfg['representatives']:
        if c['label']=='matched':continue
        params={k:c[k] for k in ['eta','delta_P','xi']}
        system=PreparedSystem(parameters(**params))
        full,_,_=full_response(system,c['eta']);red,_=reduced_response(system,c['eta'],4)
        coefficients.append(dict(label=c['label'],**params,H0=system.H0,**red,full_H=full['H'],full_S=full['S'],full_T=full['T'],full_shape_fraction=full['shape_fraction']))
    save('final_reduced_mode_coefficients.csv',coefficients)
    # Only aggregate diagnostics are retained, never the large per-point ledger.
    checks=[]
    for name in sorted(set(r['check'] for r in gates)):
        rr=[r for r in gates if r['check']==name]
        checks.append(dict(check=name,count=len(rr),maximum=max(r['value'] for r in rr),threshold=rr[0]['threshold'],passed=all(r['passed'] for r in rr)))
    write_csv(run/'validation_summary.csv',checks)
    comparison=compare_tables(ROOT,proc)
    write_json(run/'report.json',dict(status='PASS',canonical_run='20260919T021436Z-40be24e6',
               configuration_sha256=sha(ROOT/'config/study.json'),wall_seconds=time.perf_counter()-start,
               principal_points=len(principal),positive_points=sum(r['eta']>0 for r in principal),
               invariant_checks=len(gates),comparison=comparison))
    print('Core reconstruction and canonical numerical comparison PASS',len(principal),'principal points',flush=True)

if __name__=='__main__':
    main()
