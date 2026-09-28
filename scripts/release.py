#!/usr/bin/env python3
"""Fail-closed checks, canonical comparisons and portable manuscript builds."""
import argparse, json, re, shutil, subprocess, sys, math
from pathlib import Path
import numpy as np
from nonlinear_hall.io import sha, read_csv, write_json, write_csv
ROOT=Path(__file__).resolve().parents[1]

def verify_data(root=ROOT):
    inventory=json.loads((root/'data/manifest.json').read_text())
    expected=set(inventory)
    actual={str(p.relative_to(root)) for p in (root/'data/processed').glob('*.csv')}
    if actual!=expected:raise RuntimeError('Dataset inventory differs from files')
    for path,record in inventory.items():
        if sha(root/path)!=record['sha256']:raise RuntimeError('Dataset checksum mismatch: '+path)
    return inventory

# Key fields sort independently of scan traversal. Floating-point comparison tolerances
# accommodate LAPACK/BLAS changes, but categorical decisions must match exactly.
KEYS={
 'final_mode_ablation.csv':['sample','label','kernel','delta_P','xi','eta'],
 'final_invariant_subspace_residual.csv':['sample','label','kernel','delta_P','xi','eta'],
 'final_tolerance_cancellation.csv':['sample','label','kernel','delta_P','xi','eta','tolerance'],
 'final_basis_transformation.csv':['kernel','delta_P','xi','raw_mode','orthonormal_mode'],
 'final_mode_characterization.csv':['kernel','delta_P','xi','mode'],
 'final_mass_scaling.csv':['label','m_P','eta'],
 'final_reduced_mode_coefficients.csv':['label'],
 'final_cancellation_statistics.csv':['subset','tolerance','category'],
 'final_reduced_model_validation.csv':['subset','dimension','metric'],
 'final_leakage_error_correlations.csv':['subset','xi','metric'],
 'final_mode_decision_agreement.csv':['dimension','tolerance'],
 'final_leakage_summary.csv':['label'],
 'final_leakage_spectrum.csv':['label','eigenmode']}

def compare_tables(root, produced):
    report=[]
    for name in KEYS:
        reference=read_csv(root/'data/processed'/name);new=read_csv(produced/name)
        if len(reference)!=len(new):raise RuntimeError('Row count mismatch: '+name)
        keys=[k for k in KEYS[name] if k in reference[0]]
        # Numeric key ordering is lexical but identical on both sides, including exact grids.
        key=lambda r:tuple(str(r[k]) for k in keys)
        a=sorted(reference,key=key);b=sorted(new,key=key)
        maxdiff=0.;compared=0
        for old,now in zip(a,b):
            if old.keys()!=now.keys():raise RuntimeError('Schema mismatch: '+name)
            for col,val in old.items():
                v=now[col]
                if isinstance(val,float) and isinstance(v,float):
                    # Individual eigenvector signs and degenerate mixing are not physical.
                    # The canonical environment is expected to match here; grouped invariants
                    # and all signed sums are additionally checked by the reconstruction.
                    if not np.isclose(val,v,atol=1e-9,rtol=1e-7,equal_nan=True):
                        raise RuntimeError(f'Canonical mismatch {name} {key(old)} {col}: {val} vs {v}')
                    if np.isfinite(val) and np.isfinite(v):maxdiff=max(maxdiff,abs(val-v))
                elif val!=v:raise RuntimeError(f'Categorical mismatch {name} {key(old)} {col}')
                compared+=1
        report.append(dict(dataset=name,rows=len(a),cells=compared,max_absolute_difference=maxdiff,byte_identical=sha(root/'data/processed'/name)==sha(produced/name),status='PASS'))
    return report

def validate():
    inventory=verify_data()
    decisions=read_csv(ROOT/'data/processed/final_tolerance_cancellation.csv')
    for r in decisions:
        if r['decision']=='CANCELLATION_DEPENDENT_SAFE':assert r['B']<1 and r['opposite_signs']=='True'
    from nonlinear_hall.core import Parameters, contours
    from nonlinear_hall.validation import band_checks,tau_checks,solver_checks
    rows=[];p=Parameters(n=16)
    band_checks(p,rows);tau_checks(contours(p),rows);solver_checks(p,rows)
    from nonlinear_hall.referee_numerics import RevisionParameters,PreparedSystem
    from nonlinear_hall.mode_ablation import hierarchy
    r=next(r for r in read_csv(ROOT/'data/processed/final_mode_ablation.csv') if r['sample']=='representative' and r['label']=='decision_cancellation')
    sysm=PreparedSystem(RevisionParameters(n=256,eta=r['eta'],delta_P=r['delta_P'],xi=r['xi']))
    full,hs,defect,_=hierarchy(sysm,r['eta'])
    assert np.isclose(full['H'],r['full_H'],rtol=1e-10,atol=1e-11)
    assert defect['adjoint_identity_defect']<1e-10
    for m,h in enumerate(hs,1):assert np.isclose(h['H'],r[f'M{m}_H'],rtol=1e-10,atol=1e-11)
    write_csv(ROOT/'build/fast_validation.csv',rows)
    print(f'PASS: {len(inventory)} dataset hashes, {len(rows)} fresh independent checks, N=256 representative M1-M4 regression')

def build():
    verify_data()
    env=None
    for folder,stem,bib in [('manuscript','main',True),('supplement','supplement',False)]:
        wd=ROOT/folder
        for source in wd.rglob('*.tex'):
            for name in re.findall(r'\\(?:input|includegraphics)(?:\[[^]]*\])?\{([^}]+)\}',source.read_text()):
                target=wd/name
                if not target.suffix:target=target.with_suffix('.tex')
                if not target.is_file():raise RuntimeError('Missing TeX input '+name)
        cmds=[['pdflatex','-interaction=nonstopmode','-halt-on-error',stem+'.tex']]
        if bib:cmds.append(['bibtex',stem])
        cmds += [['pdflatex','-interaction=nonstopmode','-halt-on-error',stem+'.tex']]*2
        for cmd in cmds:
            p=subprocess.run(cmd,cwd=wd,capture_output=True,text=True)
            if p.returncode:raise RuntimeError(p.stdout[-5000:]+p.stderr[-1000:])
        log=(wd/(stem+'.log')).read_text()
        bad=re.findall(r'^.*(?:Overfull \\[hv]box|undefined references|Citation .* undefined|Reference .* undefined).*$' ,log,re.M)
        if bad:raise RuntimeError('\n'.join(bad))
        print(folder, re.findall(r'Output written on .*',log)[-1])

def figures():
    verify_data()
    from nonlinear_hall.mode_ablation_figures import make_figures
    from nonlinear_hall.supporting_figures import make_supporting_figures
    make_figures(ROOT);make_supporting_figures(ROOT)
    from paper_tables import tables
    tables()
    subprocess.run([sys.executable,str(ROOT/'scripts/claims.py')],check=True,capture_output=True)
    inputs={
        'main/figure_1_physical_modes':['final_mode_characterization.csv','final_basis_transformation.csv'],
        'main/figure_2_ablation':['final_mode_ablation.csv'],
        'main/figure_3_leakage':['final_invariant_subspace_residual.csv','final_leakage_spectrum.csv'],
        'main/figure_4_scalar_shape':['final_reduced_mode_coefficients.csv'],
        'main/figure_5_decision_map':['final_tolerance_cancellation.csv','final_cancellation_statistics.csv'],
        'supplement/figure_S1_spectral_weights':['final_leakage_spectrum.csv'],
        'supplement/figure_S2_weak_control':['weak_coupling_remainder.csv'],
        'supplement/figure_S3_mass':['final_mass_scaling.csv'],
        'supplement/figure_S4_successive_convergence':['successive_resolution.csv']}
    records={}
    for name,files in inputs.items():
        generator='supporting_figures.py' if name.startswith('supplement/') and 'S1_' not in name else 'mode_ablation_figures.py'
        records[name]=dict(inputs={f'data/processed/{f}':sha(ROOT/'data/processed'/f) for f in files},
            configuration_sha256=sha(ROOT/'config/study.json'),generator=generator,
            generator_sha256=sha(ROOT/'src/nonlinear_hall'/generator),
            outputs={ext:sha(ROOT/'figures'/(name+'.'+ext)) for ext in ['pdf','png']})
    write_json(ROOT/'figures/manifest.json',records)
    for p in ['main_input_manifest.json','supplement_input_manifest.json']:
        (ROOT/'figures'/p).unlink(missing_ok=True)
    print('Regenerated all five main and four supplementary figures and supplemental tables')

def clean():
    for folder in ['build','.pytest_cache']:
        if (ROOT/folder).exists():shutil.rmtree(ROOT/folder)
    for folder in ['src','tests','scripts']:
        for p in (ROOT/folder).rglob('__pycache__'):shutil.rmtree(p)
    for folder in ['manuscript','supplement']:
        for p in (ROOT/folder).rglob('*'):
            if p.is_file() and (p.suffix in ['.aux','.log','.out','.toc','.fls','.bbl','.blg','.fdb_latexmk'] or p.name.endswith('.synctex.gz')):p.unlink()

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['validate','figures','build','clean']);args=ap.parse_args()
    globals()[args.action]()
