"""Generate supplemental tables from canonical data."""
from pathlib import Path
from nonlinear_hall.io import read_csv
ROOT=Path(__file__).resolve().parents[1]
def tables():
    ab=read_csv(ROOT/'data/processed/final_mode_ablation.csv')
    inv=read_csv(ROOT/'data/processed/final_invariant_subspace_residual.csv')
    tex=[r'\begin{table}[ht]\centering\small',r'\begin{tabular}{llrrrr}\toprule',r'Case & Model & Hall error (\%) & $\sigma_{xx}$ error (\%) & $Q_A$ error (\%) & Distribution (\%)\\\midrule']
    for r in ab:
        if r['sample']!='representative':continue
        name={'decision_cancellation':'Fixed-budget','matched':'Matched range'}.get(r['label'],r['label'].capitalize())
        for m in range(1,5):tex.append(name+f' & M{m} & '+' & '.join(f"{100*r[f'M{m}_'+k]:.4g}" for k in ['error_H','error_sigma_xx','error_Q_A','distribution_error'])+r'\\')
    tex += [r'\bottomrule\end{tabular}',r'\caption{Extended ordered ablation at the five fixed configurations. All errors use the full observable or weighted distribution norm as denominator. Tiny contact errors are roundoff. No monotonic Hall-error ordering is imposed.}\end{table}']
    tex += [r'\begin{table}[ht]\centering\small',r'\begin{tabular}{lrrrr}\toprule',r'Kernel / range & Min. $\rho_{\rm rel}$ & Max. $\rho_{\rm rel}$ & Max. $\rho_{\rm abs}$ & Max. Hall error (\%)\\\midrule']
    for kernel in ['gaussian','lorentzian']:
        for xi in [0,.05,.1]:
            rr=[r for r in inv if r['sample']=='kernel_comparison' and r['eta']>0 and r['kernel']==kernel and r['xi']==xi]
            tex.append(f'{kernel.capitalize()} / {xi:.2f} & '+' & '.join(f'{v:.4g}' for v in [min(r['rho_rel'] for r in rr),max(r['rho_rel'] for r in rr),max(r['rho_abs'] for r in rr),100*max(r['H_error'] for r in rr)])+r'\\')
    tex += [r'\bottomrule\end{tabular}',r'\caption{Invariant residuals and M4 accuracy on the two-kernel grid, restricted to positive coupling. Contact entries coincide because the kernels are identical. Absolute residuals have collision-rate units; relative residuals are dimensionless.}\end{table}']
    tex += [r'\begin{table}[ht]\centering\small',r'\begin{tabular}{rrrrrr}\toprule',r'Tolerance & Safe & Robust & Dependent & Unsafe & Dependent/safe (\%)\\\midrule']
    counts=read_csv(ROOT/'data/processed/final_cancellation_statistics.csv')
    for tol in [.01,.05,.1]:
        rr={r['category']:r for r in counts if r['subset']=='all' and r['tolerance']==tol}
        tex.append(f'{100*tol:.0f}'+r'\% & '+' & '.join(str(int(rr[k]['numerator'])) for k in ['SAFE','ROBUSTLY_SAFE','CANCELLATION_DEPENDENT_SAFE','UNSAFE'])+f" & {rr['CANCELLATION_DEPENDENT_SAFE']['percentage_of_safe']:.2f}"+r'\\')
    tex += [r'\bottomrule\end{tabular}',r'\caption{Supporting all-grid counts, including 90 decoupled controls: denominator 3600; undefined count zero. Dependent means cancellation-dependent under the fixed-budget test. The principal statistics in the main text exclude these controls.}\end{table}']
    d=ROOT/'supplement/generated';d.mkdir(exist_ok=True);(d/'mode_ablation_tables.tex').write_text('\n'.join(tex)+'\n')
