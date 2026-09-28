"""Selected technical controls from included processed tables."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .io import read_csv,write_json,sha
COLORS=['#126b88','#cb7436','#7c5aa6','#2e8b68']
LABELS={'maximum-omission':'Scalar dominated','10pct-crossover':'Shape dominated','weak-cancellation':'Weak coupling'}

def make_supporting_figures(root):
    root=Path(root);base=data=root/'data/processed';dest=root/'figures/supplement';dest.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'legend.fontsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42,'savefig.dpi':220})
    manifest=[]
    def rows(name):return read_csv(base/name)
    load=rows
    def save(fig,name,inputs):
        fig.savefig(dest/(name+'.pdf'),bbox_inches='tight',metadata={'CreationDate':None});fig.savefig(dest/(name+'.png'),bbox_inches='tight');plt.close(fig)
        manifest.append(dict(figure=name,inputs={str(p.relative_to(root)):sha(p) for p in inputs}))
    finish=save
    weak=load('weak_coupling_remainder.csv');fig,axes=plt.subplots(1,3,figsize=(7.2,3),layout='constrained')
    for i,(key,title) in enumerate([('eta',r'$\eta$'),('rho_norm',r'$\|\eta R_0K\|_2$'),('rho_drive',r'$\eta\|R_0Kg_0\|/\|g_0\|$')]):
        for kernel,col,mark in [('gaussian',COLORS[0],'o'),('lorentzian',COLORS[1],'x')]:
            r=[r for r in weak if r['kernel']==kernel and r['eta']>0 and r['first_order_error']>=1e-10];axes[i].loglog([x[key] for x in r],[x['first_order_error'] for x in r],mark,color=col,ms=2,alpha=.3,label=kernel.capitalize())
        axes[i].set(xlabel=title,ylabel='First-order relative error',title=f'({chr(97+i)}) Perturbative organization');axes[i].axhline(.01,c='.4',ls=':');axes[i].legend(markerscale=2,fontsize=6)
    finish(fig,'figure_S2_weak_control',[data/'weak_coupling_remainder.csv'])

    conv=load('successive_resolution.csv');fig,ax=plt.subplots(figsize=(6,3.1),layout='constrained')
    for label,col,mark in zip(list(LABELS)[:3],COLORS,['o','s','^']):
        ss=[r for r in conv if r['label']==label];ns=sorted(set(r['n_to'] for r in ss));vals=[max(r['error'] for r in ss if r['n_to']==n and r['metric']=='relative') for n in ns];ax.semilogy(ns,np.maximum(vals,1e-16),marker=mark,color=col,label=LABELS[label])
    ax.axhline(.01,c='.3',ls=':');ax.set(xlabel='Successive angular-grid pair',ylabel='Maximum successive relative change',xticks=[64,96,128,192,256],xticklabels=['48→64','64→96','96→128','128→192','192→256']);ax.legend();ax.grid(alpha=.12)
    finish(fig,'figure_S4_successive_convergence',[data/'successive_resolution.csv'])

    mass=[r for r in rows('final_mass_scaling.csv') if r['label']=='weak'];fig,axes=plt.subplots(2,2,figsize=(7.1,5.15),layout='constrained')
    for eta,col in [(0,COLORS[0]),(.001,COLORS[1])]:
        rr=[r for r in mass if r['eta']==eta and r['m_P']>0];x=[r['m_P'] for r in rr]
        for ax,key,title in zip(axes.flat,['H_P','H_P_over_m','valley_polarization'],['(a) Passive Hall response','(b) Finite small-mass susceptibility','(c) Passive valley displacement']):
            ax.plot(x,[r[key] for r in rr],'-o',ms=3,color=col,label=rf'$\eta={eta:g}$');ax.set_xscale('log');ax.set_xlabel(r'$m_P/m_A$');ax.set_title(title,loc='left',fontsize=9)
    for ax,yl in zip([axes[0,0],axes[0,1],axes[1,0]],[r'$H_P$',r'$H_P/m_P$',r'$N_+-N_-$']):ax.set_ylabel(yl)
    axes[0,0].legend(frameon=False)
    r=next(r for r in mass if r['eta']==.001 and r['m_P']==.01);ax=axes[1,1];keys=['DIRECT_PASSIVE','ACTIVE_COUPLING','PASSIVE_COUPLING'];values=[r[k] for k in keys]
    ax.bar(range(3),values,color=COLORS[:3]);ax.axhline(0,color='#555555',lw=.7);ax.set_xticks(range(3),['Direct\npassive','Active\ncoupling','Passive\ncoupling']);ax.set_title('(d) Exact correction at $m_P=0.01$',loc='left',fontsize=9);ax.set_ylabel(r'Contribution to $H-H_0$');ax.set_ylim(-3.5,.5)
    for i,v in enumerate(values):ax.text(i,v-(.12 if v<0 else -.1),f'{v:.3g}',ha='center',va='top' if v<0 else 'bottom',fontsize=8)
    save(fig,'figure_S3_mass',[base/'final_mass_scaling.csv'])
    write_json(root/'figures/supplement_input_manifest.json',{'figures':manifest,'generator':sha(Path(__file__))})
