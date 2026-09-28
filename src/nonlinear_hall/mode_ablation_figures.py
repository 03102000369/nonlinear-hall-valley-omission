"""Manuscript figures generated exclusively from the validated final tables."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyArrowPatch
from matplotlib.colors import ListedColormap
from .io import read_csv,write_json,sha

BLUE='#256b9b';ORANGE='#cf642b';GREEN='#287c66';GRAY='#aeb5bd';PURPLE='#77529c'

def make_figures(root):
    root=Path(root);data=root/'data/processed';out=root/'figures/main';out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':9,'axes.labelsize':9,'axes.titlesize':10,'legend.fontsize':8,'font.family':'DejaVu Sans','pdf.fonttype':42,'savefig.bbox':'tight'})
    used={}
    def load(name):used[name]=sha(data/name);return read_csv(data/name)
    def save(fig,name):
        target=out
        if name.startswith('figure_S'):
            target=root/'figures/supplement';target.mkdir(parents=True,exist_ok=True)
            name='figure_S1_spectral_weights'
        fig.savefig(target/(name+'.pdf'),metadata={'CreationDate':None});fig.savefig(target/(name+'.png'),dpi=200);plt.close(fig)
    ablation=load('final_mode_ablation.csv');inv=load('final_invariant_subspace_residual.csv');spec=load('final_leakage_spectrum.csv')
    decisions=load('final_tolerance_cancellation.csv');counts=load('final_cancellation_statistics.csv')
    coeff=load('final_reduced_mode_coefficients.csv');reps={r['label']:r for r in ablation if r['sample']=='representative'}
    fig=plt.figure(figsize=(7.2,4.1));gs=fig.add_gridspec(2,1,height_ratios=[1.25,1],hspace=.28)
    ax=fig.add_subplot(gs[0]);ax.set_xlim(-15,15);ax.set_ylim(-1.7,2);ax.axis('off')
    ax.annotate('',xy=(14,0),xytext=(-14,0),arrowprops={'arrowstyle':'->','color':'#777777'});ax.text(14,.15,r'$k_x$')
    for x,r,s in [(-12,'P','-'),(-6,'A','-'),(6,'A','+'),(12,'P','+')]:
        col=BLUE if r=='A' else GREEN
        ax.add_patch(Ellipse((x,0),2.5,1.45,edgecolor=col,facecolor=col,alpha=.18))
        ax.text(x,0,f'{r}{s}',ha='center',va='center',color=col,weight='bold')
        ax.text(x,-1.17,r'$\Omega\ne0$' if r=='A' else r'$\Omega=0$',ha='center',fontsize=9)
    for a,b in [(-12,-6),(6,12)]:
        ax.add_patch(FancyArrowPatch((a,.8),(b,.8),connectionstyle='arc3,rad=-.35',arrowstyle='<->',mutation_scale=10,color=ORANGE))
    ax.text(0,1.5,r'(a) Intersector rate $\eta\Gamma_0$; intrasector rate scale $\Gamma_0$',ha='center')
    ax.text(0,-1.65,'Both pairs conduct; only A has direct Hall readout at zero passive mass.',ha='center',fontsize=9)
    ax=fig.add_subplot(gs[1]);ax.axis('off');ax.set_xlim(0,4);ax.set_ylim(0,1)
    labels=[('1  Active reference\ncurrent',r'$g_{A0}/\|g_{A0}\|$','Full uncoupled\ncurrent response'),('2  Active valley\ndistortion',r'$(I-u_1u_1^T)p_A$','Orthogonal\nvalley change'),('3  Passive\ncurrent',r'$b_P/\|b_P\|$','Driven longitudinal\nmotion'),('4  Passive valley\npopulation',r'$p_P/\|p_P\|$','Opposite partner\npopulations')]
    for i,(title,formula,desc) in enumerate(labels):
        col=BLUE if i<2 else GREEN
        ax.text(i+.5,.94,title,ha='center',va='top',fontsize=9,weight='bold',color=col)
        ax.text(i+.5,.49,formula,ha='center',va='center',fontsize=11)
        ax.text(i+.5,.04,desc,ha='center',va='bottom',fontsize=8)
        if i<3:ax.axvline(i+1,ymin=.1,ymax=.95,color='#dddddd')
    fig.text(.06,.44,'(b) Fixed order of added physical degrees of freedom',fontsize=9)
    save(fig,'figure_1_physical_modes')
    fig,axes=plt.subplots(2,2,figsize=(7.1,4.9),layout='constrained')
    for ax,label,title in zip(axes.flat,['weak','scalar','shape','decision_cancellation'],['Weak / small budget','Strong / scalar dominated','Contact / shape dominated','Finite range / fixed-budget dependence']):
        r=reps[label]
        for key,col,mark,leg in [('H',ORANGE,'o',r'$H_m/H$'),('sigma_xx',BLUE,'s',r'$\sigma_{xx,m}/\sigma_{xx}$'),('Q_A',GREEN,'^',r'$Q_{A,m}/Q_A$')]:
            values=[r[f'M{m}_{key}']/r['full_'+key] for m in range(1,5)]+[1.]
            ax.plot(range(1,6),values,marker=mark,color=col,label=leg,lw=1.3,ms=4)
        ax.axhline(1,color='#888888',ls=':',lw=.8);ax.set_xticks(range(1,6),['M1','M2','M3','M4','Full'])
        ax.set_title(title,fontsize=9);ax.set_ylabel('Response / full value');ax.set_ylim(.15,1.65)
        ax.text(.04,.07,f'M4 Hall error: {100*r["M4_error_H"]:.3g}%',transform=ax.transAxes,va='bottom',fontsize=8,bbox={'facecolor':'white','edgecolor':'none','alpha':.9})
        ax.grid(alpha=.15)
    axes[0,0].legend(loc='lower right',fontsize=7)
    save(fig,'figure_2_ablation')
    fig,axes=plt.subplots(2,2,figsize=(7.1,5.4),layout='constrained')
    colors={0:GRAY,.05:BLUE,.1:ORANGE}
    for kernel,ls in [('gaussian','-'),('lorentzian','--')]:
        for xi in [0,.05,.1]:
            rr=sorted([r for r in inv if r['sample']=='kernel_comparison' and r['kernel']==kernel and r['delta_P']==.4 and r['xi']==xi and r['eta']>0],key=lambda r:r['eta'])
            for ax,key in [(axes[0,0],'rho_rel'),(axes[0,1],'H_error')]:ax.loglog([r['eta'] for r in rr],[max(r[key],1e-15) for r in rr],ls,color=colors[xi],lw=1.3,label=rf'$\xi={xi:.2f}$' if kernel=='gaussian' else None)
    axes[0,0].set(title=r'(a) Invariant leakage, $\delta_P=0.4$',ylabel=r'$\rho_{\rm rel}$',xlabel=r'$\eta$');axes[0,0].legend()
    axes[0,1].set(title='(b) Hall projection error',ylabel=r'$|H_4-H|/|H|$',xlabel=r'$\eta$')
    axes[0,1].text(.04,.06,'Solid: Gaussian\nDashed: Lorentzian\nDisplay floor: $10^{-15}$',transform=axes[0,1].transAxes,fontsize=7)
    for xi in [.05,.1]:
        rr=[r for r in inv if r['sample']=='principal' and r['eta']>0 and r['xi']==xi]
        axes[1,0].loglog([r['rho_rel'] for r in rr],[max(r['H_error'],1e-15) for r in rr],'.',color=colors[xi],ms=2,alpha=.5,label=rf'$\xi={xi}$')
    from matplotlib.ticker import NullFormatter, FuncFormatter
    axes[1,0].set_xticks([.005,.01,.02,.05]);axes[1,0].xaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{v:g}'));axes[1,0].xaxis.set_minor_formatter(NullFormatter())
    axes[1,0].set(title='(c) Leakage is not an error bound',xlabel=r'$\rho_{\rm rel}$',ylabel=r'$|H_4-H|/|H|$')
    rr=sorted([r for r in spec if r['label']=='decision_cancellation'],key=lambda r:r['rate']);rate=np.array([r['rate'] for r in rr]);term=np.array([r['signed_error_contribution'] for r in rr])
    axes[1,1].plot(rate,np.cumsum(term),color=PURPLE,lw=1.5);axes[1,1].axhline(0,color=GRAY,lw=.8)
    axes[1,1].set(title='(d) Signed spectral error cancels',xlabel=r'Collision rate $\lambda$',ylabel=r'$\sum_{\lambda_\nu\leq\lambda} h_\nu r_\nu/\lambda_\nu$')
    axes[1,1].ticklabel_format(axis='both',style='sci',scilimits=(-2,2));axes[1,1].text(.97,.96,f'Final: {sum(term):.3g}\nSum of magnitudes: {sum(abs(term)):.3g}',transform=axes[1,1].transAxes,ha='right',va='top',fontsize=7)
    for ax in axes.flat:ax.grid(alpha=.15)
    save(fig,'figure_3_leakage')
    fig,axes=plt.subplots(2,2,figsize=(7.1,4.7),layout='constrained')
    for ax,label,title in zip(axes.flat,['scalar','shape','weak','decision_cancellation'],['Scalar dominated','Shape dominated','Small-budget safe at 1%','Cancellation-dependent at 1%']):
        r=next(r for r in coeff if r['label']==label);H0=r['H0'];S=r['full_S'];T=r['full_T'];a=abs(S)+abs(T);B=a/abs(H0);actual=abs(S+T)/abs(r['full_H']);env=B/(1-B) if B<1 else np.inf
        ax.bar(['Scalar S','Shape T','Sum'],[S/abs(H0),T/abs(H0),(S+T)/abs(H0)],color=[BLUE,ORANGE,GREEN],width=.6)
        ax.axhline(0,color='#555555',lw=.7);ax.set_title(title,fontsize=9);ax.set_ylabel(r'Correction / $|H_0|$');ax.margins(y=.3)
        ax.text(.98,.96,f'Budget B = {100*B:.3g}%\nActual error = {100*actual:.3g}%\nFixed-budget envelope = {100*env:.3g}%',transform=ax.transAxes,ha='right',va='top',fontsize=7,bbox={'facecolor':'white','edgecolor':'none','alpha':.9})
    save(fig,'figure_4_scalar_shape')
    fig,axes=plt.subplots(3,3,figsize=(7.1,6.1),sharex=True,sharey=True,layout='constrained')
    cmap=ListedColormap([BLUE,'#e3ae44','#d3d6dc']);codes={'ROBUSTLY_SAFE':0,'CANCELLATION_DEPENDENT_SAFE':1,'UNSAFE':2,'UNDEFINED':3}
    for i,tol in enumerate([.01,.05,.1]):
        for j,xi in enumerate([0,.05,.1]):
            rr=[r for r in decisions if r['tolerance']==tol and r['xi']==xi and r['eta']>0]
            x=sorted(set(r['eta'] for r in rr));y=sorted(set(r['delta_P'] for r in rr));lookup={(r['eta'],r['delta_P']):codes[r['decision']] for r in rr}
            z=np.array([[lookup[a,b] for a in x] for b in y]);axes[i,j].pcolormesh(x,y,z,cmap=cmap,vmin=0,vmax=2,shading='nearest',rasterized=True);axes[i,j].set_xscale('log')
            if i==0:axes[i,j].set_title(rf'$\xi={xi:.2f}$')
            if j==0:axes[i,j].set_ylabel(rf'{int(100*tol)}% tolerance'+'\n'+r'$\delta_P$')
            if i==2:axes[i,j].set_xlabel(r'$\eta>0$')
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(color=BLUE,label='Robustly safe'),Patch(color='#e3ae44',label='Cancellation-dependent under fixed-budget test'),Patch(color='#d3d6dc',label='Unsafe')],loc='outside upper center',ncol=1,fontsize=8)
    save(fig,'figure_5_decision_map')
    fig,axes=plt.subplots(2,2,figsize=(7.1,5.2),layout='constrained')
    for j,label in enumerate(['weak','decision_cancellation']):
        rr=[r for r in spec if r['label']==label];rates=[r['rate'] for r in rr]
        for key,col,leg in [('drive_overlap',BLUE,'Drive overlap'),('Hall_overlap',ORANGE,'Hall overlap')]:
            v=np.abs([r[key] for r in rr]);v=v/np.linalg.norm(v);axes[0,j].loglog(rates,np.maximum(v,1e-16),'.',color=col,label=leg,ms=3)
        axes[0,j].set(title=label.replace('_',' '),ylabel='Normalized absolute overlap');axes[0,j].legend(fontsize=7)
        axes[1,j].loglog(rates,np.maximum([r['leakage_fraction'] for r in rr],1e-20),'.',color=GREEN,ms=3,label='Fraction of Frobenius leakage power')
        axes[1,j].set(xlabel='Collision rate',ylabel='Leakage power fraction');axes[1,j].grid(alpha=.15)
    for ax in axes[:,1]:
        ax.set_xticks([.0003,.001,.002]);ax.set_xticklabels([r'$3\times10^{-4}$',r'$10^{-3}$',r'$2\times10^{-3}$']);ax.xaxis.set_minor_formatter(NullFormatter())
    save(fig,'figure_S5_spectral_weights')
    write_json(root/'figures/main_input_manifest.json',{'datasets':used,'generator':sha(Path(__file__)),'figures':{p.name:sha(p) for p in out.glob('*.pdf')}})
