exec(open('_ortak.py').read())
import os; os.makedirs('sekiller',exist_ok=True)
MD=S.groupby(['yontem','n_dagilimi','delta','tau2','k']).ret_orani.mean()
sel=[s for s in ST if s[0] in ('Fisher','Stouffer','MAT','PI','HKSJ')]
DN={'D1':'D$_1$ (homojen)','D2':'D$_2$ (sağa çarpık)','D3':'D$_3$ (iki büyük çalışma baskın)'}
def t2s(t): return '0' if t==0 else str(t).replace('.',',')
for dl in (0.1,0.2):
  for pair,tag in (((0,0.01),'a'),((0.04,0.16),'b')):
    fig,axs=plt.subplots(2,3,figsize=(8.9,6.2))
    for i,t in enumerate(pair):
        for j,D in enumerate(['D1','D2','D3']):
            ax=axs[i,j]
            for key,lab,kw in sel:
                kk=dict(kw); kk.pop('mew',None); kk['lw']=kk['lw']*0.8; kk['ms']=4.5
                ax.plot(X,[MD[(key,D,dl,t,k)] for k in KS],label=lab,mec='white',mew=0.6,**kk)
            ax.set_ylim(-0.05,1.04); ax.set_yticks([0,0.2,0.4,0.6,0.8,1.0]); ax.yaxis.set_major_formatter(fmt)
            ax.set_xticks(X); ax.set_xticklabels([str(k) for k in KS],fontsize=9); ax.set_xlim(-0.4,len(KS)-0.6)
            ax.tick_params(axis='y',labelsize=9.5)
            ax.grid(color='#e6e6e6',lw=0.9); ax.set_axisbelow(True)
            for s in ('top','right'): ax.spines[s].set_visible(False)
            if i==0: ax.set_title(DN[D],fontsize=11.5,pad=8)
            if j>0: ax.set_yticklabels([])
        axs[i,0].set_ylabel(f'τ² = {t2s(t)}\n\nGüç',fontsize=11)
    for ax in axs[-1,:]: ax.set_xlabel('Çalışma sayısı (k)',fontsize=10.5)
    h,l=axs[0,0].get_legend_handles_labels()
    fig.legend(h,l,loc='upper center',ncol=5,frameon=False,fontsize=11.5,handlelength=3.0,columnspacing=1.8,bbox_to_anchor=(0.5,1.0))
    fig.tight_layout(rect=(0,0,1,0.94),h_pad=1.2,w_pad=0.6)
    fig.savefig(f'sekiller/guc_d{str(dl).replace(".","")}_{tag}.png',dpi=300); plt.close(fig)
print('ok')
