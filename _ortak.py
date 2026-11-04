import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
plt.rcParams.update({'font.family':'Liberation Serif','font.size':11.5,'axes.edgecolor':'#555555','axes.linewidth':1.0})
S=pd.read_csv('sim_v7_R2000.csv'); S=S[S.grup=='p_birlestirici'][['k','n_dagilimi','delta','tau2','yontem','ret_orani']]
H=pd.read_csv('hibrit_R2000.csv').rename(columns={'Hibrit':'ret_orani'}); H['yontem']='Hibrit'
S=pd.concat([S,H[['k','n_dagilimi','delta','tau2','yontem','ret_orani']]])
M=S.groupby(['yontem','delta','tau2','k']).ret_orani.mean()
KS=[3,4,5,7,10,15,20,30,50]; X=np.arange(len(KS))
ST=[('Fisher','Fisher',dict(color='#4d4d4d',ls=(0,(5,2.5)),marker='o',lw=2.2,ms=6)),
    ('Stouffer','Stouffer',dict(color='#9a9a9a',ls=(0,(5,2.5)),marker='v',lw=2.0,ms=6)),
    ('n-Stouffer','Ağırlıklı Z (√n)',dict(color='#2b7bd6',ls=':',marker=None,lw=2.4)),
    ('GRPC','GRPC',dict(color='#8a6d1f',ls='-.',marker='+',lw=2.0,ms=8,mew=1.8)),
    ('Hibrit','Hibrit GRPC × Shannon',dict(color='#c9a227',ls='-.',marker='*',lw=2.0,ms=6)),
    ('MAT','MAT',dict(color='#e8603c',ls='-',marker='s',lw=3.0,ms=6)),
    ('PI','PI',dict(color='#1aaf7a',ls='-',marker='^',lw=3.0,ms=6)),
    ('HKSJ','HKSJ',dict(color='#7b4fbf',ls=':',marker='D',lw=2.4,ms=6))]
fmt=FuncFormatter(lambda v,_: f'{v:.2f}'.replace('.',','))
def draw(panels,ylabel,fname,nominal):
    fig,axs=plt.subplots(2,2,figsize=(8.9,6.97),sharex=False)
    for ax,(dl,t,title) in zip(axs.flat,panels):
        for key,lab,kw in ST:
            y=[M[(key,dl,t,k)] for k in KS]
            kk=dict(kw); mew=kk.pop('mew',0.8)
            if key=='GRPC': ax.plot(X,y,label=lab,mec=kk['color'],mew=mew,**kk)
            else: ax.plot(X,y,label=lab,mec='white',mew=mew,**kk)
        if nominal: ax.axhline(0.05,color='#d40000',lw=1.4,zorder=1)
        ax.set_ylim(-0.05,1.04); ax.set_yticks([0,0.2,0.4,0.6,0.8,1.0]); ax.yaxis.set_major_formatter(fmt)
        ax.set_xticks(X); ax.set_xticklabels([str(k) for k in KS]); ax.set_xlim(-0.4,len(KS)-0.6)
        ax.grid(color='#e6e6e6',lw=0.9); ax.set_axisbelow(True)
        for s in ('top','right'): ax.spines[s].set_visible(False)
        ax.set_title(title,fontsize=12)
    for ax in axs[:,0]: ax.set_ylabel(ylabel)
    for ax in axs[1,:]: ax.set_xlabel('Çalışma sayısı (k)')
    h,l=axs[0,0].get_legend_handles_labels()
    fig.legend(h,l,loc='upper center',ncol=4,frameon=False,fontsize=11.5,handlelength=3.2,columnspacing=2.2,bbox_to_anchor=(0.5,1.0))
    fig.tight_layout(rect=(0,0,1,0.915))
    fig.savefig(fname,dpi=300); plt.close(fig)
for k,l,kw in ST: kw.setdefault('mew',0.8)
