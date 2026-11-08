"""Küçük örneklemli duyarlılık analizi: tau2 düzeylerine göre I. tip hata oranları (Şekil 4-10)."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
plt.rcParams.update({'font.family':'Liberation Serif','font.size':11.5,'axes.edgecolor':'#555555','axes.linewidth':1.0})
S=pd.read_csv('sim_kucuk_R2000.csv'); S=S[(S.grup=='p_birlestirici')&(S.delta==0)]
M=S.set_index(['yontem','tau2','k']).ret_orani
KS=[3,5,10,20,50]; X=np.arange(len(KS))
ST=[('Fisher','Fisher',dict(color='#4d4d4d',ls=(0,(5,2.5)),marker='o',lw=2.2,ms=6)),
    ('Stouffer','Stouffer',dict(color='#9a9a9a',ls=(0,(5,2.5)),marker='v',lw=2.0,ms=6)),
    ('n-Stouffer','Ağırlıklı Z (√n)',dict(color='#2b7bd6',ls=':',marker=None,lw=2.4)),
    ('GRPC','GRPC',dict(color='#8a6d1f',ls='-.',marker='+',lw=2.0,ms=8)),
    ('Hibrit','Hibrit GRPC × Shannon',dict(color='#c9a227',ls='-.',marker='*',lw=2.0,ms=6)),
    ('MAT','MAT',dict(color='#e8603c',ls='-',marker='s',lw=3.0,ms=6)),
    ('PI','PI',dict(color='#1aaf7a',ls='-',marker='^',lw=3.0,ms=6)),
    ('HKSJ','HKSJ',dict(color='#7b4fbf',ls=':',marker='D',lw=2.4,ms=6))]
fmt=FuncFormatter(lambda v,_: f'{v:.2f}'.replace('.',','))
fig,axs=plt.subplots(2,2,figsize=(8.9,6.97))
for ax,t,title in zip(axs.flat,[0.0,0.01,0.04,0.16],['τ² = 0','τ² = 0,01','τ² = 0,04','τ² = 0,16']):
    for key,lab,kw in ST:
        y=[M[(key,t,k)] for k in KS]
        mec=kw['color'] if key=='GRPC' else 'white'
        ax.plot(X,y,label=lab,mec=mec,mew=1.8 if key=='GRPC' else 0.8,**kw)
    ax.axhline(0.05,color='#d40000',lw=1.4,zorder=1)
    ax.set_ylim(-0.02,0.56); ax.set_yticks([0,0.1,0.2,0.3,0.4,0.5]); ax.yaxis.set_major_formatter(fmt)
    ax.set_xticks(X); ax.set_xticklabels([str(k) for k in KS]); ax.set_xlim(-0.3,len(KS)-0.7)
    ax.grid(color='#e6e6e6',lw=0.9); ax.set_axisbelow(True)
    for s in ('top','right'): ax.spines[s].set_visible(False)
    ax.set_title(title,fontsize=12)
for ax in axs[:,0]: ax.set_ylabel('I. tip hata oranı')
for ax in axs[1,:]: ax.set_xlabel('Çalışma sayısı (k)')
h,l=axs[0,0].get_legend_handles_labels()
fig.legend(h,l,loc='upper center',ncol=4,frameon=False,fontsize=11.5,handlelength=3.2,columnspacing=2.2,bbox_to_anchor=(0.5,1.0))
fig.tight_layout(rect=(0,0,1,0.915))
import os; os.makedirs('sekiller',exist_ok=True); fig.savefig('sekiller/Sekil_4-10_kucuk_orneklem.png',dpi=300)
