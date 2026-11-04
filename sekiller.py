import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Rectangle
plt.rcParams.update({'font.family':'Liberation Serif','font.size':11.5,'axes.edgecolor':'#555555','axes.linewidth':1.0})
exec(open('_ortak.py').read())
for k,l,kw in ST: kw.setdefault('mew',0.8)
OUT='sekiller/'
import os; os.makedirs(OUT,exist_ok=True)
from decimal import Decimal, ROUND_HALF_UP
def f3(x): return str(Decimal(repr(round(float(x),6))).quantize(Decimal('0.001'),ROUND_HALF_UP)).replace('.',',')
# ---- old figures redrawn
draw([(0,0,'τ² = 0'),(0,0.01,'τ² = 0,01'),(0,0.04,'τ² = 0,04'),(0,0.16,'τ² = 0,16')],'I. tip hata oranı',OUT+'s_tip1.png',True)
# ---- new 1: real data decisions
meth=['Fisher','Stouffer','Tippett','Ağırlıklı Z (√n)','Lancaster','wFisher','TPM','RTP','ordmeta','GRPC','Hibrit GRPC × Shannon','MAT','PI','HKSJ']
P={ 'Fisher':['<0,0001','<0,0001','<0,0001','<0,0001'],'Stouffer':['<0,0001','<0,0001','<0,0001','<0,0001'],
 'Tippett':['<0,0001','0,0027','0,0592','0,0002'],'Ağırlıklı Z (√n)':['<0,0001','<0,0001','0,6589','<0,0001'],
 'Lancaster':['<0,0001','<0,0001','0,6382','<0,0001'],'wFisher':['<0,0001','<0,0001','0,3116','<0,0001'],
 'TPM':['<0,0001','0,0001','<0,0001','<0,0001'],'RTP':['<0,0001','<0,0001','<0,0001','<0,0001'],
 'ordmeta':['<0,0001','0,0019','<0,0001','<0,0001'],'GRPC':['<0,0001','0,0048','<0,0001','0,0019'],
 'Hibrit GRPC × Shannon':['0,0001','0,0027','0,0004','0,0076'],'MAT':['0,0014','0,0045','<0,0001','0,0004'],
 'PI':['0,0136','0,0015','0,0915','0,1711'],'HKSJ':['<0,0001','<0,0001','0,0009','0,0274']}
truth=[True,True,False,False]
def pv(s): return 0.00001 if s.startswith('<') else float(s.replace(',','.'))
fig,ax=plt.subplots(figsize=(8.9,7.4))
rows=[]; 
groups=[('Yerleşik yöntemler',meth[:9]),('Keşifsel yaklaşımlar',meth[9:11]),('Önerilen yöntemler',meth[11:13]),('Karşılaştırma ölçütü',meth[13:])]
y=0; ylab=[]
for g,ms in groups:
    ylab.append((y,g,True)); y+=1
    for m in ms: ylab.append((y,m,False)); y+=1
H=y
cols=['Streptokinaz','Beta-bloker','Magnezyum','Diüretik']
GOOD='#d7efe2'; BAD='#f6d3cc'
for yy,lab,isg in ylab:
    Y=H-1-yy
    if isg:
        ax.text(-0.05,Y+0.5,lab,ha='right',va='center',fontsize=11,fontweight='bold'); continue
    ax.text(-0.05,Y+0.5,'   '+lab,ha='right',va='center',fontsize=11)
    ok=0
    for j in range(4):
        rej=pv(P[lab][j])<0.05; corr=(rej==truth[j]); ok+=corr
        ax.add_patch(Rectangle((j+0.04,Y+0.06),0.92,0.88,color=GOOD if corr else BAD,lw=0))
        ax.text(j+0.5,Y+0.5,P[lab][j]+('' if corr else r'$^{\dagger}$'),ha='center',va='center',fontsize=10.5,color='#1a1a1a' if corr else '#8a1c0c')
    ax.text(4.45,Y+0.5,f'{ok}/4',ha='center',va='center',fontsize=11,fontweight='bold' if ok==4 else 'normal')
for j,c in enumerate(cols):
    ax.text(j+0.5,H+0.85,c,ha='center',va='center',fontsize=11,fontweight='bold')
    ax.text(j+0.5,H+0.3,'Etki var' if truth[j] else 'Etki yok',ha='center',va='center',fontsize=10,color='#444444')
ax.text(4.45,H+0.55,'Doğru\nkarar',ha='center',va='center',fontsize=10.5,fontweight='bold')
ax.set_xlim(-0.05,4.8); ax.set_ylim(-0.1,H+1.2); ax.axis('off')
ax.add_patch(Rectangle((0.6,-1.05),0.25,0.5,color=GOOD,clip_on=False)); ax.text(0.9,-0.8,'Doğru karar',va='center',fontsize=10.5)
ax.add_patch(Rectangle((2.1,-1.05),0.25,0.5,color=BAD,clip_on=False)); ax.text(2.4,-0.8,r'Yanlış karar ($^{\dagger}$)',va='center',fontsize=10.5)
fig.subplots_adjust(left=0.27,right=0.99,top=0.98,bottom=0.07)
fig.savefig(OUT+'s_gercek.png',dpi=300); plt.close(fig)
# ---- new 2: D3 Fisher/MAT/PI/HKSJ type I
MD=S.groupby(['yontem','n_dagilimi','delta','tau2','k']).ret_orani.mean()
sel=[s for s in ST if s[0] in ('Fisher','MAT','PI','HKSJ')]
fig,axs=plt.subplots(2,2,figsize=(8.9,6.97))
for ax,(t,title) in zip(axs.flat,[(0,'τ² = 0'),(0.01,'τ² = 0,01'),(0.04,'τ² = 0,04'),(0.16,'τ² = 0,16')]):
    for key,lab,kw in sel:
        yv=[MD[(key,'D3',0,t,k)] for k in KS]; kk=dict(kw); mew=kk.pop('mew',0.8)
        ax.plot(X,yv,label=lab,mec='white',mew=mew,**kk)
    ax.axhline(0.05,color='#d40000',lw=1.4,zorder=1)
    ax.set_ylim(-0.05,1.04); ax.set_yticks([0,0.2,0.4,0.6,0.8,1.0]); ax.yaxis.set_major_formatter(fmt)
    ax.set_xticks(X); ax.set_xticklabels([str(k) for k in KS]); ax.set_xlim(-0.4,len(KS)-0.6)
    ax.grid(color='#e6e6e6',lw=0.9); ax.set_axisbelow(True)
    for s in ('top','right'): ax.spines[s].set_visible(False)
    ax.set_title(title,fontsize=12)
for ax in axs[:,0]: ax.set_ylabel('I. tip hata oranı')
for ax in axs[1,:]: ax.set_xlabel('Çalışma sayısı (k)')
h,l=axs[0,0].get_legend_handles_labels()
fig.legend(h,l,loc='upper center',ncol=4,frameon=False,fontsize=11.5,handlelength=3.2,columnspacing=2.2,bbox_to_anchor=(0.5,1.0))
fig.tight_layout(rect=(0,0,1,0.93)); fig.savefig(OUT+'s_d3.png',dpi=300); plt.close(fig)
# ---- new 3: distribution of type I error over all 108 null scenarios
N=S[S.delta==0]
order=['Fisher','Stouffer','Tippett','n-Stouffer','Lancaster','wFisher','TPM','RTP','ordmeta','GRPC','Hibrit','MAT','PI','HKSJ']
names=dict(zip(order,meth))
fig,ax=plt.subplots(figsize=(8.9,6.2))
colmap={k:kw['color'] for k,l,kw in ST}
for i,m in enumerate(order[::-1]):
    v=N[N.yontem==m].ret_orani.values
    c=colmap.get(m,'#7a7a7a')
    jit=(np.random.default_rng(1).random(len(v))-0.5)*0.36
    ax.scatter(v,i+jit,s=11,color=c,alpha=0.55,lw=0)
    ax.plot([np.median(v)]*2,[i-0.32,i+0.32],color='#1a1a1a',lw=2.2)
    ax.text(1.02,i,f3(np.max(v)),va='center',fontsize=10.5)
ax.axvline(0.05,color='#d40000',lw=1.4,zorder=0)
ax.set_yticks(range(len(order))); ax.set_yticklabels([names[m] for m in order[::-1]])
ax.set_xlim(-0.02,1.0); ax.xaxis.set_major_formatter(fmt); ax.set_xticks([0,0.2,0.4,0.6,0.8,1.0])
ax.set_xlabel('I. tip hata oranı'); ax.text(1.02,len(order)-0.3,'En\nyüksek',fontsize=10.5,va='bottom')
ax.grid(axis='x',color='#e6e6e6',lw=0.9); ax.set_axisbelow(True)
for s in ('top','right'): ax.spines[s].set_visible(False)
fig.subplots_adjust(left=0.24,right=0.89,top=0.93,bottom=0.09)
fig.savefig(OUT+'s_ozet.png',dpi=300); plt.close(fig)
print(N.groupby('yontem').ret_orani.agg(['median','max']).round(3))
