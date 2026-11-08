"""Küçük örneklemli çalışmalarla duyarlılık analizi.
Her çalışmada toplam örneklem n ~ U(10, 40) (yuvarlanmış), iki eşit grup.
Gerçek etkiler theta_i ~ N(delta, tau2); gözlenen etki Hedges g, varyansı veriden kestirilir.
Çalışma düzeyindeki tek yönlü p değerleri iki örneklem t testinden (sd = n - 2) hesaplanır.
Faktörler: k {3,5,10,20,50} x delta {0, 0.2} x tau2 {0, 0.01, 0.04, 0.16} = 40 senaryo, R = 2000.
Tohum = 20270000 + senaryo sırası.
"""
import sys, itertools, math, time
sys.argv=[sys.argv[0]]+sys.argv[1:2]
import numpy as np, pandas as pd
from scipy import stats, special
import sim_v7 as S
from hypo import hypo_sf
S.BASE_SEED=20270000
_P={}
def gen_meta(rng,k,D,delta,tau2,Sx):
    shape=(S.R,k)
    n=rng.uniform(10,40,shape).round(); n1=np.floor(n/2); n2=n-n1
    theta=delta+np.sqrt(tau2)*rng.standard_normal(shape)
    se=np.sqrt(1/n1+1/n2)
    diff=theta+se*rng.standard_normal(shape)
    df=n-2; s=np.sqrt(rng.chisquare(df,shape)/df)
    t=diff/(s*se); d=diff/s
    J=np.exp(special.gammaln(df/2)-0.5*np.log(df/2)-special.gammaln((df-1)/2))
    g=J*d; v=(n1+n2)/(n1*n2)+g**2/(2*(n1+n2))
    _P['p']=np.clip(stats.t.sf(t,df),1e-300,1-1e-16)
    return g,v,n
S.gen_meta=gen_meta
def pcomb(y,v,n,t2_dl,mu_dl,se_dl,hk_se):
    Rr,k=y.shape; p=_P['p']; zc=stats.norm.isf(p); lp=-np.log(p); out={}
    out['Fisher']=stats.chi2.sf(2*lp.sum(1),2*k)
    out['Stouffer']=stats.norm.sf(zc.sum(1)/math.sqrt(k))
    w=np.sqrt(n); out['n-Stouffer']=stats.norm.sf((w*zc).sum(1)/np.sqrt((w**2).sum(1)))
    ps=np.sort(p,1); rtp_null,ord_null=S.null_tables(k); i=np.arange(1,k+1)
    out['ordmeta']=S.emp_left_p(ord_null,special.betainc(i,k-i+1,ps).min(1))
    wg=np.exp(-(n-800.0)**2/(2*400.0**2)); Tg=(wg*lp).sum(1)
    out['GRPC']=np.array([hypo_sf(Tg[r],wg[r]) for r in range(Rr)])
    sh=n/n.sum(1,keepdims=True); wh=wg*(-sh*np.log(sh)); Th=(wh*lp).sum(1)
    out['Hibrit']=np.array([hypo_sf(Th[r],wh[r]) for r in range(Rr)])
    wm=np.minimum(1.0,np.median(n,1,keepdims=True)/n); Tm=(wm*lp).sum(1)
    out['MAT']=np.array([hypo_sf(Tm[r],wm[r]) for r in range(Rr)])
    out['PI']=stats.t.sf(mu_dl/np.sqrt(t2_dl+se_dl**2),k-2)
    out['HKSJ']=stats.t.sf(mu_dl/hk_se,k-1)
    return out,p
S.pcomb_all=pcomb
S.METHODS2=['Fisher','Stouffer','n-Stouffer','ordmeta','GRPC','Hibrit','MAT','PI','HKSJ']
K=[3,5,10,20,50]; DEL=[0.0,0.2]; TAU=[0.0,0.01,0.04,0.16]
def job(a):
    idx,(k,dl,t)=a
    return S.run_scenario(idx,k,'DK',dl,t,'S0')
if __name__=='__main__':
    from multiprocessing import Pool
    grid=list(itertools.product(K,DEL,TAU)); t0=time.time()
    with Pool() as pool: res=pool.map(job,list(enumerate(grid)))
    pd.DataFrame([r for rr in res for r in rr]).to_csv('sim_kucuk_R%d.csv'%S.R,index=False)
    print('bitti',len(grid),round(time.time()-t0))
