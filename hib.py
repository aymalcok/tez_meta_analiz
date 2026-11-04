import itertools, time, sys, numpy as np, pandas as pd
from scipy import stats
import sim_v7 as S
from hypo import hypo_sf
def job(a):
    idx,(k,D,delta,tau2,Sx)=a
    rng=np.random.default_rng(S.BASE_SEED+idx)
    y,v,n=S.gen_meta(rng,k,D,delta,tau2,Sx)
    Rr=y.shape[0]
    p=np.clip(stats.norm.sf(y/np.sqrt(v)),1e-300,1-1e-16); lp=-np.log(p)
    wg=np.exp(-(n-800.0)**2/(2*400.0**2))
    s=n/n.sum(1,keepdims=True); we=-s*np.log(s); wh=wg*we
    Th=(wh*lp).sum(1); ph=np.array([hypo_sf(Th[r],wh[r]) for r in range(Rr)])
    Tg=(wg*lp).sum(1); pg=np.array([hypo_sf(Tg[r],wg[r]) for r in range(Rr)])
    return dict(senaryo=idx,k=k,n_dagilimi=D,delta=delta,tau2=tau2,Hibrit=(ph<S.ALPHA).mean(),GRPC_chk=(pg<S.ALPHA).mean())
if __name__=='__main__':
    from multiprocessing import Pool
    grid=list(itertools.product(S.K_LEV,S.D_LEV,S.DELTA_LEV,S.TAU2_LEV,S.S_LEV))
    t0=time.time(); rows=[]
    with Pool(2) as pool:
        for i,r in enumerate(pool.imap(job,list(enumerate(grid)),chunksize=1)):
            rows.append(r)
            if i%24==0: print(i,len(grid),round(time.time()-t0),flush=True)
    pd.DataFrame(rows).to_csv('hibrit_R2000.csv',index=False); print('bitti',time.time()-t0)
