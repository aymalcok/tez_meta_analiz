"""
sim_bin.py — İkili veriyle duyarlılık analizi (2x2 tablodan türetilen benzetim)
Faktörler: k {3,5,10,20,50} x D {D1,D2,D3} x p_c {0.10,0.03} x RR {1,0.9,0.8} x tau2 {0,0.01,0.04,0.16} = 360 senaryo
Veri üretimi (çalışma i):
  n_i ~ D (ana benzetimle aynı), n1 = floor(n_i/2), n2 = n_i - n1
  ln RR_i ~ N(ln RR, tau2);  p_t = min(p_c * RR_i, 0.999)
  a ~ Bin(n1, p_t) (tedavi), c ~ Bin(n2, p_c) (kontrol)
Analiz:
  iki kolda da sıfır olay -> çalışma çıkarılır
  herhangi bir hücre sıfır -> tüm hücrelere 0,5 (a+0.5, c+0.5, n1+1, n2+1)
  y = -ln RR_hat (tedavi lehine pozitif), v = 1/a - 1/n1 + 1/c - 1/n2
  p_i = 1 - Phi(y/sqrt(v)) (tek yönlü, tedavi lehine)
Yöntemler: Fisher, Stouffer, Ağırlıklı Z (sqrt n), GRPC, Hibrit, MAT, PI (t_{k-2}), HKSJ (t_{k-1})
k_eff < 3 olan tekrarlar: hiçbir yöntem H0'ı reddetmemiş sayılır (oranı raporlanır).
Tohum = 20261020 + senaryo sırası. R = 2000.
"""
import itertools, math, sys, time
import numpy as np, pandas as pd
from scipy import stats
from hypo import hypo_sf
R = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
OUT = sys.argv[2] if len(sys.argv) > 2 else f"sim_bin_R{R}.csv"
BASE_SEED = 20261020; ALPHA = 0.05
K_LEV=[3,5,10,20,50]; D_LEV=["D1","D2","D3"]; PC_LEV=[0.10,0.03]; RR_LEV=[1.0,0.9,0.8]; TAU2_LEV=[0.0,0.01,0.04,0.16]
M=["Fisher","Stouffer","n-Stouffer","GRPC","Hibrit","MAT","PI","HKSJ"]
def draw_n(rng, D, shape):
    Rr, k = shape
    if D == "D1": return rng.uniform(80, 120, shape).round()
    if D == "D2": return np.clip(rng.lognormal(np.log(100), 1.0, shape), 20, 5000).round()
    n = np.clip(rng.lognormal(np.log(100), 0.8, shape), 20, 2000).round()
    n[:, -2:] = rng.uniform(10000, 17000, (Rr, 2)).round()
    return n
def one(y, v, n):
    k = y.size; out = {}
    z = y/np.sqrt(v); p = np.clip(stats.norm.sf(z), 1e-300, 1-1e-16); zc = stats.norm.isf(p); lp = -np.log(p)
    out["Fisher"] = stats.chi2.sf(2*lp.sum(), 2*k)
    out["Stouffer"] = stats.norm.sf(zc.sum()/math.sqrt(k))
    w = np.sqrt(n); out["n-Stouffer"] = stats.norm.sf((w*zc).sum()/math.sqrt((w**2).sum()))
    wg = np.exp(-(n-800.0)**2/(2*400.0**2)); out["GRPC"] = hypo_sf((wg*lp).sum(), wg)
    s = n/n.sum(); wh = wg*(-s*np.log(s)); out["Hibrit"] = hypo_sf((wh*lp).sum(), wh)
    wm = np.minimum(1.0, np.median(n)/n); out["MAT"] = hypo_sf((wm*lp).sum(), wm)
    wf = 1/v; mu = (wf*y).sum()/wf.sum(); Q = (wf*(y-mu)**2).sum()
    c = wf.sum() - (wf**2).sum()/wf.sum(); t2 = max(0.0, (Q-(k-1))/c)
    wr = 1/(v+t2); mur = (wr*y).sum()/wr.sum(); se = math.sqrt(1/wr.sum())
    out["PI"] = stats.t.sf(mur/math.sqrt(t2+se**2), k-2)
    hk = (wr*(y-mur)**2).sum()/((k-1)*wr.sum()); out["HKSJ"] = stats.t.sf(mur/math.sqrt(max(hk,1e-300)), k-1)
    I2 = max(0.0, (Q-(k-1))/Q) if Q > 0 else 0.0
    return out, I2
def run(idx, k, D, pc, RR, tau2):
    rng = np.random.default_rng(BASE_SEED+idx)
    n = draw_n(rng, D, (R, k)); n1 = np.floor(n/2); n2 = n-n1
    lrr = math.log(RR) + math.sqrt(tau2)*rng.standard_normal((R, k))
    pt = np.minimum(pc*np.exp(lrr), 0.999)
    a = rng.binomial(n1.astype(int), pt).astype(float); cc = rng.binomial(n2.astype(int), pc).astype(float)
    dz = (a == 0) & (cc == 0)
    zc_ = (a == 0) | (cc == 0) | (a == n1) | (cc == n2)
    a2 = np.where(zc_, a+0.5, a); c2 = np.where(zc_, cc+0.5, cc); m1 = np.where(zc_, n1+1, n1); m2 = np.where(zc_, n2+1, n2)
    y = -np.log((a2/m1)/(c2/m2)); v = 1/a2 - 1/m1 + 1/c2 - 1/m2
    rej = {m: 0 for m in M}; I2s = []; few = 0
    for r in range(R):
        keep = ~dz[r]
        if keep.sum() < 3: few += 1; continue
        out, I2 = one(y[r, keep], v[r, keep], n[r, keep]); I2s.append(I2)
        for m in M: rej[m] += out[m] < ALPHA
    base = dict(senaryo=idx, k=k, n_dagilimi=D, p_c=pc, RR=RR, tau2=tau2, R=R,
                sifir_hucre_oran=float(zc_[~dz].mean()) if (~dz).any() else np.nan,
                cift_sifir_oran=float(dz.mean()), keff_3alti_oran=few/R,
                I2_ortanca=float(np.median(I2s)) if I2s else np.nan)
    return [dict(base, yontem=m, ret_orani=rej[m]/R, MCSE_ret=math.sqrt((rej[m]/R)*(1-rej[m]/R)/R)) for m in M]
def _job(a): return run(a[0], *a[1])
if __name__ == "__main__":
    from multiprocessing import Pool
    grid = list(itertools.product(K_LEV, D_LEV, PC_LEV, RR_LEV, TAU2_LEV))
    t0 = time.time(); rows = []
    with Pool(int(sys.argv[3]) if len(sys.argv) > 3 else 4) as pool:
        for i, rr in enumerate(pool.imap_unordered(_job, list(enumerate(grid)), chunksize=1)):
            rows += rr
            if i % 10 == 0: print(f"{i+1}/{len(grid)} {time.time()-t0:.0f}s", flush=True)
    pd.DataFrame(rows).sort_values(["senaryo","yontem"]).to_csv(OUT, index=False); print("bitti", time.time()-t0)
