"""
sim_anlamsiz.py
===============
Anlamsız çalışma oranlarına ilişkin ek benzetim (Tablo 4-15 ve Tablo 4-16).

Veri üretim mekanizması ana benzetimle (sim_v7.py) aynıdır:
    theta_i ~ N(delta, tau2), v_i = 4 / n_i, y_i ~ N(theta_i, v_i),
    p_i = 1 - Phi(y_i / sqrt(v_i))  (tek yönlü)
Faktörler: k (9 düzey) x D (D1, D2, D3) x delta (0; 0,1; 0,2) x tau2 (4 düzey)
= 324 senaryo, senaryo başına 20.000 tekrar. Yayın yanlılığı yoktur.

Bütün senaryolar tek bir rastgele sayı üretecinden (tohum 7), sim_v7.py ile
aynı senaryo sırasında, ardışık olarak üretilir. Bu nedenle betik tek süreçte
çalıştırılmalıdır.

Her senaryoda kaydedilenler:
    hic_anlamsiz_yok : bütün çalışmaları anlamlı (p_i <= 0,05) olan meta analiz oranı
    hepsi_anlamsiz   : bütün çalışmaları anlamsız (p_i > 0,05) olan meta analiz oranı
    n_hepsi_anlamsiz : bu meta analizlerin sayısı
    ret_<yöntem>     : bu meta analizler arasında birleşik p <= 0,05 olanların sayısı
                       (Fisher, Stouffer, ağırlıklı Z (sqrt n), MAT, PI, HKSJ)

Tablo 4-15: hic_anlamsiz_yok, tau2 düzeyleri üzerinden ortalama (delta x D x k).
Tablo 4-16 A: hepsi_anlamsiz, D ve tau2 üzerinden ortalama (delta x k).
Tablo 4-16 B: ret sayılarının D ve tau2 üzerinden toplamı / n_hepsi_anlamsiz
             toplamı; toplam 100'ün altındaysa ret oranı verilmez.

Kullanım:  python3 sim_anlamsiz.py            # sim_anlamsiz_R20000.csv
           python3 tablolar_anlamsiz.py       # Tablo_4-15.csv, Tablo_4-16.csv
"""
import itertools, math, sys, time
import numpy as np
import pandas as pd
from scipy import stats
from hypo import hypo_sf
from sim_v7 import K_LEV, D_LEV, DELTA_LEV, TAU2_LEV, draw_n

R = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
SEED = 7
ALPHA = 0.05
YONTEMLER = ["Fisher", "Stouffer", "n-Stouffer", "MAT", "PI", "HKSJ"]


def pvals(y, v, n):
    """Seçilen meta analizlerde altı yöntemin birleşik p değerleri."""
    Rr, k = y.shape
    z = y / np.sqrt(v)
    p = np.clip(stats.norm.sf(z), 1e-300, 1 - 1e-16)
    zc = stats.norm.isf(p)
    lp = -np.log(p)
    out = {}
    out["Fisher"] = stats.chi2.sf(2 * lp.sum(1), 2 * k)
    out["Stouffer"] = stats.norm.sf(zc.sum(1) / math.sqrt(k))
    w = np.sqrt(n)
    out["n-Stouffer"] = stats.norm.sf((w * zc).sum(1) / np.sqrt((w ** 2).sum(1)))
    wm = np.minimum(1.0, np.median(n, 1, keepdims=True) / n)
    Tm = (wm * lp).sum(1)
    out["MAT"] = np.array([hypo_sf(Tm[r], wm[r]) for r in range(Rr)])
    # DerSimonian-Laird, PI ve HKSJ (sim_v7.py ile aynı)
    wf = 1.0 / v
    mu_fe = (wf * y).sum(1) / wf.sum(1)
    Q = (wf * (y - mu_fe[:, None]) ** 2).sum(1)
    c = wf.sum(1) - (wf ** 2).sum(1) / wf.sum(1)
    t2 = np.maximum(0.0, (Q - (k - 1)) / c)
    wr = 1.0 / (v + t2[:, None])
    mu = (wr * y).sum(1) / wr.sum(1)
    se = np.sqrt(1.0 / wr.sum(1))
    out["PI"] = stats.t.sf(mu / np.sqrt(t2 + se ** 2), k - 2)
    hk = (wr * (y - mu[:, None]) ** 2).sum(1) / ((k - 1) * wr.sum(1))
    out["HKSJ"] = stats.t.sf(mu / np.sqrt(np.maximum(hk, 1e-300)), k - 1)
    return out


def main():
    rng = np.random.default_rng(SEED)
    grid = list(itertools.product(K_LEV, D_LEV, DELTA_LEV, TAU2_LEV))
    rows, t0 = [], time.time()
    for idx, (k, D, delta, tau2) in enumerate(grid):
        n = draw_n(rng, D, (R, k))
        theta = delta + np.sqrt(tau2) * rng.standard_normal((R, k))
        v = 4.0 / n
        y = theta + np.sqrt(v) * rng.standard_normal((R, k))
        p = stats.norm.sf(y / np.sqrt(v))
        anlamli = p <= ALPHA
        hepsi_ns = ~anlamli.any(1)
        row = dict(senaryo=idx, k=k, n_dagilimi=D, delta=delta, tau2=tau2, R=R,
                   hic_anlamsiz_yok=float(anlamli.all(1).mean()),
                   hepsi_anlamsiz=float(hepsi_ns.mean()),
                   n_hepsi_anlamsiz=int(hepsi_ns.sum()))
        if hepsi_ns.any():
            pv = pvals(y[hepsi_ns], v[hepsi_ns], n[hepsi_ns])
            for m in YONTEMLER:
                row["ret_" + m] = int((pv[m] <= ALPHA).sum())
        else:
            for m in YONTEMLER:
                row["ret_" + m] = 0
        rows.append(row)
        if idx % 36 == 0:
            print(f"{idx + 1}/{len(grid)}  {time.time() - t0:.0f}s", flush=True)
    pd.DataFrame(rows).to_csv(f"sim_anlamsiz_R{R}.csv", index=False)
    print("bitti", round(time.time() - t0), "s")


if __name__ == "__main__":
    main()
