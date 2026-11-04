"""
sim_rmse_kapsama.py
====================
Tez benzetim modülü: etki kestiricileri için yanlılık, göreli yanlılık, RMSE ve
kapsama oranı; p birleştirme yöntemleri için I. tip hata ve güç.

TİK 12.05.2026 önerisi (RMSE, relatif bias, tahmin aralığı kapsama oranı) ve tez
önerisindeki "anlamsız çalışma oranı / hiç anlamsız çalışma yok" senaryolarını
karşılar.

ADEMP
-----
Veri üretimi (her çalışma i):
    theta_i ~ N(delta, tau2)             gerçek çalışma etkisi
    v_i     = 4 / n_i                    SMD örnekleme varyansı (iki eşit kol)
    y_i     ~ N(theta_i, v_i)            gözlenen etki
    p_i     = 1 - Phi(y_i / sqrt(v_i))   tek yönlü p (H1: etki > 0)

Faktörler (4 x 3 x 3 x 4 x 3 = 432 senaryo):
    k      : 5, 10, 20, 33
    D      : D1 homojen  n ~ U(80,120)
             D2 tipik    n ~ lognormal(ln 100, 1.0), [20, 5000] aralığına kırpılmış
             D3 baskın   k-2 küçük n ~ lognormal(ln 100, 0.8) [20, 2000]
                         + 2 mega-çalışma n ~ U(10000, 17000)
    delta  : 0, 0.1, 0.2
    tau2   : 0, 0.01, 0.04, 0.16
    S      : S0 seçici yayın yok
             S1 anlamsız (p >= 0.05) çalışma 0.5 olasılıkla yayımlanır
             S2 yalnızca p < 0.05 olan çalışmalar yayımlanır (hiç anlamsız yok)
    Seçici yayında her çalışma "yuvası" kabul edilene dek yeniden üretilir;
    mega-çalışma yuvaları mega-çalışma olarak kalır.

Hedef nicelik: delta (popülasyon ortalama etkisi). Grup 2 için H0: delta = 0.

Yöntemler
    Grup 1 (kestirici): FE (ters varyans), RE-DL, RE-REML, HKSJ (DL tabanlı) GA,
                        PI (DL, t_{k-2}) — PI'nin merkezi RE-DL kestirimidir.
    Grup 2 (p birleştirici, tek yönlü):
        Fisher, Stouffer, n-Stouffer (w=sqrt n), Liptak (w=n), Tippett, Simes,
        Lancaster (df=n_i), wFisher (Yoon 2021), TPM (tau=0.05),
        RTP (K=ceil(k/2)), ordmeta, CCT, HMP (Landau),
        GRPC (Gauss ağırlık mu=800, sigma=400; moment eşlemeli gama),
        MAT (saddlepoint), PI testi (tek yönlü t_{k-2}), HKSJ testi (tek yönlü t_{k-1}).

Performans
    Grup 1: yanlılık, göreli yanlılık (delta != 0), RMSE, ampirik SH,
            GA kapsama (delta), PI kapsama (yeni çalışmanın gerçek etkisi
            theta_new ~ N(delta, tau2); koşullu olasılık analitik olarak
            hesaplanır), ortalama aralık genişliği, MCSE.
    Grup 2: ret oranı (delta=0 iken I. tip hata, delta>0 iken güç), MCSE.
    Ek: yayımlanan çalışmalarda anlamsız çalışma oranı, gerçekleşen I² ortancası.

R = 2000 tekrar / senaryo. Tohum = 20260924 + senaryo sırası.
"""
import itertools, math, sys, time
import numpy as np
from hypo import hypo_sf
import pandas as pd
from scipy import stats, special
from scipy.optimize import brentq

R = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
BASE_SEED = 20260924
ALPHA = 0.05
Z975 = stats.norm.ppf(0.975)

K_LEV = [3, 4, 5, 7, 10, 15, 20, 30, 50]
D_LEV = ["D1", "D2", "D3"]
DELTA_LEV = [0.0, 0.1, 0.2]
TAU2_LEV = [0.0, 0.01, 0.04, 0.16]
S_LEV = ["S0"]

# --------------------------------------------------------------- veri üretimi
def draw_n(rng, D, shape):
    """shape = (R, k) -> örneklem büyüklükleri (yuva tipine göre)."""
    Rr, k = shape
    if D == "D1":
        return rng.uniform(80, 120, shape).round()
    if D == "D2":
        return np.clip(rng.lognormal(np.log(100), 1.0, shape), 20, 5000).round()
    n = np.clip(rng.lognormal(np.log(100), 0.8, shape), 20, 2000).round()
    n[:, -2:] = rng.uniform(10000, 17000, (Rr, 2)).round()
    return n


def gen_meta(rng, k, D, delta, tau2, S):
    shape = (R, k)
    n = draw_n(rng, D, shape)
    theta = delta + np.sqrt(tau2) * rng.standard_normal(shape)
    v = 4.0 / n
    y = theta + np.sqrt(v) * rng.standard_normal(shape)
    if S == "S0":
        return y, v, n
    keep_prob_ns = 0.5 if S == "S1" else 0.0
    z = y / np.sqrt(v)
    sig = z > stats.norm.isf(ALPHA)
    acc = sig | (rng.random(shape) < keep_prob_ns)
    it = 0
    while not acc.all():
        it += 1
        if it > 20000:
            raise RuntimeError("seçici yayın döngüsü yakınsamadı")
        idx = ~acc
        n_new = draw_n(rng, D, shape)
        th_new = delta + np.sqrt(tau2) * rng.standard_normal(shape)
        v_new = 4.0 / n_new
        y_new = th_new + np.sqrt(v_new) * rng.standard_normal(shape)
        n[idx], v[idx], y[idx] = n_new[idx], v_new[idx], y_new[idx]
        z = y / np.sqrt(v)
        sig_new = z > stats.norm.isf(ALPHA)
        acc_new = sig_new | (rng.random(shape) < keep_prob_ns)
        acc = acc | (idx & acc_new)
    return y, v, n

# --------------------------------------------------------------- Grup 1
def fe(y, v):
    w = 1.0 / v
    mu = (w * y).sum(1) / w.sum(1)
    return mu, np.sqrt(1.0 / w.sum(1)), w


def dl_tau2(y, v):
    mu, _, w = fe(y, v)
    Q = (w * (y - mu[:, None]) ** 2).sum(1)
    k = y.shape[1]
    c = w.sum(1) - (w ** 2).sum(1) / w.sum(1)
    return np.maximum(0.0, (Q - (k - 1)) / c), Q


def reml_tau2(y, v, t0, iters=100):
    t = t0.copy()
    for _ in range(iters):
        w = 1.0 / (v + t[:, None])
        mu = (w * y).sum(1) / w.sum(1)
        num = (w ** 2 * ((y - mu[:, None]) ** 2 - v)).sum(1) / (w ** 2).sum(1) + 1.0 / w.sum(1)
        t_new = np.maximum(0.0, num)
        if np.max(np.abs(t_new - t)) < 1e-10:
            t = t_new
            break
        t = t_new
    return t


def re_est(y, v, t2):
    w = 1.0 / (v + t2[:, None])
    mu = (w * y).sum(1) / w.sum(1)
    se = np.sqrt(1.0 / w.sum(1))
    return mu, se, w

# --------------------------------------------------------------- Grup 2 yardımcıları
NULL_CACHE = {}

def null_tables(k, B=2_000_000, seed=7):
    """RTP ve ordmeta için H0 altındaki dağılım (yalnızca k'ye bağlı)."""
    if k in NULL_CACHE:
        return NULL_CACHE[k]
    rng = np.random.default_rng(seed + k)
    K = math.ceil(k / 2)
    i = np.arange(1, k + 1)
    rtp, ordm = [], []
    for _ in range(B // 200_000):
        u = np.sort(rng.random((200_000, k)), axis=1)
        rtp.append(np.log(u[:, :K]).sum(1))
        ordm.append(special.betainc(i, k - i + 1, u).min(1))
    NULL_CACHE[k] = (np.sort(np.concatenate(rtp)), np.sort(np.concatenate(ordm)))
    return NULL_CACHE[k]


def emp_left_p(null_sorted, obs):
    """P(null <= obs)."""
    return (np.searchsorted(null_sorted, obs, side="right") + 1) / (null_sorted.size + 1)


def tpm_p(p, t=0.05):
    Rr, k = p.shape
    lw = np.where(p <= t, np.log(p), 0.0).sum(1)
    has = (p <= t).any(1)
    tot = np.zeros(Rr)
    lt = math.log(t)
    for r in range(1, k + 1):
        c = math.comb(k, r) * (1 - t) ** (k - r)
        cond = lw <= r * lt
        x = r * lt - lw
        s = np.zeros(Rr)
        term = np.ones(Rr)
        for j in range(r):
            if j > 0:
                term = term * x / j
            s += term
        tot += np.where(cond, c * np.exp(lw) * s, c * t ** r)
    return np.where(has, np.minimum(tot, 1.0), 1.0)


def satt_saddle(T, w):
    """Ağırlıklı üstel toplam için Lugannani–Rice (MAT, eşitlik 11)."""
    t = float(T); mean = w.sum()
    if abs(t - mean) < 1e-9 * max(1.0, mean):
        k2 = float(np.sum(w ** 2)); k3 = 2.0 * float(np.sum(w ** 3))
        return 0.5 - k3 / (6.0 * math.sqrt(2 * math.pi) * k2 ** 1.5)
    Kp = lambda z: np.sum(w / (1.0 - w * z))
    K = lambda z: -np.sum(np.log1p(-w * z))
    Kpp = lambda z: np.sum(w ** 2 / (1.0 - w * z) ** 2)
    hi = 1.0 / w.max() - 1e-12
    lo = -1.0
    while Kp(lo) > t:
        lo *= 2.0
    zh = brentq(lambda z: Kp(z) - t, lo, hi, xtol=1e-14)
    r = np.sign(zh) * math.sqrt(max(2.0 * (zh * t - K(zh)), 0.0))
    u = zh * math.sqrt(Kpp(zh))
    if abs(r) < 1e-8 or abs(u) < 1e-8:
        return 0.5
    return float(stats.norm.sf(r) - stats.norm.pdf(r) * (1.0 / r - 1.0 / u))


def pcomb_all(y, v, n, t2_dl, mu_dl, se_dl, hk_se):
    Rr, k = y.shape
    z = y / np.sqrt(v)
    p = np.clip(stats.norm.sf(z), 1e-300, 1 - 1e-16)
    zc = stats.norm.isf(p)
    out = {}
    lp = -np.log(p)
    out["Fisher"] = stats.chi2.sf(2 * lp.sum(1), 2 * k)
    out["Stouffer"] = stats.norm.sf(zc.sum(1) / math.sqrt(k))
    w = np.sqrt(n)
    out["n-Stouffer"] = stats.norm.sf((w * zc).sum(1) / np.sqrt((w ** 2).sum(1)))
    w = n
    out["Liptak"] = stats.norm.sf((w * zc).sum(1) / np.sqrt((w ** 2).sum(1)))
    pmin = p.min(1)
    out["Tippett"] = -np.expm1(k * np.log1p(-pmin))
    ps = np.sort(p, 1)
    out["Simes"] = np.minimum(1, (k * ps / np.arange(1, k + 1)).min(1))
    out["Lancaster"] = stats.chi2.sf(stats.chi2.isf(p, n).sum(1), n.sum(1))
    wy = k * n / n.sum(1, keepdims=True)
    out["wFisher"] = stats.gamma.sf(stats.gamma.isf(p, wy, scale=2).sum(1), k, scale=2)
    out["TPM"] = tpm_p(p, 0.05)
    rtp_null, ord_null = null_tables(k)
    K = math.ceil(k / 2)
    crtp = np.array([min(j, K) / j for j in range(1, k + 1)])
    Trtp = -np.log(ps[:, :K]).sum(1)
    out["RTP"] = np.array([hypo_sf(Trtp[r], crtp) for r in range(Rr)])
    i = np.arange(1, k + 1)
    out["ordmeta"] = emp_left_p(ord_null, special.betainc(i, k - i + 1, ps).min(1))
    Tc = np.tan((0.5 - p) * np.pi).mean(1)
    out["CCT"] = 0.5 - np.arctan(Tc) / np.pi
    ph = k / (1.0 / p).sum(1)
    out["HMP"] = np.minimum(1, stats.landau.sf(1.0 / ph, loc=np.log(k) + 0.874, scale=np.pi / 2))
    # GRPC (hipoüstel, tam)
    wg = np.exp(-(n - 800.0) ** 2 / (2 * 400.0 ** 2))
    Tg = (wg * lp).sum(1)
    out["GRPC"] = np.array([hypo_sf(Tg[r], wg[r]) for r in range(Rr)])
    # MAT (hipoüstel, tam)
    wm = np.minimum(1.0, np.median(n, 1, keepdims=True) / n)
    Tm = (wm * lp).sum(1)
    out["MAT"] = np.array([hypo_sf(Tm[r], wm[r]) for r in range(Rr)])
    # PI testi (tek yönlü, t_{k-2})
    out["PI"] = stats.t.sf(mu_dl / np.sqrt(t2_dl + se_dl ** 2), k - 2)
    out["HKSJ"] = stats.t.sf(mu_dl / hk_se, k - 1)
    return out, p

METHODS2 = ["Fisher", "Stouffer", "n-Stouffer", "Liptak", "Tippett", "Simes",
            "Lancaster", "wFisher", "TPM", "RTP", "ordmeta", "CCT", "HMP",
            "GRPC", "MAT", "PI", "HKSJ"]

# --------------------------------------------------------------- tek senaryo
def run_scenario(idx, k, D, delta, tau2, S):
    rng = np.random.default_rng(BASE_SEED + idx)
    y, v, n = gen_meta(rng, k, D, delta, tau2, S)
    # Grup 1
    mu_fe, se_fe, w_fe = fe(y, v)
    t2_dl, Q = dl_tau2(y, v)
    mu_dl, se_dl, w_dl = re_est(y, v, t2_dl)
    t2_re = reml_tau2(y, v, t2_dl)
    mu_re, se_re, _ = re_est(y, v, t2_re)
    hk_var = (w_dl * (y - mu_dl[:, None]) ** 2).sum(1) / ((k - 1) * w_dl.sum(1))
    hk_se = np.sqrt(np.maximum(hk_var, 1e-300))
    tq1 = stats.t.ppf(0.975, k - 1)
    tq2 = stats.t.ppf(0.975, k - 2)
    I2 = np.where(Q > 0, np.maximum(0, (Q - (k - 1)) / np.where(Q > 0, Q, 1)), 0)

    rows = []
    base = dict(senaryo=idx, k=k, n_dagilimi=D, delta=delta, tau2=tau2, secici_yayin=S, R=R)
    est = {
        "FE": (mu_fe, se_fe, Z975),
        "RE-DL": (mu_dl, se_dl, Z975),
        "RE-REML": (mu_re, se_re, Z975),
        "HKSJ": (mu_dl, hk_se, tq1),
    }
    for name, (mu, se, q) in est.items():
        err = mu - delta
        bias = err.mean()
        lo, hi = mu - q * se, mu + q * se
        cov = ((lo <= delta) & (delta <= hi)).mean()
        rows.append(dict(base, grup="kestirici", yontem=name,
                         yanlilik=bias,
                         goreli_yanlilik_yuzde=(100 * bias / delta) if delta != 0 else np.nan,
                         RMSE=math.sqrt((err ** 2).mean()),
                         ampirik_SH=err.std(ddof=1),
                         MCSE_yanlilik=err.std(ddof=1) / math.sqrt(R),
                         GA_kapsama=cov,
                         MCSE_kapsama=math.sqrt(cov * (1 - cov) / R),
                         GA_genislik=(hi - lo).mean()))
    # PI kapsama (yeni çalışmanın gerçek etkisi)
    half = tq2 * np.sqrt(t2_dl + se_dl ** 2)
    lo, hi = mu_dl - half, mu_dl + half
    if tau2 > 0:
        s = math.sqrt(tau2)
        pcov = stats.norm.cdf((hi - delta) / s) - stats.norm.cdf((lo - delta) / s)
    else:
        pcov = ((lo <= delta) & (delta <= hi)).astype(float)
    err = mu_dl - delta
    rows.append(dict(base, grup="kestirici", yontem="PI (tahmin araligi)",
                     yanlilik=err.mean(),
                     goreli_yanlilik_yuzde=(100 * err.mean() / delta) if delta != 0 else np.nan,
                     RMSE=math.sqrt((err ** 2).mean()), ampirik_SH=err.std(ddof=1),
                     MCSE_yanlilik=err.std(ddof=1) / math.sqrt(R),
                     PI_kapsama=pcov.mean(), MCSE_kapsama=pcov.std(ddof=1) / math.sqrt(R),
                     PI_genislik=(hi - lo).mean()))
    # Grup 2
    pv, p = pcomb_all(y, v, n, t2_dl, mu_dl, se_dl, hk_se)
    for m in METHODS2:
        rr = float((pv[m] < ALPHA).mean())
        rows.append(dict(base, grup="p_birlestirici", yontem=m, ret_orani=rr,
                         MCSE_ret=math.sqrt(rr * (1 - rr) / R),
                         nominalden_sapma=(rr - ALPHA) if delta == 0 else np.nan))
    extra = dict(anlamsiz_oran=float((p >= ALPHA).mean()),
                 I2_ortanca=float(np.median(I2)),
                 tau2_DL_ort=float(t2_dl.mean()),
                 en_buyuk_cal_payi=float((n.max(1) / n.sum(1)).mean()))
    for r_ in rows:
        r_.update(extra)
    return rows


def _job(args):
    idx, sc = args
    return run_scenario(idx, *sc)

def main():
    from multiprocessing import Pool
    grid = list(itertools.product(K_LEV, D_LEV, DELTA_LEV, TAU2_LEV, S_LEV))
    t0 = time.time()
    allrows = []
    with Pool(2) as pool:
        for i, rows in enumerate(pool.imap(_job, list(enumerate(grid)), chunksize=1)):
            allrows += rows
            if i % 24 == 0:
                print(f"{i+1}/{len(grid)}  {time.time()-t0:.0f}s", flush=True)
    df = pd.DataFrame(allrows)
    df.to_csv(f"sim_v7_R{R}.csv", index=False)
    print("bitti", time.time() - t0)


if __name__ == "__main__":
    main()
