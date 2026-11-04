"""Hipoüstel (ağırlıklı üstel toplam) dağılımın tam sağ kuyruk olasılığı.
T = sum w_i E_i, E_i ~ Exp(1) bağımsız. Faz tipi gösterim (Neuts 1981):
P(T >= t) = e1' expm(A t) 1, A: köşegen -1/w_i, üst köşegen 1/w_i.
Eşit ağırlıklar dahil her durumda tam; expm = Al-Mohy & Higham (2009)."""
import numpy as np
from scipy.linalg import expm

REL = 1e-12  # en büyük ağırlığın 1e-12'sinden küçük ağırlıklar ihmal edilir

def hypo_sf(t, w):
    w = np.asarray(w, float)
    if w.size == 0 or w.max() <= 0:
        return 1.0
    w = w[w > REL * w.max()]
    if t <= 0:
        return 1.0
    lam = 1.0 / w
    A = np.diag(-lam) + np.diag(lam[:-1], 1)
    s = float(expm(A * t)[0].sum())
    return min(max(s, 0.0), 1.0)
