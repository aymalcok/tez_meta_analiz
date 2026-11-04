"""İkili veri benzetiminin tablolarını (Tablo 4-17 – 4-22) CSV olarak üretir.
Ret oranları, analizde en az üç çalışma kalan tekrarlar üzerinden hesaplanır:
ret_orani / (1 - keff_3alti_oran)."""
import pandas as pd
d = pd.read_csv('sim_bin_R2000.csv')
d['ret'] = (d.ret_orani / (1 - d.keff_3alti_oran)).round(3)
SPEC = {'Tablo_4-17': (1.0, (0.0, 0.04)), 'Tablo_4-18': (0.8, (0.0, 0.04)),
        'Tablo_4-19': (1.0, (0.01, 0.16)), 'Tablo_4-20': (0.8, (0.01, 0.16)),
        'Tablo_4-21': (0.9, (0.0, 0.04)), 'Tablo_4-22': (0.9, (0.01, 0.16))}
for ad, (rr, t2) in SPEC.items():
    s = d[(d.RR == rr) & (d.tau2.isin(t2))]
    w = s.pivot_table(index=['p_c', 'n_dagilimi', 'yontem'], columns=['tau2', 'k'], values='ret')
    w = w.sort_index(level=['p_c', 'n_dagilimi'], ascending=[False, True])
    w.to_csv(f'{ad}.csv')
    print(ad, w.shape)
