"""sim_anlamsiz_R20000.csv dosyasından Tablo 4-15 ve Tablo 4-16'yı CSV olarak üretir."""
import pandas as pd
S = pd.read_csv("sim_anlamsiz_R20000.csv")
K = [3, 4, 5, 7, 10, 15, 20, 30, 50]
YON = {"Fisher": "Fisher", "Stouffer": "Stouffer", "n-Stouffer": "Ağırlıklı Z (√n)",
       "MAT": "MAT", "PI": "PI", "HKSJ": "HKSJ"}

def f(v, d):
    if v == 0:
        return "0"
    if v < 0.5 * 10 ** -d:
        return "<0," + "0" * (d - 1) + "1"
    return f"{v:.{d}f}".replace(".", ",")

# Tablo 4-15: anlamsız çalışma bulunmayan meta analiz oranı (tau2 üzerinden ortalama)
g = S.groupby(["delta", "n_dagilimi", "k"]).hic_anlamsiz_yok.mean()
rows = [[d, D] + [f(g[(d, D, k)], 4) for k in K] for d in (0.0, 0.1, 0.2) for D in ("D1", "D2", "D3")]
pd.DataFrame(rows, columns=["delta", "dagilim"] + K).to_csv("Tablo_4-15.csv", index=False)

# Tablo 4-16: A oran (D ve tau2 üzerinden ortalama); B ret oranı (havuzlanmış, m < 100 ise "–")
g = S.groupby(["delta", "k"])
A, N = g.hepsi_anlamsiz.mean(), g.n_hepsi_anlamsiz.sum()
rows = [["A", d, "oran"] + [f(A[(d, k)], 3) for k in K] for d in (0.0, 0.1, 0.2)]
for d in (0.0, 0.1, 0.2):
    for m, ad in YON.items():
        r = g["ret_" + m].sum()
        rows.append(["B", d, ad] + [f(r[(d, k)] / N[(d, k)], 3) if N[(d, k)] >= 100 else "–" for k in K])
rows += [["m", d, "meta analiz sayısı"] + [int(N[(d, k)]) for k in K] for d in (0.0, 0.1, 0.2)]
pd.DataFrame(rows, columns=["bolum", "delta", "satir"] + K).to_csv("Tablo_4-16.csv", index=False)
print("Tablo_4-15.csv ve Tablo_4-16.csv yazıldı")
