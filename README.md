# Tez benzetim kodları

Aydan Malçok Demirtaş. *Çalışmaların Örnek Büyüklükleri ve p Değerlerinin Meta Analiz Sonuçlarına Etkileri ve Karşılaşılan Sorunlara İlişkin Yeni Çözüm Önerisi.* Doktora tezi, İstanbul Üniversitesi Sağlık Bilimleri Enstitüsü, Biyoistatistik ve Tıp Bilişimi Anabilim Dalı.

Bu depo, tezdeki benzetim çalışmalarının kodlarını, ham sonuçlarını ve şekil üreten betikleri içermektedir. MAT yönteminin uygulama kodu ayrı bir depoda (github.com/aymalcok/MAT) yer almaktadır.

## Dosyalar

| Dosya | İçerik | Tezdeki karşılığı |
|---|---|---|
| `hypo.py` | Hipoüstel dağılımın tam sağ kuyruk olasılığı (matris üstel fonksiyonu) | Bölüm 2.8; GRPC, hibrit, MAT ve RTP p değerleri |
| `sim_v7.py` | Ana benzetim (324 senaryo, 2.000 tekrar, tohum 20260924 + senaryo sırası) | Bölüm 3.5.1–3.5.5; Tablo 4-7 – 4-14 |
| `hib.py` | Ana benzetimin aynı verisinde hibrit GRPC × Shannon yöntemi | Tablo 4-7 – 4-13'teki hibrit satırı |
| `sim_bin.py` | İkili veriyle (2×2 tablo) duyarlılık analizi (360 senaryo, 2.000 tekrar, tohum 20261020 + senaryo sırası) | Bölüm 3.5.6; Tablo 4-17 – 4-22 |
| `tablolar_ikili.py` | İkili veri benzetiminin tablolarını CSV olarak üretir | Tablo 4-17 – 4-22 |
| `sim_anlamsiz.py` | Anlamsız çalışma oranlarına ilişkin ek benzetim (324 senaryo, 20.000 tekrar, tek üreteç, tohum 7) | Bölüm 3.5.5; Tablo 4-15 ve 4-16 |
| `tablolar_anlamsiz.py` | Ek benzetimin tablolarını CSV olarak üretir | Tablo 4-15 ve 4-16 |
| `sim_kucuk.py` | Küçük örneklemli çalışmalarla duyarlılık analizi (40 senaryo, 2.000 tekrar, tohum 20270000 + senaryo sırası) | Bölüm 3.5.7; Tablo 4-23 |
| `sekiller.py`, `sekiller_guc.py`, `_ortak.py` | Bulgular bölümündeki benzetim şekilleri | Şekil 4-1 – 4-8 |
| `compare.py`, `compare_bin.py` | Sonuçları tezdeki tablolarla hücre hücre karşılaştıran kontrol betikleri | — |
| `*_R2000.csv` | Ham sonuçlar | — |

## Çalıştırma

```
pip install -r requirements.txt
python3 sim_v7.py 2000                       # sim_v7_R2000.csv
python3 hib.py                               # hibrit_R2000.csv (sim_v7.py ve hypo.py gerekir)
python3 sim_bin.py 2000 sim_bin_R2000.csv 2  # sim_bin_R2000.csv
python3 tablolar_ikili.py                    # Tablo_4-17.csv ... Tablo_4-22.csv
python3 sim_anlamsiz.py                      # sim_anlamsiz_R20000.csv (tek süreçte çalışır)
python3 tablolar_anlamsiz.py                 # Tablo_4-15.csv, Tablo_4-16.csv
python3 sim_kucuk.py 2000                    # sim_kucuk_R2000.csv
python3 sekiller.py; python3 sekiller_guc.py # sekiller/ klasörü
```

## Notlar

- İkili veri benzetiminde ret oranları, analizde en az üç çalışma kalan tekrarlar üzerinden hesaplanır: `ret_orani / (1 - keff_3alti_oran)`. `sim_bin_R2000.csv` dosyasındaki `ret_orani` sütunu bütün tekrarlar üzerinden hesaplanmıştır; `tablolar_ikili.py` bu düzeltmeyi uygular.
- Ana benzetim ve ikili veri benzetiminin tohum değeri aralıkları kısmen örtüşmektedir (20261020–20261247). İki benzetimin veri üretim mekanizmaları farklıdır ve sonuçları birlikte bir çıkarımda kullanılmamaktadır.
- Rastgele sayı üretimi NumPy sürümüne bağlıdır. Sonuçlar aşağıdaki sürümlerle üretilmiştir.
- Anlamsız çalışma oranlarına ilişkin ek benzetimde bütün senaryolar tek bir rastgele sayı üretecinden (tohum 7), `sim_v7.py` ile aynı senaryo sırasında ardışık olarak üretilir. Bu nedenle betik paralel çalıştırılmamalıdır. Tablo 4-16'nın B bölümünde ret oranı, yalnızca anlamsız çalışmalardan oluşan meta analizlerin sayısı (m) 100'ün altında kaldığında verilmez. m değerleri `Tablo_4-16.csv` dosyasının son satırlarında yer almaktadır.
- v1.1.0: Tablo 4-15 ve 4-16'yı üreten ek benzetim ile küçük örneklemli duyarlılık analizi (Tablo 4-23) pakete eklenmiştir. Tablo 4-15 ve 4-16 tezde bu sürümdeki kodun çıktısıyla verilmektedir.
- Küçük örneklemli duyarlılık analizinin tohum değerleri (20270000–20270039) diğer benzetimlerin tohum değerleriyle örtüşmemektedir.

## Yazılım sürümleri

Python 3.13, NumPy 2.5.3, SciPy 1.18.1, pandas 3.0.5, matplotlib.
