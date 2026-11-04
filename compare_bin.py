import docx, re, pandas as pd
from decimal import Decimal, ROUND_HALF_UP
f3=lambda x: str(Decimal(repr(round(x,6))).quantize(Decimal('0.001'),rounding=ROUND_HALF_UP)).replace('.',',')
W='{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
def TX(e): return ''.join(t.text or '' for t in e.iter(W+'t'))
d=docx.Document('/mnt/user-data/outputs/tez/Tez_Tam_Metin_v34.docx'); B=list(d.element.body.iterchildren())
def table_after(cap):
    for i,e in enumerate(B):
        if e.tag==W+'p' and TX(e).startswith(cap):
            j=i+1
            while B[j].tag!=W+'tbl': j+=1
            return B[j]
X=pd.read_csv('sim_bin_R2000.csv'); X['ret_c']=X.ret_orani/(1-X.keff_3alti_oran)
LAB={'Fisher':'Fisher','Stouffer':'Stouffer','Ağırlıklı Z (√n)':'n-Stouffer','GRPC':'GRPC','Hibrit GRPC × Shannon':'Hibrit','MAT':'MAT','PI':'PI','HKSJ':'HKSJ','I² ortancası':'I2'}
K5=[3,5,10,20,50]
def val(m,**c):
    col='ret_c'
    if m=='I2': m='Fisher'; col='I2_ortanca'
    s=X[X.yontem==m]
    for a,b in c.items(): s=s[s[a]==b]
    return float(s[col].iloc[0])
tot=bad=0
for cap,RR in [('Tablo 4.15.',1.0),('Tablo 4.16.',0.8)]:
    tb=table_after(cap); pc=D=None
    for r in tb.findall(W+'tr'):
        cells=[TX(c).strip() for c in r.findall(W+'tc')]
        if len(cells)==1:
            c0=cells[0].replace('₁','1').replace('₂','2').replace('₃','3')
            m=re.search(r'riski = (0,\d+), D([123])',c0)
            if m: pc=float(m.group(1).replace(',','.')); D='D'+m.group(2)
            continue
        if cells[0] not in LAB or D is None: continue
        exp=[f3(val(LAB[cells[0]],RR=RR,p_c=pc,n_dagilimi=D,tau2=t,k=k)) for t in [0,0.04] for k in K5]
        for g,e in zip(cells[1:],exp):
            tot+=1
            if g!=e: bad+=1; print(cap,pc,D,cells[0],g,e)
print('binary tables: compared',tot,'mismatch',bad)
