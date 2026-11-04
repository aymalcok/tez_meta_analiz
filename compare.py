import docx, re, pandas as pd, sys
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
S=pd.read_csv('sim_v7_R2000.csv'); S=S[S.grup=='p_birlestirici'][['k','n_dagilimi','delta','tau2','yontem','ret_orani']]
H=pd.read_csv('hibrit_R2000.csv').rename(columns={'Hibrit':'ret_orani'}); H['yontem']='Hibrit'
S=pd.concat([S,H[['k','n_dagilimi','delta','tau2','yontem','ret_orani']]])
LAB={'Fisher':'Fisher','Stouffer':'Stouffer','Tippett':'Tippett','Ağırlıklı Z (√n)':'n-Stouffer','Lancaster':'Lancaster','wFisher':'wFisher','TPM':'TPM','RTP':'RTP','ordmeta':'ordmeta','GRPC':'GRPC','Hibrit GRPC × Shannon':'Hibrit','MAT':'MAT','PI':'PI','HKSJ':'HKSJ'}
def val(m,**c):
    s=S[S.yontem==m]
    for a,b in c.items(): s=s[s[a]==b]
    return float(s.ret_orani.iloc[0])
KS=[3,4,5,7,10,15,20,30,50]; K5=[3,5,10,20,50]
specs={'Tablo 4.7.':('k9',0,[0]),'Tablo 4.8.':('k9',0,[0.04]),'Tablo 4.9.':('k9',0,[0.01,0.16]),'Tablo 4.10.':('k5',0.1,[0,0.04]),'Tablo 4.11.':('k5',0.2,[0,0.04])}
tot=bad=0
for cap,(kind,dl,taus) in specs.items():
    tb=table_after(cap); D=None; panel=-1
    for r in tb.findall(W+'tr'):
        cells=[TX(c).strip() for c in r.findall(W+'tc')]
        if len(cells)==1:
            m=re.search(r'D([123])',cells[0].replace('₁','1').replace('₂','2').replace('₃','3'))
            if m and cells[0][:2] in ('A.','B.','C.','D.','E.','F.'): D='D'+m.group(1); panel+=1
            continue
        lab=cells[0]
        if lab not in LAB or D is None: continue
        if kind=='k9':
            t=taus[0] if len(taus)==1 else taus[panel//3]
            exp=[f3(val(LAB[lab],delta=dl,tau2=t,n_dagilimi=D,k=k)) for k in KS]
        else:
            exp=[f3(val(LAB[lab],delta=dl,tau2=t,n_dagilimi=D,k=k)) for t in taus for k in K5]
        got=cells[1:]
        for g,e in zip(got,exp):
            tot+=1
            if g!=e: bad+=1; print(cap,D,lab,g,e)
print('main tables: compared',tot,'mismatch',bad)
