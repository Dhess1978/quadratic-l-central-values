# Úplnost: množina fundamentálních diskriminantů obou znamének s |D| v [lo, hi) proti datům; duplicity; rozsahy bloků.
import json, gzip, os, sys, numpy as np
base=sys.argv[1]; m=json.load(open(os.path.join(base,'manifest.json')))
def sqfree(lo,hi):
    a=np.ones(hi-lo,bool); r=int(hi**.5)+2
    s=np.ones(r+1,bool); s[:2]=False
    for p in range(2,int(r**.5)+1):
        if s[p]: s[p*p::p]=False
    for p in np.nonzero(s)[0]:
        q=int(p)*int(p); st=(-lo)%q; a[st::q]=False
    return a
def fund(lo,hi):
    out=set()
    sf=sqfree(lo//4-2, hi//4+2); base4=lo//4-2
    sfa=sqfree(lo,hi)
    for sgn in (1,-1):
        for x in range(lo,hi):
            D=sgn*x
            if D%4==1 and sfa[x-lo]: out.add(D)
            elif D%4==0:
                mm=D//4
                if mm%4 in (2,3) and sf[abs(mm)-base4]: out.add(D)
    return out
rep={}
for fam,v in m['rodiny'].items():
    lo,hi=v['spec']['lo'],v['spec']['hi']
    seen=[]; nL_man=0; mimo=0
    for k,b in v['bloky'].items():
        nL_man+=b['nL']
        for l in gzip.open(os.path.join(base,fam,b['soubor']),'rt'):
            r=json.loads(l); seen.append(r['D'])
    S=set(seen); F=fund(lo,hi)
    rep[fam]=dict(lo=lo,hi=hi,zaznamu=len(seen),unikatnich=len(S),duplicit=len(seen)-len(S),
                  ocekavano=len(F),chybi=len(F-S),navic=len(S-F),nL_manifest=nL_man,
                  min_absD=min(map(abs,S)),max_absD=max(map(abs,S)),kladnych=sum(1 for d in S if d>0))
    print(fam,rep[fam],flush=True)
json.dump(rep,open('/home/claude/br78/data/uplnost_D100bc.json','w'),indent=1)
