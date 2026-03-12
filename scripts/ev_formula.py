#!/usr/bin/env python3
"""CANONICAL Ev FORMULA — import this, never redefine.
P[0,0]=0.01904482 — use as fingerprint check."""
import numpy as np
from scipy.stats import skew as skewness

PRIMES=[p for p in range(2,500) if all(p%i!=0 for i in range(2,int(p**0.5)+1))][:95]
P=np.random.default_rng(42).standard_normal((256,500))/np.sqrt(256)
Z1_T, Z3_T = 0.382, -1.471
NULL={'slope':-5.938,'intercept':4.471,'std':0.456}
BASES='ACGT'
KI={a+b+c+d:i*64+j*16+k*4+l for i,a in enumerate(BASES) for j,b in enumerate(BASES)
    for k,c in enumerate(BASES) for l,d in enumerate(BASES)}

assert abs(P[0,0]-(0.01904482))<1e-6, "WRONG RNG — use default_rng(42) not seed(42)"

def ev_resid(seg):
    seg=seg.upper().replace('N','')
    if len(seg)<1000: return None
    v=np.zeros(256)
    for i in range(len(seg)-3):
        k=seg[i:i+4]
        if k in KI: v[KI[k]]+=1
    if v.sum()==0: return None
    f=v/v.sum()
    gc=(seg.count('G')+seg.count('C'))/len(seg)
    e=abs(skewness((f@P)[PRIMES]))*6.07+0.10
    return (e-(NULL['slope']*gc+NULL['intercept']))/NULL['std']
