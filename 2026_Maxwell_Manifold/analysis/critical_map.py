"""Non-exhaustive numerical map of stationary discriminant-root curves.

Only the two separately certified points have rigorous existence certificates.
The finite search is not a theorem about the first possible critical size.
"""
from pathlib import Path
import csv,json
import numpy as np
from scipy.optimize import root
from modal_coefficients import coefficients
import mpmath as mp
from refine_critical import disc,abc

ROOT=Path(__file__).resolve().parents[1]

def delta_coeff(k,rho):
    a,b,c=coefficients(1,'TM',k,rho);d,e,f=coefficients(1,'TE',k,rho)
    return b*b-4*a*c,2*b*e-4*(a*f+c*d),e*e-4*d*f

def weight_roots(k,rho):
    d0,d1,d2=delta_coeff(k,rho);det=d1*d1-4*d2*d0
    with np.errstate(invalid='ignore',divide='ignore'):
        sq=np.sqrt(np.where(det>=0,det,np.nan))
        return [(-d1-sq)/(2*d2),(-d1+sq)/(2*d2)]

def dval(k,rho,w):
    a,b,c=coefficients(1,'TM',k,rho)+w*coefficients(1,'TE',k,rho)
    return b*b-4*a*c,b*b+abs(4*a*c)+1e-25

def equations(z,rho):
    k,w=z;h=2e-4
    if k<=0 or w<=0:return [100+abs(k),100+abs(w)]
    d,scale=dval(k,rho,w)
    deriv=(dval(k-2*h,rho,w)[0]-8*dval(k-h,rho,w)[0]+8*dval(k+h,rho,w)[0]-dval(k+2*h,rho,w)[0])/(12*h)
    return [d/scale,deriv/scale]

def one_rho(rho):
    ks=np.linspace(.05,14,4001);ws=weight_roots(ks,rho);found=[]
    for ww in ws:
        mask=np.isfinite(ww)&(ww>0)&(ww<1e7)
        h=2e-4
        dp=(dval(ks+h,rho,ww)[0]-dval(ks-h,rho,ww)[0])/(2*h)
        for i in np.flatnonzero(mask[:-1]&mask[1:]&(dp[:-1]*dp[1:]<0)):
            guess=[(ks[i]+ks[i+1])/2,(ww[i]+ww[i+1])/2]
            sol=root(lambda z:equations(z,rho),guess,tol=1e-10)
            k,w=sol.x
            if not (0.05<k<14 and w>0 and np.linalg.norm(equations([k,w],rho))<1e-6):continue
            if any(abs(k-p['k'])<1e-5 and abs(w-p['w'])<1e-4*max(1,w) for p in found):continue
            try:
                mp.mp.dps=42;rr=mp.mpf(str(rho))
                kk,ww0=mp.findroot(lambda x,y:(disc(x,rr,y),mp.diff(lambda u:disc(u,rr,y),x)),(mp.mpf(str(k)),mp.mpf(str(w))),tol=mp.mpf('1e-33'),maxsteps=20)
                a,b,c=abc(kk,rr,ww0);mu=-b/(2*a);ep=kk*kk/(rr*rr*mu)
                if not(.05<kk<14 and ww0>0 and mp.mpf('.1')<=abs(mu)<=100 and mp.mpf('.1')<=abs(ep)<=100):continue
                dkk=mp.diff(lambda x:disc(x,rr,ww0),kk,2);dw=mp.diff(lambda y:disc(kk,rr,y),ww0)
                if abs(dw)<mp.mpf('1e-25') or abs(dkk)<mp.mpf('1e-25'):continue
                found.append(dict(rho=float(rr),k=float(kk),w=float(ww0),epsilon=float(ep),mu=float(mu),
                                  Delta_kk=float(dkk),Delta_w=float(dw),type='saddle' if dkk>0 else 'definite',
                                  status='numerical candidate; no interval certificate'))
            except (ValueError,ZeroDivisionError):continue
    return sorted(found,key=lambda p:p['k'])

def main():
    rhos=np.arange(.5,3.0001,.125);points=[];counts=[]
    for rho in rhos:
        pp=one_rho(float(rho));points+=pp
        cc=dict(rho=float(rho),total=len(pp),saddle=sum(p['type']=='saddle' for p in pp),definite=sum(p['type']=='definite' for p in pp))
        counts.append(cc);print(cc,flush=True)
    for name,rows in [('fig8_critical_map',points),('fig8_candidate_counts',counts)]:
        with (ROOT/'data'/f'{name}.csv').open('w',newline='') as f:
            wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
    curve=[]
    for rho in [1.,1.5,2.,3.]:
        ks=np.linspace(.05,14,4001)
        for branch,ws in enumerate(weight_roots(ks,rho)):
            for k,w in zip(ks,ws):
                if np.isfinite(w) and 0<w<200:curve.append(dict(rho=rho,k=float(k),w=float(w),branch=branch))
    with (ROOT/'data/fig8_discriminant_root_curves.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(curve[0]));wr.writeheader();wr.writerows(curve)
    (ROOT/'audit/critical_map.json').write_text(json.dumps(dict(search='non-exhaustive; t>0, .05<k<14, .1<=abs(epsilon),abs(mu)<=100',counts=counts,candidates=points),indent=2))

if __name__=='__main__':main()
