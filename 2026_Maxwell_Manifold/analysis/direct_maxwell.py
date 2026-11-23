"""Independent direct Maxwell channel including all four real material quadrants."""
import numpy as np
from scipy.special import spherical_jn,spherical_yn
from numpy.polynomial.legendre import leggauss
from functools import lru_cache
@lru_cache(maxsize=None)
def quad(n):
 x,w=leggauss(n);r=(1+x)/2;return r,w*r*r/2

def waves(l,k,r,outgoing=False):
 z=k*np.asarray(r);j=spherical_jn(l,z);p=spherical_jn(l,z,True)
 if outgoing:j=j+1j*spherical_yn(l,z);p=p+1j*spherical_yn(l,z,True)
 pp=-2*p/z-(1-l*(l+1)/z**2)*j;d=p+j/z;dp=pp+p/z-j/z**2
 return [np.array([j]),np.array([np.sqrt(l*(l+1))*j/z,d]),np.array([np.asarray(r)*p]),np.array([np.sqrt(l*(l+1))*np.asarray(r)*(p/z-j/z**2),np.asarray(r)*dp])]

def channel(e,m,rho,l=1,pol='TM',n=120):
 assert e!=0 and m!=0 and rho>0
 k=rho*np.sqrt(complex(e*m));k=float(k.real) if abs(k.imag)==0 else k
 ui,vi=(0,1) if pol=='TE' else (1,0);b=waves(l,k,1.);out=waves(l,rho,1.,True)
 U,V,Up,Vp=b[ui][-1],b[vi][-1],b[ui+2][-1],b[vi+2][-1]
 mat=np.array([[U,-out[ui][-1]],[k*V/(1j*rho*m),1j*out[vi][-1]]])
 h,B=np.linalg.solve(mat,[-1j*rho*m/(2*k)*Up,-(V+k*Vp)/(2*k)])
 rr,ww=quad(n);a=waves(l,k,rr);u,uk=a[ui],a[ui+2];E=1j*rho*m/(2*k)*uk+h*u
 inner=lambda x,y:np.sum(ww*np.sum(np.conj(x)*y,axis=0))
 nu=inner(u,u).real;P=-inner(u,E)/2
 return {'q':P.imag/abs(B)**2,'resource':nu/(2*abs(B)**2),'power_balance':P.real/(abs(B)**2/(2*rho*rho)),'B':B,'nu':nu}
