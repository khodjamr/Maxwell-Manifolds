"""Exact interval arithmetic shared by the three proof replayers.

This module uses only the standard library. All arithmetic endpoints remain
fractions; decimal strings are for display, never interval-proof inputs.
The formula definitions are unchanged from the supplied certificate engine.
"""
from fractions import Fraction as Q
from math import factorial
from decimal import Decimal, localcontext

if not __debug__:
    raise RuntimeError("Run the certificate checker without Python -O or -OO.")

BITS=360
GRID=1<<BITS
N=90

def decimal(value):
    with localcontext() as context:
        context.prec=65
        return str(Decimal(value.numerator)/Decimal(value.denominator))

def down(x):return Q((x.numerator*GRID)//x.denominator,GRID)

def up(x):return -down(-x)

def frac(x):
 if isinstance(x,Q):return x
 if isinstance(x,int):return Q(x)
 return Q(str(x))

class I:
 def __init__(self,lo=0,hi=None):
  if isinstance(lo,I) and hi is None:self.lo,self.hi=lo.lo,lo.hi;return
  self.lo=down(frac(lo));self.hi=up(frac(lo if hi is None else hi));assert self.lo<=self.hi
 def __add__(self,b):b=I(b);return I(self.lo+b.lo,self.hi+b.hi)
 __radd__=__add__
 def __neg__(self):return I(-self.hi,-self.lo)
 def __sub__(self,b):return self+-I(b)
 def __rsub__(self,b):return I(b)+-self
 def __mul__(self,b):
  if isinstance(b,J):return b*self
  b=I(b);v=[self.lo*b.lo,self.lo*b.hi,self.hi*b.lo,self.hi*b.hi];return I(min(v),max(v))
 __rmul__=__mul__
 def inv(self):
  assert not self.lo<=0<=self.hi, 'interval division by zero'
  return I(1/self.hi,1/self.lo)
 def __truediv__(self,b):return self*I(b).inv()
 def __rtruediv__(self,b):return I(b)*self.inv()
 def __pow__(self,n):
  if n<0:return self.inv()**(-n)
  out=I(1)
  for _ in range(n):out=out*self
  return out
 def absmax(self):return max(abs(self.lo),abs(self.hi))
 def bounds(self):return [decimal(self.lo),decimal(self.hi)]

def trig(iv,cosine=False):
 c=(iv.lo+iv.hi)/2;rad=(iv.hi-iv.lo)/2;v=I(c);v2=v*v
 term=I(1) if cosine else v;total=term
 for n in range(1,N+1):
  term= -term*v2/((2*n-1)*(2*n) if cosine else (2*n)*(2*n+1));total=total+term
 degree=2*N if cosine else 2*N+1
 rem=abs(c)**(degree+1)/factorial(degree+1)
 return total+I(-rem-rad,rem+rad)

class J:
 """Formal Taylor coefficients in k, through order three."""
 def __init__(self,v=0):
  if isinstance(v,J):self.c=v.c;return
  self.c=[I(v),I(0),I(0),I(0)]
 @classmethod
 def make(cls,c):z=cls();z.c=[I(x) for x in c];return z
 def __add__(self,b):b=J(b);return J.make([a+x for a,x in zip(self.c,b.c)])
 __radd__=__add__
 def __neg__(self):return J.make([-x for x in self.c])
 def __sub__(self,b):return self+-J(b)
 def __rsub__(self,b):return J(b)+-self
 def __mul__(self,b):b=J(b);return J.make([sum((self.c[j]*b.c[n-j] for j in range(n+1)),I(0)) for n in range(4)])
 __rmul__=__mul__
 def inv(self):
  out=[self.c[0].inv()]
  for n in range(1,4):out.append(-sum((self.c[j]*out[n-j] for j in range(1,n+1)),I(0))/self.c[0])
  return J.make(out)
 def __truediv__(self,b):return self*J(b).inv()
 def __rtruediv__(self,b):return J(b)*self.inv()
 def __pow__(self,n):
  if n<0:return self.inv()**(-n)
  z=J(1)
  for _ in range(n):z=z*self
  return z
 def deriv(self):return J.make([self.c[1],2*self.c[2],3*self.c[3],0])
 def sincos(self):
  s=trig(self.c[0]);c=trig(self.c[0],True);h=J.make([0]+self.c[1:]);
  return (J(s)+h*c-h*h*s/2-h*h*h*c/6,J(c)-h*s-h*h*c/2+h*h*h*s/6)

def coeff_jets(kbox,rho):
 k=J.make([kbox,1,0,0]);s,c=k.sincos();j0=s/k;j=s/k**2-c/k;p=j0-2*j/k;d=p+j/k
 vm=(j*j-j0*(3*j/k-j0))/2;vn=vm+j*d/k
 aa=Q(1,rho**2)+Q(1,rho**4);bb=Q(1,rho**2)-Q(1,rho**4)+Q(1,rho**6);cross=-Q(1,rho**5)
 result=[]
 for pol in ['TM','TE']:
  U,V,nu,exta,extb=(d,j,vn,bb,aa) if pol=='TM' else (j,d,vm,aa,bb)
  r11=U*U/nu;r22=k*k*V*V/nu;r12=k*U*V/nu
  A=rho*extb*r11.deriv()/(8*k);B=-cross*r12.deriv()/(4*k);C=exta*r22.deriv()/(8*rho*k)
  result.append((A,B,C))
 return result

def interval_system(kbox,wbox,rho):
 (ae,be,ce),(am,bm,cm)=coeff_jets(kbox,rho)
 A=ae+am*wbox;B=be+bm*wbox;C=ce+cm*wbox
 D=B*B-4*A*C;Dw=2*B*bm-4*(am*C+A*cm)
 return [D.c[0],D.c[1]],[[D.c[1],Dw.c[0]],[2*D.c[2],Dw.c[1]]],(A.c[0],B.c[0],C.c[0]),(2*D.c[2],Dw.c[0])
