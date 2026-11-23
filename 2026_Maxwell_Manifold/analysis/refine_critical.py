from pathlib import Path
import json
import mpmath as mp
mp.mp.dps=65

def radial(k):
    j0=mp.sin(k)/k;j=mp.sin(k)/k**2-mp.cos(k)/k;p=j0-2*j/k
    d=p+j/k;dp=-p/k-j+j/k**2
    vm=(j*j-j0*(3*j/k-j0))/2;vmp=(j*j-3*vm)/k
    vn=vm+j*d/k;vnp=vmp+(p*d+j*dp)/k-j*d/k**2
    return j,d,p,dp,vm,vn,vmp/2,vnp/2

def coeff(pol,k,rho):
    j,d,p,dp,vm,vn,bm,bn=radial(k)
    jm=1/rho**2+1/rho**4;jn=1/rho**2-1/rho**4+1/rho**6;cross=-1/rho**5
    U,V,Up,Vp,nu,beta,aa,bb=(j,d,p,dp,vm,bm,jm,jn) if pol=='TE' else (d,j,dp,p,vn,bn,jn,jm)
    a=rho/k*U*(nu*Up-beta*U)*bb/(4*nu**2)
    c=k/rho*V*(nu*(V+k*Vp)/k-beta*V)*aa/(4*nu**2)
    b=(2*beta*U*V-nu*Up*V-nu*U*(V+k*Vp)/k)*cross/(4*nu**2)
    return a,b,c

def abc(k,rho,w):return tuple(a+w*b for a,b in zip(coeff('TM',k,rho),coeff('TE',k,rho)))
def disc(k,rho,w):a,b,c=abc(k,rho,w);return b*b-4*a*c

def refine(rho,k0,w0):
    rho=mp.mpf(str(rho))
    k,w=mp.findroot(lambda k,w:(disc(k,rho,w),mp.diff(lambda x:disc(x,rho,w),k)),(mp.mpf(str(k0)),mp.mpf(str(w0))),tol=mp.mpf('1e-58'),maxsteps=30)
    a,b,c=abc(k,rho,w);mu=-b/(2*a);eps=k*k/(rho*rho*mu)
    F=lambda ee,mm,ww: sum(v*z for v,z in zip(abc(rho*mp.sqrt(ee*mm),rho,ww),(mm,1,1/mm)))
    f=F(eps,mu,w);fe=mp.diff(lambda x:F(x,mu,w),eps);fm=mp.diff(lambda x:F(eps,x,w),mu)
    H=mp.matrix([[mp.diff(lambda x:F(x,mu,w),eps,2),mp.diff(lambda x:mp.diff(lambda y:F(x,y,w),mu),eps)],[mp.diff(lambda x:mp.diff(lambda y:F(x,y,w),mu),eps),mp.diff(lambda y:F(eps,y,w),mu,2)]])
    ev=mp.eigsy(H,eigvals_only=True)
    dkk=mp.diff(lambda x:disc(x,rho,w),k,2);dw=mp.diff(lambda x:disc(k,rho,x),w);fw=mp.diff(lambda x:F(eps,mu,x),w)
    o={'rho':rho,'k':k,'t':k*k,'w':w,'lambda':mp.sqrt(w),'epsilon':eps,'mu':mu,'F':f,'F_epsilon':fe,'F_mu':fm,'D_kk':dkk,'D_w':dw,'F_w':fw,'hessian':[[H[i,j] for j in range(2)] for i in range(2)],'hessian_eigenvalues':list(ev),'classification':'saddle' if ev[0]*ev[1]<0 else 'definite'}
    def strings(v):
      if isinstance(v,dict):return {k:strings(x) for k,x in v.items()}
      if isinstance(v,list):return [strings(x) for x in v]
      if isinstance(v,str):return v
      return mp.nstr(v,58)
    return strings(o)

if __name__=='__main__':
    cases=[(3,6.003022377754758,3.458302594955618),(3,8.449848767106866,1.0658307471841615)]
    out=[]
    for c in cases:
      r=refine(*c);out.append(r);print(json.dumps(r,indent=2),flush=True)
    Path(__file__).resolve().parents[1].joinpath('data/high_precision_critical_points.json').write_text(json.dumps(out,indent=2))
