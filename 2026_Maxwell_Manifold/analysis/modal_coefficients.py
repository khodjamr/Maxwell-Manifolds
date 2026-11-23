"""Exact spherical reactive-power reduction, a=epsilon0=mu0=1.

Fixed outgoing coefficients; real nonzero material parameters.
F(k,mu)=a(k) mu+b(k)+c(k)/mu, k=rho sqrt(epsilon mu).
Analytic real-axis continuation is used for derivative checks only.
"""
import numpy as np
from scipy.special import spherical_jn as jn,spherical_yn as yn

def radial(l,k):
    j=jn(l,k);p=jn(l,k,True);L=l*(l+1)
    pp=-2*p/k-(1-L/k**2)*j
    d=p+j/k;dp=pp+p/k-j/k**2
    vm=.5*(j*j-jn(l-1,k)*jn(l+1,k))
    vmp=(j*j-3*vm)/k
    vn=vm+j*d/k
    vnp=vmp+(p*d+j*dp)/k-j*d/k**2
    return j,d,p,dp,vm,vn,vmp/2,vnp/2

def coefficients(l,pol,k,rho):
    j,d,p,dp,vm,vn,bm,bn=radial(l,k)
    je=jn(l,rho);ye=yn(l,rho)
    jde=jn(l,rho,True)+je/rho;yde=yn(l,rho,True)+ye/rho
    if pol=='TE':U,V,Up,Vp,nu,beta=j,d,p,dp,vm,bm;aa=je*je+ye*ye;bb=jde*jde+yde*yde
    else:U,V,Up,Vp,nu,beta=d,j,dp,p,vn,bn;aa=jde*jde+yde*yde;bb=je*je+ye*ye
    cross=je*jde+ye*yde
    a=rho/k*U*(nu*Up-beta*U)*bb/(4*nu**2)
    c=k/rho*V*(nu*(V+k*Vp)/k-beta*V)*aa/(4*nu**2)
    b=(2*beta*U*V-nu*Up*V-nu*U*(V+k*Vp)/k)*cross/(4*nu**2)
    return np.array([a,b,c])

def mixed_coeff(k,rho,w):
    return coefficients(1,'TM',k,rho)+w*coefficients(1,'TE',k,rho)

def power(e,m,rho,w):
    product=np.asarray(e*m)
    k=rho*np.sqrt(product if np.all(product>0) else product.astype(complex))
    a,b,c=mixed_coeff(k,rho,w)
    return np.real(a*m+b+c/m)

def discriminant(k,rho,w):
    a,b,c=mixed_coeff(k,rho,w)
    return b*b-4*a*c

if __name__=='__main__':
    from direct_maxwell import channel
    for e,m,r in [(2,3,.3),(-3,-1.7,.4),(-4,-3,.12),(4,2,2.)]:
        k=r*np.sqrt(e*m)
        for l in [1,2,3]:
            for pol in ['TE','TM']:
                abc=coefficients(l,pol,k,r);qq=abc@[m,1,1/m];ref=channel(e,m,r,l,pol,n=150)['q']
                print(e,m,r,l,pol,float(qq),float(ref),float(abs(qq-ref)/max(1,abs(ref))))
    for r in [.1,.4,1.,2.,3.,5.]:
        ks=np.linspace(.15,25,1500)
        for w in [.01,.1,1.,1.69,10,100.]:
            vals=np.array([mixed_coeff(k,r,w) for k in ks])
            D=vals[:,1]**2-4*vals[:,0]*vals[:,2]
            dn=D/(vals[:,1]**2+abs(4*vals[:,0]*vals[:,2])+1e-99)
            print('SCAN',r,w,'min signed normalized D',float(min(dn)),'negative count',int(np.sum(D<0)))
