"""Independent volume-work Maxwell evaluation at arbitrary precision.

No material derivative or discriminant formula is used. A direct transmission
system is solved and the interior E.J* integral is evaluated by quadrature.
"""
import mpmath as mp
from direct_maxwell import quad

def channel(e,m,rho,l=1,pol='TM',n=64,dps=55):
    with mp.workdps(dps):
        e,m,rho=map(lambda x:mp.mpf(str(x)),[e,m,rho]);k=rho*mp.sqrt(e*m)
        def sj(ell,z):return mp.sqrt(mp.pi/(2*z))*mp.besselj(ell+mp.mpf('.5'),z)
        def sy(ell,z):return mp.sqrt(mp.pi/(2*z))*mp.bessely(ell+mp.mpf('.5'),z)
        def waves(q,r,out=False):
            z=q*r;j=sj(l,z);jm=sj(l-1,z)
            if out:j+=1j*sy(l,z);jm+=1j*sy(l-1,z)
            p=jm-(l+1)*j/z;pp=-2*p/z-(1-l*(l+1)/z**2)*j
            d=p+j/z;dp=pp+p/z-j/z**2
            return [[j],[mp.sqrt(l*(l+1))*j/z,d],[r*p],
                    [mp.sqrt(l*(l+1))*r*(p/z-j/z**2),r*dp]]
        ui,vi=(0,1) if pol=='TE' else (1,0)
        boundary=waves(k,mp.mpf(1));out=waves(rho,mp.mpf(1),True)
        U,V,Up,Vp=boundary[ui][-1],boundary[vi][-1],boundary[ui+2][-1],boundary[vi+2][-1]
        mat=mp.matrix([[U,-out[ui][-1]],[k*V/(1j*rho*m),1j*out[vi][-1]]])
        h,B=mp.lu_solve(mat,mp.matrix([-1j*rho*m*Up/(2*k),-(V+k*Vp)/(2*k)]))
        rs,ws=quad(n);nu=mp.mpf(0);work=mp.mpc(0)
        for rv,wv in zip(rs,ws):
            rv,wv=mp.mpf(str(rv)),mp.mpf(str(wv));v=waves(k,rv);u,uk=v[ui],v[ui+2]
            E=[1j*rho*m*x/(2*k)+h*y for x,y in zip(uk,u)]
            nu+=wv*sum(abs(x)**2 for x in u)
            work+=wv*sum(x*mp.conj(y) for x,y in zip(E,u))
        P=-work/2
        return {'q':float(mp.im(P)/abs(B)**2),'resource':float(nu/(2*abs(B)**2)),
                'power_balance':float(mp.re(P)/(abs(B)**2/(2*rho*rho)))}
