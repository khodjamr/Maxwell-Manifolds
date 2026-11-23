"""Check the predicted O(eta^2) shift after adding an exact TM quadrupole."""
from pathlib import Path
import json,csv
import mpmath as mp
from refine_critical import abc
mp.mp.dps=55
ROOT=Path(__file__).resolve().parents[1]

def quadrupole(k,rho):
    l=2
    def j(n,z):return mp.sqrt(mp.pi/(2*z))*mp.besselj(n+mp.mpf('.5'),z)
    def y(n,z):return mp.sqrt(mp.pi/(2*z))*mp.bessely(n+mp.mpf('.5'),z)
    J=j(l,k);jp=j(l-1,k)-(l+1)*J/k
    jpp=-2*jp/k-(1-l*(l+1)/k**2)*J
    d=jp+J/k;dp=jpp+jp/k-J/k**2
    vm=(J**2-j(l-1,k)*j(l+1,k))/2;vmp=(J**2-3*vm)/k
    vn=vm+J*d/k;vnp=vmp+(jp*d+J*dp)/k-J*d/k**2
    h=j(l,rho)+1j*y(l,rho)
    dh=j(l-1,rho)+1j*y(l-1,rho)-l*h/rho
    U,V,Up,Vp,nu,beta=d,J,dp,jp,vn,vnp/2
    aa,bb,cross=abs(dh)**2,abs(h)**2,mp.re(dh*mp.conj(h))
    return (rho*U*(nu*Up-beta*U)*bb/(4*k*nu**2),
            (2*beta*U*V-nu*Up*V-nu*U*(V+k*Vp)/k)*cross/(4*nu**2),
            k*V*(nu*(V+k*Vp)/k-beta*V)*aa/(4*rho*nu**2))

def run(case):
    k,w,rho=[mp.mpf(case[s]) for s in ['k','w','rho']]
    def coeff(x,y,p):return tuple(a+p*b for a,b in zip(abc(x,rho,y),quadrupole(x,rho)))
    def D(x,y,p):a,b,c=coeff(x,y,p);return b*b-4*a*c
    G=lambda x,y,p:mp.matrix([D(x,y,p),mp.diff(lambda z:D(z,y,p),x)])
    Jac=mp.matrix([[mp.diff(lambda z:G(z,w,0)[i],k),mp.diff(lambda z:G(k,z,0)[i],w)] for i in range(2)])
    gp=mp.matrix([mp.diff(lambda p:G(k,w,p)[i],0) for i in range(2)])
    deriv=mp.lu_solve(Jac,-gp);rows=[]
    for eta in [mp.mpf('0.01'),mp.mpf('0.005'),mp.mpf('0.0025'),mp.mpf('0.00125')]:
        p=eta**2
        kn,wn=mp.findroot(lambda x,y:tuple(G(x,y,p)),(k+deriv[0]*p,w+deriv[1]*p),tol=mp.mpf('1e-43'))
        err=mp.sqrt((kn-k-deriv[0]*p)**2+(wn-w-deriv[1]*p)**2)
        rows.append([case['classification'],str(eta),mp.nstr(kn-k,35),mp.nstr(wn-w,35),mp.nstr(deriv[0]*p,35),mp.nstr(deriv[1]*p,35),mp.nstr(err,35)])
    return rows

if __name__=='__main__':
    cases=json.loads((ROOT/'data/high_precision_critical_points.json').read_text())
    rows=sum([run(c) for c in cases],[])
    with (ROOT/'data/higher_multipole_persistence.csv').open('w',newline='') as f:
        wr=csv.writer(f);wr.writerow(['type','quadrupole_amplitude_eta','delta_k','delta_w','predicted_delta_k','predicted_delta_w','joint_linearized_shift_error']);wr.writerows(rows)
    print(json.dumps(rows,indent=2))
