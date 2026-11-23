"""Weak-loss Maxwell dipoles: conjugate-profile minimum-norm currents.

The real-parameter formula is not analytically continued as a power formula.
Hermitian norms and the conjugated interior radial source are retained. The
reactive power is independently checked against Cauchy-matched volume work.
Finite-loss critical points are high-precision numerical continuations, not
new interval-certified roots. Rigorous persistence follows locally by IFT.
"""
from pathlib import Path
from functools import lru_cache
import csv,json
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
mp.mp.dps=55

@lru_cache(None)
def quad(n=40):
    with mp.workdps(70):
        x,w=mp.gauss_quadrature(n,'legendre')
        return tuple((mp.mpf((x[i]+1)/2),mp.mpf(w[i]/2)*((x[i]+1)/2)**2) for i in range(n))

def profiles(t,r,pol):
    k=mp.sqrt(t);z=k*r;j0=mp.sin(z)/z;j=mp.sin(z)/z**2-mp.cos(z)/z
    jp=j0-2*j/z;jpp=-2*jp/z-(1-2/z**2)*j
    d=jp+j/z;dp=jpp+jp/z-j/z**2
    if pol=='TE':
        u=[j/k];ut=[(r*jp-j/k)/(2*k*k)]
        cu=[mp.sqrt(2)*j/z,d]
        cut=[mp.sqrt(2)*r*(jp/z-j/z**2)/(2*k),r*dp/(2*k)]
    else:
        u=[mp.sqrt(2)*j/z,d]
        ut=[mp.sqrt(2)*r*(jp/z-j/z**2)/(2*k),r*dp/(2*k)]
        cu=[k*j];cut=[(j+k*r*jp)/(2*k)]
    return u,ut,cu,cut

def boundary(t,pol):
    u,ut,cu,cut=profiles(t,mp.mpf(1),pol)
    return u[-1],ut[-1],cu[-1],cut[-1]

def exterior(rho,pol):
    h0=(mp.sin(rho)-1j*mp.cos(rho))/rho
    h1=mp.sin(rho)/rho**2-mp.cos(rho)/rho+1j*(-mp.cos(rho)/rho**2-mp.sin(rho)/rho)
    dh=h0-h1/rho
    return (h1,dh) if pol=='TE' else (dh,h1)

def channel(e,m,rho,pol,n=40,direct=False):
    t=rho*rho*e*m;U,Ut,W,Wt=boundary(t,pol);Ae,Be=exterior(rho,pol)
    nu=mp.mpf(0);moment=mp.mpc(0)
    samples=[]
    for r,weight in quad(n):
        u,ut,cu,cut=profiles(t,r,pol)
        nu+=weight*sum(abs(x)**2 for x in u)
        moment+=weight*sum(x*mp.conj(y) for x,y in zip(ut,u))
        if direct:samples.append((r,weight,u,ut,cu,cut))
    Z=U*Be-W*Ae/(rho*m);Ze=rho*rho*m*Ut*Be-rho*Wt*Ae
    nue=2*mp.re(rho*rho*m*moment)
    nueta=-2*mp.im(rho*rho*m*moment)
    R=abs(Z)**2/(2*nu)
    F=(mp.re(Ze*mp.conj(Z))/nu-abs(Z)**2*nue/(2*nu*nu))/(2*rho)
    Preal=(-mp.im(Ze*mp.conj(Z))/nu-abs(Z)**2*nueta/(2*nu*nu))/(2*rho)
    out=dict(F=F,Preal=Preal,resource=R)
    if not direct:return out
    # Solve both traces, using the conjugate-profile source and divided
    # difference particular field. This does not differentiate a resource.
    tb=mp.conj(t);dt=tb-t;sigma=1 if pol=='TE' else -1
    if dt:
        Ub,_,Wb,_=boundary(tb,pol);Ud=(Ub-U)/dt;Wd=(Wb-W)/dt
    else:Ud,Wd=Ut,Wt
    c,h=mp.lu_solve(mp.matrix([[1j*rho*m*Ud,U],[Wd,W/(1j*rho*m)]]),mp.matrix([Ae,-1j*Be]))
    work=mp.mpc(0);enorm=mp.mpf(0);hnorm=mp.mpf(0);jnorm=mp.mpf(0)
    for r,weight,u,ut,cu,cut in samples:
        v,_,cv,_=profiles(tb,r,pol)
        dd=[(a-b)/dt for a,b in zip(v,u)] if dt else ut
        cd=[(a-b)/dt for a,b in zip(cv,cu)] if dt else cut
        E=[1j*rho*m*c*x+h*y for x,y in zip(dd,u)]
        H=[c*x+h*y/(1j*rho*m) for x,y in zip(cd,cu)]
        J=[c*x for x in v]
        work-=weight*sum(x*mp.conj(y) for x,y in zip(E,J))/2
        enorm+=weight*sum(abs(x)**2 for x in E);hnorm+=weight*sum(abs(x)**2 for x in H)
        jnorm+=weight*sum(abs(x)**2 for x in J)/2
    D=Ud*W-U*Wd
    out.update(direct_F=mp.im(work),direct_Preal=mp.re(work),direct_resource=jnorm,
               boundary_determinant_error=abs(D-sigma*nu)/nu,
               c_formula_error=abs(c-1j*Z/(sigma*nu))/max(1,abs(c)),
               power_balance_error=abs(mp.re(work)-1/(2*rho*rho)-rho*(mp.im(e)*enorm+mp.im(m)*hnorm)/2),
               reactive_error=abs(F-mp.im(work)),real_error=abs(Preal-mp.re(work)))
    return out

def mixed(x,y,w,etae=0,etam=0,rho=3,n=40):
    e=mp.mpc(x,etae);m=mp.mpc(y,etam);rho=mp.mpf(rho)
    tm=channel(e,m,rho,'TM',n);te=channel(e,m,rho,'TE',n)
    return tm['F']+w*te['F']

def continuation(case,etae,etam,n=40):
    x0,y0,w0=map(mp.mpf,[case['epsilon'],case['mu'],case['w']]);rho=mp.mpf(case['rho'])
    fn=lambda x,y,w:mixed(x,y,w,etae,etam,rho,n)
    # Internal real t and permeability give a much better conditioned chart
    # than raw epsilon/mu near the very anisotropic definite example.
    fr=lambda tr,y,w:fn((tr/(rho*rho)+etae*etam)/y,y,w)
    sysr=lambda tr,y,w:(fr(tr,y,w),mp.diff(lambda z:fr(z,y,w),tr),mp.diff(lambda z:fr(tr,z,w),y))
    tr0=rho*rho*(x0*y0-etae*etam)
    sol=mp.findroot(sysr,(tr0,y0,w0),tol=mp.mpf('1e-37'),maxsteps=50)
    tr,y,w=sol;x=(tr/(rho*rho)+etae*etam)/y
    H=mp.matrix([[mp.diff(lambda z:fn(z,y,w),x,2),mp.diff(lambda z:mp.diff(lambda q:fn(z,q,w),y),x)],
                 [mp.diff(lambda z:mp.diff(lambda q:fn(z,q,w),y),x),mp.diff(lambda z:fn(x,z,w),y,2)]])
    residual=max(abs(fn(x,y,w)),abs(mp.diff(lambda z:fn(z,y,w),x)),abs(mp.diff(lambda z:fn(x,z,w),y)))
    ev=mp.eigsy(H,eigvals_only=True)
    return dict(type=case['classification'],epsilon_real=x,mu_real=y,w=w,eta_e=mp.mpf(etae),eta_m=mp.mpf(etam),
                epsilon_shift=x-x0,mu_shift=y-y0,w_shift=w-w0,residual=residual,
                hessian_eigenvalues=list(ev),hessian_determinant=mp.det(H))

def serial(v):
    if isinstance(v,dict):return {k:serial(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [serial(x) for x in v]
    if isinstance(v,str):return v
    return mp.nstr(v,42)

def main():
    cases=json.loads((ROOT/'data/high_precision_critical_points.json').read_text());audits=[];rows=[];derivs=[]
    for case in cases:
        x,y,w,r=map(mp.mpf,[case['epsilon'],case['mu'],case['w'],case['rho']])
        # Local implicit continuation: controls (epsilon_R,mu_R,w) may retune.
        for control in ['electric','magnetic']:
            ge=lambda x,y,eta:mixed(x,y,w,eta if control=='electric' else 0,eta if control=='magnetic' else 0,r)
            Feta=mp.diff(lambda z:ge(x,y,z),mp.mpf(0))
            g=[mp.diff(lambda xx:mp.diff(lambda z:ge(xx,y,z),mp.mpf(0)),x),
               mp.diff(lambda yy:mp.diff(lambda z:ge(x,yy,z),mp.mpf(0)),y)]
            H=mp.matrix([[mp.mpf(v) for v in rr] for rr in case['hessian']]);v=mp.matrix([mp.mpf(0),mp.mpf(0)])
            dw=-Feta/mp.mpf(case['F_w'])
            fwgrad=[mp.diff(lambda xx:mixed(xx,y,mp.mpf(1),0,0,r)-mixed(xx,y,mp.mpf(0),0,0,r),x),
                    mp.diff(lambda yy:mixed(x,yy,mp.mpf(1),0,0,r)-mixed(x,yy,mp.mpf(0),0,0,r),y)]
            v=mp.lu_solve(H,-mp.matrix([g[i]+fwgrad[i]*dw for i in range(2)]))
            dr=dict(type=case['classification'],loss=control,F_loss_derivative=Feta,epsilon_rate=v[0],mu_rate=v[1],w_rate=dw)
            derivs.append(dr);print(json.dumps(serial(dr)),flush=True)
        for eta in [mp.mpf('0.00001'),mp.mpf('0.0001'),mp.mpf('0.001')]:
            for control in ['electric','magnetic']:
                row=continuation(case,eta if control=='electric' else 0,eta if control=='magnetic' else 0)
                row['loss']=control;rows.append(row);print(json.dumps(serial(row)),flush=True)
        for pol in ['TM','TE']:
            out=channel(mp.mpc(x,'0.0001'),mp.mpc(y,'0.0001'),r,pol,40,True)
            out.update(type=case['classification'],polarization=pol);audits.append(out)
            assert out['reactive_error']<mp.mpf('1e-25') and out['real_error']<mp.mpf('1e-25')
            assert out['power_balance_error']<mp.mpf('1e-25') and out['c_formula_error']<mp.mpf('1e-25')
    (ROOT/'audit/weak_loss.json').write_text(json.dumps(serial(dict(derivatives=derivs,continuations=rows,direct_checks=audits,
       status='Analytic local persistence theorem; finite-loss coordinates are numerical continuations, not new interval certificates.')),indent=2))
    with (ROOT/'data/fig10_weak_loss.csv').open('w',newline='') as f:
        wr=csv.writer(f);wr.writerow(['type','loss','eta_e','eta_m','epsilon_real','mu_real','w','epsilon_shift','mu_shift','w_shift','residual','hessian_determinant'])
        for x in rows:wr.writerow([x['type'],x['loss']]+[mp.nstr(x[k],36) for k in ['eta_e','eta_m','epsilon_real','mu_real','w','epsilon_shift','mu_shift','w_shift','residual','hessian_determinant']])
    print('PASS: lossy Cauchy fields, source norm, resource derivatives, and power balance.',flush=True)

if __name__=='__main__':main()
