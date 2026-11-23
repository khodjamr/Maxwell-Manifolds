"""Independent Cauchy-boundary computation for a prolate spheroid.

This is a converged regular-wave collocation experiment, not a certified
continuum discretization. Electric and magnetic traces are fitted separately
at each material value. No affine-current law is imposed in the solver.
"""
from pathlib import Path
import csv, json
import numpy as np
from scipy.special import spherical_jn, spherical_yn, eval_legendre, factorial2
from numpy.polynomial.legendre import leggauss

ROOT=Path(__file__).resolve().parents[1]

def geometry(x, aspect):
    r=((1-x*x)+x*x/aspect**2)**(-.5)
    q=x*np.sqrt(1-x*x)*(1/aspect**2-1)*r*r
    nr=1/np.sqrt(1+q*q); nt=-q*nr
    return r,nr,nt

def wave(l,k,r,x,pol,outgoing=False,normalise=True):
    """Return u, du/dt, curl u, d(curl u)/dt in spherical components."""
    r,x=np.broadcast_arrays(r,x);z=k*r;L=l*(l+1)
    y=np.sqrt((2*l+1)/(4*np.pi))*eval_legendre(l,x)
    # Nodes avoid the poles, where this quotient has a removable limit.
    pl=eval_legendre(l,x);pm=eval_legendre(l-1,x)
    p=np.sqrt((2*l+1)/(4*np.pi))*l*(x*pl-pm)/(np.sqrt(L)*np.sqrt(1-x*x))
    j=spherical_jn(l,z);jp=spherical_jn(l,z,True)
    if outgoing:
        j=j+1j*spherical_yn(l,z);jp=jp+1j*spherical_yn(l,z,True)
    jpp=-2*jp/z-(1-L/z**2)*j
    d=jp+j/z;dp=jpp+jp/z-j/z**2
    zero=np.zeros_like(j)
    m=np.array([zero,zero,-p*j]);mk=np.array([zero,zero,-p*r*jp])
    n=np.array([np.sqrt(L)*y*j/z,p*d,zero])
    nk=np.array([np.sqrt(L)*y*r*(jp/z-j/z**2),p*r*dp,zero])
    u,uk,v,vk=(m,mk,n,nk) if pol=='TE' else (n,nk,m,mk)
    s=(l if pol=='TE' else l-1) if normalise else 0
    factor=factorial2(2*l+1)/k**s if normalise else 1
    return factor*u,factor*(uk-s*u/k)/(2*k),factor*k*v,factor*((1-s)*v+k*vk)/(2*k)

def tangential(v,nr,nt,pol):
    return v[2] if pol=='TE' else nr*v[1]-nt*v[0]

def solve_channel(rho,k,mu,aspect,pol,lmax,ntrain=100,ncheck=173,nrad=60,nangle=110):
    orders=list(range(1,lmax+1,2));x,wx=leggauss(ntrain)
    r,nr,nt=geometry(x,aspect);surf=2*np.pi*wx*r*r/nr
    out=wave(1,rho,r,x,pol,True,False)
    et=tangential(out[0],nr,nt,pol)
    other='TM' if pol=='TE' else 'TE'
    ht=tangential(out[2]/(1j*rho),nr,nt,other)
    scales=[np.sqrt(np.sum(surf*abs(v)**2)) for v in [et,ht]]
    cols=[]
    for l in orders:
        u,ut,cu,cut=wave(l,k,r,x,pol)
        cols.append((1j*rho*mu*tangential(ut,nr,nt,pol),tangential(cut,nr,nt,other)))
    for l in orders:
        u,ut,cu,cut=wave(l,k,r,x,pol)
        cols.append((tangential(u,nr,nt,pol),tangential(cu/(1j*rho*mu),nr,nt,other)))
    mat=np.array([np.concatenate([np.sqrt(surf)*c[i]/scales[i] for i in [0,1]]) for c in cols]).T
    rhs=np.concatenate([np.sqrt(surf)*v/sc for v,sc in zip([et,ht],scales)])
    cs=np.linalg.norm(mat,axis=0)
    sol,_,rank,sv=np.linalg.lstsq(mat/cs,rhs,rcond=1e-14);sol=sol/cs
    c,h=sol[:len(orders)],sol[len(orders):]
    # Held-out trace residuals use different nodes and more sampling points.
    xx,ww=leggauss(ncheck);rr,nrr,ntt=geometry(xx,aspect);sw=2*np.pi*ww*rr*rr/nrr
    out=wave(1,rho,rr,xx,pol,True,False)
    eout=tangential(out[0],nrr,ntt,pol);hout=tangential(out[2]/(1j*rho),nrr,ntt,other)
    ee=np.zeros_like(eout);hh=np.zeros_like(hout)
    for l,cc,hhh in zip(orders,c,h):
        u,ut,cu,cut=wave(l,k,rr,xx,pol)
        ee+=tangential(1j*rho*mu*cc*ut+hhh*u,nrr,ntt,pol)
        hh+=tangential(cc*cut+hhh*cu/(1j*rho*mu),nrr,ntt,other)
    er=np.sqrt(np.sum(sw*abs(ee-eout)**2)/np.sum(sw*abs(eout)**2))
    hr=np.sqrt(np.sum(sw*abs(hh-hout)**2)/np.sum(sw*abs(hout)**2))
    xa,wa=leggauss(nangle);sr,wr=leggauss(nrad);s=(sr+1)/2;wr=wr/2
    ra,_,_=geometry(xa,aspect);rv=ra[:,None]*s[None,:];xv=xa[:,None]+np.zeros_like(rv)
    vw=2*np.pi*wa[:,None]*wr[None,:]*ra[:,None]**3*s[None,:]**2
    j=np.zeros((3,nangle,nrad),complex);e=j.copy()
    for l,cc,hhh in zip(orders,c,h):
        u,ut,_,_=wave(l,k,rv,xv,pol)
        j+=cc*u;e+=1j*rho*mu*cc*ut+hhh*u
    inner=lambda a,b:np.sum(vw*np.sum(a*np.conj(b),axis=0))
    power=-inner(e,j)/2;res=inner(j,j).real/2
    return dict(c=c,h=h,J=j,E=e,weights=vw,q=power.imag,resource=res,
                balance_error=abs(power.real*2*rho*rho-1),e_trace=float(er),h_trace=float(hr),
                condition=float(sv[0]/sv[-1]),rank=int(rank))

def main():
    rho=.8;k=1.2;aspect=1.3;w=0.1;mus=[.7,1.0,1.6,2.3]
    rows=[];final={}
    for lmax in [1,3,5,7,9,11,13,15,17,19,21]:
        solutions={pol:[solve_channel(rho,k,mu,aspect,pol,lmax) for mu in mus] for pol in ['TM','TE']}
        fit_errors=[]
        qmix=sum(weight*np.array([s['q'] for s in solutions[pol]]) for pol,weight in [('TM',1),('TE',w)])
        # Three independently solved powers determine a test polynomial; fourth is held out.
        X=np.array([[mu,1,1/mu] for mu in mus]);abc=np.linalg.solve(X[:3],qmix[:3])
        qerror=abs(X[3]@abc-qmix[3])/max(1,abs(qmix[3]))
        for pol in solutions:
            ss=solutions[pol];a,b=ss[0],ss[2];ma,mb=mus[0],mus[2]
            jc=(mb*b['J']-ma*a['J'])/(mb-ma);jr=ma*mb*(a['J']-b['J'])/(mb-ma)
            for mu,ss0 in zip(mus,ss):
                diff=jc+jr/mu-ss0['J']
                fit_errors.append(np.sqrt(np.sum(ss0['weights']*np.sum(abs(diff)**2,axis=0))/(2*ss0['resource'])))
        allsol=[s for ss in solutions.values() for s in ss]
        row=dict(lmax=lmax,rho=rho,k=k,aspect=aspect,w=w,
                 trace_E=max(s['e_trace'] for s in allsol),trace_H=max(s['h_trace'] for s in allsol),
                 real_power_error=max(s['balance_error'] for s in allsol),
                 affine_current_error=max(fit_errors),quadratic_power_error=qerror,
                 A=abc[0],B=abc[1],C=abc[2],Delta=abc[1]**2-4*abc[0]*abc[2],
                 max_condition=max(s['condition'] for s in allsol))
        rows.append(row);print(json.dumps(row),flush=True)
        final=solutions
    with (ROOT/'data/fig1_spheroid_convergence.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    roots=np.roots([rows[-1]['A'],rows[-1]['B'],rows[-1]['C']])
    rootchecks=[]
    for mu in roots:
        if abs(mu.imag)<1e-12 and mu.real>0:
            mu=float(mu.real);ss=[solve_channel(rho,k,mu,aspect,p,21) for p in ['TM','TE']]
            rootchecks.append(dict(mu=mu,epsilon=k*k/(rho*rho*mu),q=ss[0]['q']+w*ss[1]['q'],
                                   resource=ss[0]['resource']+w*ss[1]['resource'],
                                   trace_E=max(s['e_trace'] for s in ss),trace_H=max(s['h_trace'] for s in ss)))
    per_pol={pol:np.linalg.solve(X[:3],[s['q'] for s in final[pol]][:3]) for pol in final}
    exclusion=per_pol['TM']+per_pol['TE']
    # Independent sphere reference checks the angular signs and volume weights.
    from modal_coefficients import coefficients
    benchmarks=[]
    for pol in ['TM','TE']:
        ss=solve_channel(rho,k,1.6,1.0,pol,5)
        qref=coefficients(1,pol,k,rho)@np.array([1.6,1,1/1.6])
        benchmarks.append(dict(polarization=pol,relative_q_error=abs(ss['q']-qref)/max(1,abs(qref)),
                               e_trace=ss['e_trace'],h_trace=ss['h_trace']))
    # A denser quadrature/collocation repeats a solved material point.
    refined=[]
    for pol in ['TM','TE']:
        ss=solve_channel(rho,k,1.6,aspect,pol,21,ntrain=160,ncheck=241,nrad=85,nangle=160)
        old=final[pol][2]
        refined.append(dict(polarization=pol,q_change=abs(ss['q']-old['q'])/max(1,abs(old['q'])),
                            resource_change=abs(ss['resource']-old['resource'])/old['resource'],e_trace=ss['e_trace'],h_trace=ss['h_trace']))
    out=dict(method='regular-wave Cauchy collocation; numerical evidence, not an error certificate',
             geometry='prolate spheroid with semiaxes (1,1,1.3)',mus=mus,convergence=rows,root_checks=rootchecks,
             unit_polarization_ABC={p:xs.tolist() for p,xs in per_pol.items()},
             w1_exclusion_ABC=exclusion.tolist(),w1_discriminant=float(exclusion[1]**2-4*exclusion[0]*exclusion[2]),
             sphere_benchmarks=benchmarks,quadrature_refinements=refined)
    (ROOT/'audit/spheroid_check.json').write_text(json.dumps(out,indent=2))

if __name__=='__main__':main()
