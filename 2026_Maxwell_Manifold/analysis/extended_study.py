"""Revision figures: spheroid, non-exhaustive transition map, and design reading.

Every plotted point is computed from the Maxwell equations or their explicit
coefficients. This script does not generate synthetic research observations.
"""
from pathlib import Path
import csv,json
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from build_study import save,phi0,phi2,mpF,e_static,CASES,write,ROOT,COLORS
from refine_critical import disc
mp.mp.dps=55

def figure1():
    out=json.loads((ROOT/'audit/spheroid_check.json').read_text());rows=out['convergence'];ls=[x['lmax'] for x in rows]
    fig,ax=plt.subplots(figsize=(7.4,4.5))
    for key,label,mark in [('trace_E','Electric trace error','o'),('trace_H','Magnetic trace error','s'),
                           ('affine_current_error','Affine-current discrepancy','^'),('quadratic_power_error','Held-out power discrepancy','d')]:
        ax.semilogy(ls,[max(x[key],1e-16) for x in rows],mark+'-',label=label)
    ax.set(xlabel='Largest interior angular order',ylabel='Relative error or discrepancy',
           title='Prolate spheroid: semiaxes $(1,1,1.3)$',xticks=[1,5,9,13,17,21])
    ax.grid(which='both',alpha=.15);ax.legend(fontsize=11);save(fig,'fig1_spheroid_convergence')

def figure8():
    out=json.loads((ROOT/'audit/critical_map.json').read_text());pp=out['candidates']
    fig,axs=plt.subplots(1,2,figsize=(11,4.3))
    cr=np.genfromtxt(ROOT/'data/fig8_discriminant_root_curves.csv',delimiter=',',names=True)
    for rho,col in zip([1.,1.5,2.,3.],['#185a8d','#e18a24','#8f3174','#24855d']):
        used=False
        for branch in [0,1]:
            ss=cr[(cr['rho']==rho)&(cr['branch']==branch)];k=ss['k'];w=ss['w']
            # Split jumps and missing portions; never connect a pole in w.
            split=np.where((np.diff(k)>.005)|(abs(np.diff(w))>20))[0]+1
            for ids in np.split(np.arange(len(k)),split):
                axs[0].semilogy(k[ids],w[ids],color=col,lw=1.1,label=rf'$\rho={rho:g}$' if not used else None);used=True
    for kind,col,mark in [('saddle','#185a8d','o'),('definite','#b44335','s')]:
        ss=[p for p in pp if p['type']==kind]
        axs[1].semilogy([p['rho'] for p in ss],[p['w'] for p in ss],mark,color=col,ms=4,label=kind.capitalize()+' candidate')
    for c in CASES:
        axs[1].plot(float(c['rho']),float(c['w']),'k*',ms=12)
    axs[1].plot([],[],'k*',ms=12,label='Interval-certified examples')
    axs[0].set(xlabel=r'Internal wavenumber $k$',ylabel=r'Positive discriminant root $w_\pm$',ylim=(1e-3,200),title='(a) Curves whose stationary points are critical')
    axs[1].set(xlabel=r'Exterior electrical size $\rho$',ylabel=r'Critical composition $w$',title='(b) Candidates inside the stated contrast window')
    for ax in axs:ax.legend(fontsize=10);ax.grid(which='both',alpha=.15)
    save(fig,'fig8_critical_map')

def resource(e,m,r,w):
    k=r*mp.sqrt(e*m);j0=mp.sin(k)/k;j=mp.sin(k)/k**2-mp.cos(k)/k
    d=j0-j/k;vm=(j*j-j0*(3*j/k-j0))/2;vn=vm+j*d/k
    he=mp.sin(r)/r**2-mp.cos(r)/r+1j*(-mp.cos(r)/r**2-mp.sin(r)/r)
    h0=(mp.sin(r)-1j*mp.cos(r))/r;de=h0-he/r
    vals=[]
    for U,W,nu,Ae,Be in [(d,k*j,vn,de,he),(j/k,d,vm/(k*k),he,de)]:
        nu=mp.re(nu);assert nu>0
        Z=U*Be-W*Ae/(r*m);vals.append(abs(Z)**2/(2*nu))
    return mp.re((vals[0]+w*vals[1])/(1+w))

def loop_points(delta_w,n=161):
    c=CASES[1];r,kc,wc=map(mp.mpf,[c['rho'],c['k'],c['w']]);w=wc+mp.mpf(str(delta_w))
    scale=mp.sqrt(-2*mp.mpf(c['D_w'])*mp.mpf(str(delta_w))/mp.mpf(c['D_kk']))
    roots=[mp.findroot(lambda k:disc(k,r,w),(kc+sgn*.95*scale,kc+sgn*1.05*scale),tol=mp.mpf('1e-48')) for sgn in [-1,1]]
    mid=sum(roots)/2;half=(roots[1]-roots[0])/2
    from refine_critical import abc
    vals=[]
    for i in range(n):
        phase=mp.pi*i/(n-1);k=mid+half*mp.cos(phase);a,b,_=abc(k,r,w);D=disc(k,r,w)
        if D<0:
            assert abs(D)<mp.mpf('1e-44');D=mp.mpf(0)
        for branch in [1,-1]:
            m=(-b+branch*mp.sqrt(D))/(2*a);e=k*k/(r*r*m)
            vals.append((float(phase if branch==1 else 2*mp.pi-phase),e,m,w))
    return sorted(vals,key=lambda z:z[0])

def figure9():
    c=CASES[1];ec,mc,r,wc=map(mp.mpf,[c['epsilon'],c['mu'],c['rho'],c['w']]);base=resource(ec,mc,r,wc);rows=[];widths=[]
    fig,axs=plt.subplots(1,3,figsize=(11.7,4.1))
    for dw,col in [(1e-7,COLORS[0]),(1e-5,COLORS[1])]:
        vals=loop_points(dw);ee=[];mm=[];ang=[];res=[];kaps=[]
        for ph,e,m,w in vals:
            R=resource(e,m,r,w)/base
            ge=mp.re(mp.diff(lambda x:mpF(x,m,r,w),e));gm=mp.re(mp.diff(lambda x:mpF(e,x,r,w),m))
            kap=(1+w)/(2*r*r*mp.sqrt((e*ge)**2+(m*gm)**2))
            ee.append(float(100*(e/ec-1)));mm.append(float(100*(m/mc-1)));ang.append(ph/(2*np.pi));res.append(float(R));kaps.append(float(kap))
            rows.append([dw,ph,e,m,100*(e/ec-1),100*(m/mc-1),R,kap])
        axs[0].plot(ee,mm,color=col,lw=2.5 if dw==1e-7 else 1.5,zorder=3 if dw==1e-7 else 2,label=rf'$w-w_*={dw:.0e}$')
        if dw==1e-5:
            axs[1].plot(ang,res,color=col);axs[2].semilogy(ang,kaps,color=col)
        widths.append(dict(delta_w=dw,max_relative_epsilon_percent=max(abs(x) for x in ee),
                           max_relative_mu_percent=max(abs(x) for x in mm),resource_min=min(res),resource_max=max(res),
                           relative_sensitivity_min=min(kaps),relative_sensitivity_max=max(kaps)))
    more=[]
    for dw in [1e-9,1e-8,1e-6,1e-4]:
        vs=loop_points(dw,n=81)
        more.append(dict(delta_w=dw,max_relative_epsilon_percent=max(float(abs(100*(x[1]/ec-1))) for x in vs),
                         max_relative_mu_percent=max(float(abs(100*(x[2]/mc-1))) for x in vs)))
    axs[0].set(xlabel=r'$100(\varepsilon/\varepsilon_*-1)$',ylabel=r'$100(\mu/\mu_*-1)$',title='(a) Relative material excursions')
    axs[0].legend(fontsize=12);axs[1].set(xlabel='Position around the loop',ylabel=r'$\mathcal{E}_0/\mathcal{E}_{0,*}$ at fixed radiated power',title='(b) Source resource along the larger loop')
    axs[2].set(xlabel='Position around the loop',ylabel=r'$\kappa_{\rm rel}$',title='(c) Relative retuning sensitivity')
    for ax in axs:ax.grid(alpha=.15)
    write('fig9_design_interpretation.csv',['delta_w','phase_radians','epsilon','mu','relative_epsilon_percent','relative_mu_percent','fixed_power_resource_ratio','relative_metric_inverse_sensitivity'],[[float(z) for z in row] for row in rows])
    save(fig,'fig9_design_interpretation')
    out=dict(reference='same total radiated power at every w; reference is the critical point',
             main_loops=widths,additional_widths=more)
    (ROOT/'audit/design_interpretation.json').write_text(json.dumps(out,indent=2));return out

def broadened_asymptotics():
    cases=[('original',-3.,-2.5,1.69),('near_TE_root',-2.4,-8.9,.4),
           ('coordinate_fold',e_static(-np.sqrt(18),.7),-np.sqrt(18),.7),
           ('epsilon_negative',-2.1,3.,.3),('mu_negative',3.,-2.1,.8)]
    rows=[];summary=[]
    for name,ef,mf,wf in cases:
        e,m,w=map(lambda x:mp.mpf(str(x)),[ef,mf,wf]);local=[]
        for rf in np.geomspace(.008,.10,11):
            r=mp.mpf(str(rf));p=lambda x,y:mp.re(r**5*mpF(x,y,r,w))
            p0=lambda x,y:phi0(x,y,w);p2=lambda x,y:p0(x,y)+r*r*phi2(x,y,w)
            g=[mp.diff(lambda x:p(x,m),e),mp.diff(lambda x:p(e,x),m)]
            g0=[mp.diff(lambda x:p0(x,m),e),mp.diff(lambda x:p0(e,x),m)]
            g2=[mp.diff(lambda x:p2(x,m),e),mp.diff(lambda x:p2(e,x),m)]
            er=[abs(p(e,m)-p0(e,m)),abs(p(e,m)-p2(e,m)),mp.sqrt(sum((x-y)**2 for x,y in zip(g,g0))),mp.sqrt(sum((x-y)**2 for x,y in zip(g,g2)))]
            row=[name,rf,ef,mf,wf]+[float(z) for z in er];rows.append(row);local.append(row)
        a=np.array([x[1:2]+x[5:] for x in local],float)
        slopes=[float(np.polyfit(np.log(a[:6,0]),np.log(a[:6,j]),1)[0]) for j in range(1,5)]
        summary.append(dict(case=name,epsilon=ef,mu=mf,w=wf,power_leading_slope=slopes[0],power_corrected_slope=slopes[1],gradient_leading_slope=slopes[2],gradient_corrected_slope=slopes[3]))
    write('fig4_extended_asymptotics.csv',['case','rho','epsilon','mu','w','power_leading_error','power_corrected_error','gradient_leading_error','gradient_corrected_error'],rows)
    (ROOT/'audit/asymptotic_extended_cases.json').write_text(json.dumps(summary,indent=2));return summary

def modified_bessel_audit():
    rows=[]
    for l in range(1,5):
        for kk in [mp.mpf('.4'),mp.mpf('1.2')]:
            L=l*(l+1)
            il=lambda z:mp.sqrt(mp.pi/(2*z))*mp.besseli(l+mp.mpf('.5'),z)
            dil=lambda z:mp.diff(il,z)+il(z)/z
            for pol in ['TE','TM']:
                if pol=='TE':
                    U=lambda z:il(z)/z**l;W=lambda z:z**(1-l)*dil(z)
                    nu=kk**(-2*l)*mp.quad(lambda rr:rr*rr*il(kk*rr)**2,[0,1]);sigma=1
                else:
                    U=lambda z:z**(1-l)*dil(z);W=lambda z:-z**(2-l)*il(z)
                    nu=kk**(-2*l)*mp.quad(lambda rr:L*il(kk*rr)**2+kk*kk*rr*rr*dil(kk*rr)**2,[0,1]);sigma=-1
                det=(-mp.diff(U,kk)*W(kk)+U(kk)*mp.diff(W,kk))/(2*kk)
                err=abs(det-sigma*nu)/nu
                rows.append(dict(ell=l,kappa=float(kk),polarization=pol,nu=str(nu),determinant_relative_error=str(err)))
    assert max(mp.mpf(x['determinant_relative_error']) for x in rows)<mp.mpf('1e-45')
    (ROOT/'audit/modified_bessel_audit.json').write_text(json.dumps(rows,indent=2))

if __name__=='__main__':
    figure1();figure8();print(json.dumps(figure9(),indent=2),flush=True)
    print(json.dumps(broadened_asymptotics(),indent=2),flush=True);modified_bessel_audit()
