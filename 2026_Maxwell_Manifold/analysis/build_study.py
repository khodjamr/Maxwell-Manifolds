"""Generate research figures and their explicitly labelled data tables.

The plots use the exact Maxwell coefficients. Floating-point illustrations are
separate from the standard-library rational existence proof. Run from any cwd.
"""
from pathlib import Path
import csv
import json
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from modal_coefficients import coefficients, mixed_coeff, discriminant, power
from direct_maxwell import channel
from direct_high_precision import channel as channel_hp
from refine_critical import abc as mpabc

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT/'figures'; DATA=ROOT/'data'
FIG.mkdir(exist_ok=True); DATA.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':14,'axes.labelsize':14,'legend.fontsize':13,
                     'axes.titlesize':14,'lines.linewidth':1.7,'pdf.fonttype':42,
                     'savefig.bbox':'tight','figure.constrained_layout.use':True})
COLORS=['#185a8d','#e18a24','#8f3174']
mp.mp.dps=55
CASES=json.loads((DATA/'high_precision_critical_points.json').read_text())

def write(name,header,rows):
    with (DATA/name).open('w',newline='') as f:
        w=csv.writer(f);w.writerow(header);w.writerows(rows)

def save(fig,name):
    fig.savefig(FIG/(name+'.pdf'))
    fig.savefig(FIG/(name+'.png'),dpi=150)
    plt.close(fig)

def phi0(e,m,w): return (e+2)/4-w*(m+2)*(m+9)/(14*m)
def phi2(e,m,w): return -m*(e+1)/10-e/4+w*(-2*e*m*m+105*m+34*e*m+223*e-1890/m)/1470
def e_static(m,w): return -2+(2*w/7)*(m+11+18/m)
def mpF(e,m,r,w):
    k=r*mp.sqrt(e*m)
    a,b,c=mpabc(k,r,w)
    return a*m+b+c/m

# ======================================================================
# FIGURE 2 AND CSV: explicit electric and magnetic dipole static profiles
# ======================================================================
def figure2():
    fig,axs=plt.subplots(1,2,figsize=(10.5,4.0));rows=[]
    ee=np.linspace(-5.5,-1.1,361);mm=np.linspace(-12,-1.1,461)
    for r,col in zip([.05,.15,.3],COLORS):
        vals=[]
        for e in ee:
            a,b,c=coefficients(1,'TM',r*np.sqrt(e*(-3)),r)
            val=r**5*(a*(-3)+b+c/(-3));vals.append(val)
            rows.append(['TM',r,e,-3,val,(e+2)/4])
        axs[0].plot(ee,vals,color=col,label=rf'$\rho={r:g}$')
        vals=[]
        for m in mm:
            a,b,c=coefficients(1,'TE',r*np.sqrt((-3)*m),r)
            val=r**5*(a*m+b+c/m);vals.append(val)
            rows.append(['TE',r,-3,m,val,-(m+2)*(m+9)/(14*m)])
        axs[1].plot(mm,vals,color=col,label=rf'$\rho={r:g}$')
    axs[0].plot(ee,(ee+2)/4,'k--',label='Static limit')
    axs[1].plot(mm,-(mm+2)*(mm+9)/(14*mm),'k--',label='Static limit')
    for ax in axs:
        ax.axhline(0,color='.5',lw=.7);ax.grid(alpha=.15);ax.legend(loc='best')
        ax.set_ylabel(r'$\rho^5 q$')
    axs[0].set(xlabel=r'Relative permittivity $\varepsilon$',title=r'(a) Electric-type dipole, $\mu=-3$')
    axs[1].set(xlabel=r'Relative permeability $\mu$',title=r'(b) Magnetic-type dipole, $\varepsilon=-3$')
    axs[1].plot([-9,-2],[0,0],'ko',ms=4)
    write('fig2_static_profiles.csv',['polarization','rho','epsilon','mu','scaled_exact_q','scaled_static_q'],rows)
    save(fig,'fig2_static_profiles')

# ======================================================================
# FIGURE 3 AND CSV: small-size tuning graph and regular coordinate fold
# ======================================================================
def figure3():
    w=1.69;ms=np.linspace(-12,-1.4,401);es=e_static(ms,w)
    r=.16;exact=[];second=[];grad=[];rows=[]
    for m,e0 in zip(ms,es):
        e=brentq(lambda x:float(power(x,m,r,w)),e0-.8,e0+.8,xtol=1e-13)
        e2=e0-4*r*r*phi2(e0,m,w)
        h=1e-4
        fe=r**5*(float(power(e-2*h,m,r,w))-8*float(power(e-h,m,r,w))+
                     8*float(power(e+h,m,r,w))-float(power(e+2*h,m,r,w)))/(12*h)
        exact.append(e);second.append(e2);grad.append(fe)
        rows.append([r,w,m,e,e0,e2,fe])
    fig,axs=plt.subplots(1,2,figsize=(10.5,4.0))
    axs[0].plot(ms,exact,color=COLORS[0],label=r'Exact Maxwell, $\rho=0.16$')
    axs[0].plot(ms,es,'k--',label='Leading graph')
    axs[0].plot(ms,second,color=COLORS[1],ls=':',label=r'Graph through $\rho^2$')
    mf=-np.sqrt(18);ef=e_static(mf,w)
    axs[0].plot([mf],[ef],'ko',ms=4,label='Static coordinate fold')
    axs[0].set(xlabel=r'$\mu$',ylabel=r'Tuned $\varepsilon$',title='(a) Mixed-dipole tuning graph')
    axs[0].legend()
    axs[1].plot(ms,grad,color=COLORS[0],label=r'Exact $\partial_\varepsilon(\rho^5F)$')
    axs[1].axhline(.25,color='k',ls='--',label='Static derivative $1/4$')
    axs[1].axvline(mf,color='.6',ls=':',lw=1)
    axs[1].set(xlabel=r'$\mu$',ylabel='Permittivity derivative',title='(b) A regular coordinate fold')
    axs[1].legend()
    for ax in axs:ax.grid(alpha=.2)
    write('fig3_regular_tuning_graph.csv',['rho','TE_to_TM_power_ratio_w','mu','epsilon_exact','epsilon_static','epsilon_second_order','scaled_F_epsilon_exact'],rows)
    save(fig,'fig3_regular_tuning_graph')

# ======================================================================
# FIGURE 4 AND CSV: independent value and material-gradient asymptotic tests
# ======================================================================
def figure4():
    e=mp.mpf(-3);m=mp.mpf('-2.5');w=mp.mpf('1.69');rows=[]
    for r0 in np.geomspace(.015,.24,15):
        r=mp.mpf(str(r0));p=lambda x,y:r**5*mpF(x,y,r,w)
        p0=lambda x,y:phi0(x,y,w);p2=lambda x,y:p0(x,y)+r*r*phi2(x,y,w)
        g=[mp.diff(lambda x:p(x,m),e),mp.diff(lambda x:p(e,x),m)]
        g0=[mp.diff(lambda x:p0(x,m),e),mp.diff(lambda x:p0(e,x),m)]
        g2=[mp.diff(lambda x:p2(x,m),e),mp.diff(lambda x:p2(e,x),m)]
        rows.append([r0,float(abs(p(e,m)-p0(e,m))),float(abs(p(e,m)-p2(e,m))),
                     float(mp.sqrt(sum((a-b)**2 for a,b in zip(g,g0)))),
                     float(mp.sqrt(sum((a-b)**2 for a,b in zip(g,g2))))])
    arr=np.array(rows);fig,axs=plt.subplots(1,2,figsize=(10.5,4.0))
    for ax,i,title in [(axs[0],1,'(a) Scaled-power error'),(axs[1],3,'(b) Material-gradient error')]:
        ax.loglog(arr[:,0],arr[:,i],'o-',color=COLORS[0],label='Leading approximation')
        ax.loglog(arr[:,0],arr[:,i+1],'s-',color=COLORS[1],label=r'Approximation through $\rho^2$')
        x=arr[:,0];ax.loglog(x,arr[0,i]*(x/x[0])**2,'k--',lw=1,label=r'Order $\rho^2$')
        ax.loglog(x,arr[0,i+1]*(x/x[0])**4,'k:',lw=1,label=r'Order $\rho^4$')
        ax.set(xlabel=r'Exterior electrical size $\rho$',ylabel='Absolute error',title=title)
        ax.legend();ax.grid(which='both',alpha=.15)
    write('fig4_asymptotic_convergence.csv',['rho','power_leading_error','power_second_order_error','gradient_leading_error','gradient_second_order_error'],rows)
    save(fig,'fig4_asymptotic_convergence')
    return {'power_leading_slope':float(np.polyfit(np.log(arr[:7,0]),np.log(arr[:7,1]),1)[0]),
            'power_second_order_slope':float(np.polyfit(np.log(arr[:7,0]),np.log(arr[:7,2]),1)[0]),
            'gradient_leading_slope':float(np.polyfit(np.log(arr[:7,0]),np.log(arr[:7,3]),1)[0]),
            'gradient_second_order_slope':float(np.polyfit(np.log(arr[:7,0]),np.log(arr[:7,4]),1)[0])}

def branches(case,w,khalf,n=1801):
    rho=float(case['rho']);kc=float(case['k']);ks=np.linspace(kc-khalf,kc+khalf,n)
    ds=discriminant(ks,rho,w)
    folds=[]
    for i in range(len(ks)-1):
        if ds[i]*ds[i+1]<0:folds.append(brentq(lambda k:discriminant(k,rho,w),ks[i],ks[i+1],xtol=1e-14))
    ks=np.sort(np.r_[ks,folds,kc]);a,b,c=mixed_coeff(ks,rho,w);d=b*b-4*a*c
    # Only the tiny roundoff-sized discriminant at a known fold is clamped.
    for fold in folds:d[np.argmin(abs(ks-fold))]=0
    if abs(w-float(case['w']))<1e-14:d[np.argmin(abs(ks-kc))]=0
    d=np.where(d>=0,d,np.nan)
    mus=[(-b+sgn*np.sqrt(d))/(2*a) for sgn in [1,-1]]
    return [(ks*ks/(rho*rho*m),m,ks) for m in mus]

# ======================================================================
# FIGURE 5 AND CSV: certified saddle and predicted branch reconnection
# ======================================================================
def figure5():
    c=CASES[0];w0=float(c['w']);ec=float(c['epsilon']);mc=float(c['mu']);rows=[]
    fig,axs=plt.subplots(1,3,figsize=(11,3.9),sharex=True,sharey=True)
    for panel,(ax,dw,col) in enumerate(zip(axs,[-.02,0,.02],[COLORS[0],'#222222','#a43632'])):
        for branch,(es,ms,ks) in enumerate(branches(c,w0+dw,.12)):
            ax.plot(es,ms,color=col)
            for i,(e,m,k) in enumerate(zip(es,ms,ks)):
                if np.isfinite(e):rows.append([panel,w0+dw,branch,i,k,e,m])
        if dw==0:ax.plot([ec],[mc],'ko',ms=4)
        ax.set(xlim=(ec-.045,ec+.045),ylim=(mc-.15,mc+.15),xlabel=r'$\varepsilon$',title=rf'$w-w_*={dw:+.2f}$')
        ax.grid(alpha=.15)
    axs[0].set_ylabel(r'$\mu$')
    write('fig5_saddle_reconnection.csv',['panel','w','quadratic_branch','sample','internal_k','epsilon','mu'],rows)
    save(fig,'fig5_saddle_reconnection')

# ======================================================================
# FIGURE 6 AND CSV: positive-material closed-loop birth
# ======================================================================
def figure6():
    c=CASES[1];r=float(c['rho']);w0=float(c['w']);kc=float(c['k']);ec=float(c['epsilon']);mc=float(c['mu'])
    dw0=1e-7;dw=float(c['D_w']);dkk=float(c['D_kk'])
    kscale=np.sqrt(-2*dw*dw0/dkk);dscale=np.sqrt(dw*dw0)
    fig,axs=plt.subplots(1,2,figsize=(10.5,4.0));rows=[]
    ks=np.linspace(kc-.0008,kc+.0008,1001)
    for dw,col in zip([-dw0,0,dw0],['#185a8d','#222222','#a43632']):
        ds=discriminant(ks,r,w0+dw)
        axs[0].plot((ks-kc)*1e3,ds*1e11,color=col,label=rf'$w-w_*={dw:+.0e}$')
        for i,(k,d) in enumerate(zip(ks,ds)):rows.append(['discriminant',w0+dw,0,i,k,'','',d,'',''])
    axs[0].axhline(0,color='.6',lw=.7)
    axs[0].set(xlabel=r'$10^3(k-k_*)$',ylabel=r'$10^{11}\Delta$',title='(a) A critical discriminant maximum')
    axs[0].legend()
    for branch,(es,ms,ks) in enumerate(branches(c,w0+dw0,.0008,n=2401)):
        a,b,_=mixed_coeff(ks,r,w0+dw0)
        # The exact square completion supplies a smooth local coordinate.
        # A linear Hessian coordinate is misleading at this anisotropic point.
        coord=np.stack([(ks-kc)/kscale,(2*a*ms+b)/dscale])
        axs[1].plot(coord[0],coord[1],color='#a43632',label='Exact tuning loop' if branch==0 else None)
        for i,(e,m,k,x,y) in enumerate(zip(es,ms,ks,coord[0],coord[1])):
            if np.isfinite(e):rows.append(['loop',w0+dw0,branch,i,k,e,m,'',x,y])
    theta=np.linspace(0,2*np.pi,501)
    axs[1].plot(np.cos(theta),np.sin(theta),'k--',lw=1,label='Quadratic prediction')
    axs[1].plot([0],[0],'ko',ms=4,label=r'Critical point at $w=w_*$')
    axs[1].set(xlabel='Normal-form coordinate $x$',ylabel='Normal-form coordinate $y$',
               title='(b) A closed curve for $w>w_*$',aspect='equal',xlim=(-1.15,1.15),ylim=(-1.15,1.15))
    for ax in axs:ax.grid(alpha=.15)
    write('fig6_closed_loop_birth.csv',['dataset','w','branch','sample','internal_k','epsilon','mu','Delta','normal_form_x','normal_form_y'],rows)
    save(fig,'fig6_closed_loop_birth')

# ======================================================================
# FIGURE 7 AND CSV: square-root transition and inverse conditioning law
# ======================================================================
def figure7():
    rows=[];slopes={};fig,axs=plt.subplots(1,2,figsize=(10.5,4.0))
    for case,col,exponents in [(CASES[0],COLORS[0],(-8,-4)),(CASES[1],COLORS[1],(-10,-6))]:
        ec,mc,wc,r=[mp.mpf(case[s]) for s in ['epsilon','mu','w','rho']]
        H=mp.matrix([[mp.mpf(x) for x in row] for row in case['hessian']]);ev,V=mp.eigsy(H)
        # A negative-curvature ray has two nearby roots for positive delta w.
        idx=0
        eig=ev[idx];v=[V[i,idx] for i in range(2)];fw=mp.mpf(case['F_w']);vals=[]
        for d0 in np.geomspace(10.**exponents[0],10.**exponents[1],13):
            dw=mp.mpf(str(d0));guess=mp.sqrt(-2*fw*dw/eig)
            fn=lambda s:mpF(ec+s*v[0],mc+s*v[1],r,wc+dw)
            s=mp.findroot(fn,(guess*.95,guess*1.05),tol=mp.mpf('1e-43'))
            e=ec+s*v[0];m=mc+s*v[1]
            ge=mp.diff(lambda x:mpF(x,m,r,wc+dw),e);gm=mp.diff(lambda x:mpF(e,x,r,wc+dw),m)
            kappa=1/mp.sqrt(ge*ge+gm*gm);prediction=1/mp.sqrt(2*fw*abs(eig)*dw)
            row=[case['classification'],d0,float(s),float(guess),float(kappa),float(prediction),float(e),float(m),float(fn(s))]
            rows.append(row);vals.append(row)
        a=np.array([x[1:6] for x in vals],float);x=a[:,0]
        axs[0].loglog(x,a[:,1],'o-',color=col,label=case['classification'].capitalize())
        axs[0].loglog(x,a[:,2],'k--',lw=.9)
        axs[1].loglog(x,a[:,3],'o-',color=col,label=case['classification'].capitalize())
        axs[1].loglog(x,a[:,4],'k--',lw=.9)
        slopes[case['classification']]={'distance_slope':float(np.polyfit(np.log(x[:6]),np.log(a[:6,1]),1)[0]),
                                      'kappa_slope':float(np.polyfit(np.log(x[:6]),np.log(a[:6,3]),1)[0])}
    axs[0].set(xlabel=r'Pattern change $w-w_*>0$',ylabel='Material distance on the selected ray',title='(a) Square-root displacement')
    axs[1].set(xlabel=r'Pattern change $w-w_*>0$',ylabel=r'Inverse sensitivity $1/\|\nabla F\|$',title='(b) Inverse tuning sensitivity')
    for ax in axs:ax.legend();ax.grid(which='both',alpha=.15)
    write('fig7_critical_scaling.csv',['type','delta_w','material_distance_exact','distance_quadratic_prediction','inverse_gradient_exact','inverse_gradient_quadratic_prediction','epsilon','mu','F_residual'],rows)
    save(fig,'fig7_critical_scaling')
    return slopes

def direct_audit():
    rows=[]
    for e,m,r in [(2,3,.3),(-3,-1.7,.4),(-4,2,.3),(2,-3,.4)]:
        for l in range(1,5):
            for pol in ['TM','TE']:
                ref=channel(e,m,r,l,pol,n=160)
                original_balance_error=abs(ref['power_balance']-1)
                if original_balance_error>1e-10:
                    ref=channel_hp(e,m,r,l,pol,n=64)
                k=r*np.sqrt(complex(e*m));aa,bb,cc=coefficients(l,pol,k,r)
                exact=float(np.real(aa*m+bb+cc/m))
                rows.append({'epsilon':e,'mu':m,'rho':r,'ell':l,'polarization':pol,
                             'direct_q':ref['q'],'coefficient_q':exact,
                             'relative_error':abs(exact-ref['q'])/max(1,abs(ref['q'])),
                             'initial_double_power_balance_relative_error':original_balance_error,
                             'high_precision_recheck':bool(original_balance_error>1e-10),
                             'radiated_power_relative_error':abs(ref['power_balance']-1)})
    critical=[]
    for c in CASES:
        e,m,r,w=[float(c[s]) for s in ['epsilon','mu','rho','w']]
        te=channel(e,m,r,1,'TE',n=220);tm=channel(e,m,r,1,'TM',n=220)
        critical.append({'type':c['classification'],'direct_F':tm['q']+w*te['q'],
                         'relative_radiated_power_error':max(abs(te['power_balance']-1),abs(tm['power_balance']-1))})
    out={'max_modal_relative_error':max(row['relative_error'] for row in rows),
         'max_power_balance_relative_error':max(row['radiated_power_relative_error'] for row in rows),
         'modal_checks':rows,'critical_point_direct_checks':critical}
    (ROOT/'audit/independent_maxwell_audit.json').write_text(json.dumps(out,indent=2))
    return {k:v for k,v in out.items() if k!='modal_checks'}

if __name__=='__main__':
    figure2();figure3();conv=figure4();figure5();figure6();scales=figure7();audit=direct_audit()
    out={'asymptotic_slopes':conv,'critical_scaling_slopes':scales,'independent_Maxwell_checks':audit}
    (ROOT/'audit/study_summary.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))
