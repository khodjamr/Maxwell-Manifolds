"""Figure 10 and additional independent weak-loss checks."""
from pathlib import Path
import json,csv
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from build_study import save
from weak_loss import channel,mixed,ROOT
mp.mp.dps=55

def main():
    a=json.loads((ROOT/'audit/weak_loss.json').read_text())
    fig,axs=plt.subplots(1,2,figsize=(10.5,4.2));checks=[]
    for kind,col in [('saddle','#185a8d'),('definite','#e18a24')]:
        for loss,style in [('electric','o-'),('magnetic','s--')]:
            rows=sorted([x for x in a['continuations'] if x['type']==kind and x['loss']==loss],key=lambda x:float(x['eta_e'])+float(x['eta_m']))
            eta=np.array([float(x['eta_e'])+float(x['eta_m']) for x in rows]);dw=np.array([float(x['w_shift']) for x in rows])
            axs[0].loglog(eta,abs(dw),style,color=col,label=kind.capitalize()+', '+loss+' loss')
            if loss=='magnetic':
                rate=float(next(x['w_rate'] for x in a['derivatives'] if x['type']==kind and x['loss']=='magnetic'))
                axs[1].semilogx(eta,dw/(rate*eta),'o-',color=col,label=kind.capitalize())
            row=rows[1];e=mp.mpc(row['epsilon_real'],row['eta_e']);m=mp.mpc(row['mu_real'],row['eta_m']);w=mp.mpf(row['w'])
            tm=channel(e,m,mp.mpf(3),'TM',56,True);te=channel(e,m,mp.mpf(3),'TE',56,True)
            direct=tm['direct_F']+w*te['direct_F'];resource=tm['F']+w*te['F']
            base=mixed(mp.re(e),mp.re(m),w,mp.im(e),mp.im(m),3,40)
            assert abs(direct)<mp.mpf('1e-25') and abs(resource-base)<mp.mpf('1e-25')
            checks.append(dict(type=kind,loss=loss,direct_tuning_residual=mp.nstr(direct,35),
                               quadrature_40_to_56_change=mp.nstr(abs(resource-base),35),
                               max_power_balance_error=mp.nstr(max(tm['power_balance_error'],te['power_balance_error']),35)))
    axs[0].set(xlabel=r'Loss component $\eta$',ylabel=r'Critical-control shift $|w_c-w_*|$',title='(a) Critical-composition shifts')
    axs[1].axhline(1,color='k',ls='--',lw=1)
    axs[1].set(xlabel=r'Magnetic loss $\eta_\mu$',ylabel=r'$(w_c-w_*)/(\eta_\mu\,\partial_{\eta_\mu}w_c|_0)$',title='(b) Magnetic-loss derivative check')
    handles,labels=axs[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,-.12),ncol=2,fontsize=11)
    for ax in axs:ax.grid(which='both',alpha=.15)
    axs[1].legend(fontsize=11)
    save(fig,'fig10_weak_loss')
    out=dict(checks=checks,max_direct_power_identity_error=max(float(x['reactive_error']) for x in a['direct_checks']),
             max_direct_real_power_identity_error=max(float(x['real_error']) for x in a['direct_checks']),
             max_direct_power_balance_error=max(float(x['power_balance_error']) for x in a['direct_checks']))
    (ROOT/'audit/weak_loss_independent_checks.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2));print('PASS: finite-loss critical-point volume work and quadrature refinement.')

if __name__=='__main__':main()
