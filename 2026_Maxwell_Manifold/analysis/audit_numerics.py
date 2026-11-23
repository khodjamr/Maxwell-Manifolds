"""Check the supplied CSV values against independent Maxwell evaluations.

This optional audit needs numpy, scipy and mpmath. It reads the figure data
without changing them, and writes only verification_results/numerical_consistency.json.
The certificate verifiers are separate and need only the standard library.
"""
from pathlib import Path
from collections import Counter
import csv
import json
import math
import numpy as np
import mpmath as mp
from modal_coefficients import coefficients, mixed_coeff, power, discriminant
from refine_critical import abc, coeff, radial
from weak_loss import channel as lossy_channel, mixed as lossy_mixed

ROOT = Path(__file__).resolve().parents[1]
mp.mp.dps = 55
REPORT = {"scope": "CSV consistency and independent numerical evaluation, not interval certification", "checks": {}}


def rows(name):
    with (ROOT / 'data' / (name + '.csv')).open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))


def number(value):
    return mp.mpf(value)


def mf(e, m, rho, w):
    a, b, c = abc(rho * mp.sqrt(e * m), rho, w)
    return mp.re(a * m + b + c / m)


def mg(e, m, rho, w):
    return (mp.diff(lambda x: mf(x, m, rho, w), e),
            mp.diff(lambda x: mf(e, x, rho, w), m))


def p0(e, m, w):
    return (e + 2) / 4 - w * (m + 2) * (m + 9) / (14 * m)


def p2(e, m, w):
    return -m * (e + 1) / 10 - e / 4 + w * (-2 * e * m*m + 105*m + 34*e*m + 223*e - 1890/m) / 1470


def resource(e, m, rho, w):
    k = rho * mp.sqrt(e*m)
    j, d, _, _, vm, vn, _, _ = radial(k)
    h = mp.sin(rho)/rho**2 - mp.cos(rho)/rho + 1j*(-mp.cos(rho)/rho**2-mp.sin(rho)/rho)
    h0 = (mp.sin(rho)-1j*mp.cos(rho))/rho
    dh = h0-h/rho
    ans = 0
    for U, W, nu, Ae, Be, weight in [(d, k*j, vn, dh, h, 1), (j, k*d, vm, h, dh, w)]:
        ans += weight * abs(U*Be-W*Ae/(rho*m))**2/(2*nu)
    return mp.re(ans)


def check(name, metrics, limits=None):
    for key, limit in (limits or {}).items():
        if not math.isfinite(float(metrics[key])) or not float(metrics[key]) <= limit:
            raise RuntimeError('%s: %s=%s exceeds %s' % (name, key, metrics[key], limit))
    REPORT['checks'][name] = metrics
    print('PASS:', name, flush=True)


def vector_zero_errors(table, rho):
    ee = np.array([float(x['epsilon']) for x in table]); mm = np.array([float(x['mu']) for x in table]); ww = np.array([float(x['w']) for x in table])
    kk = rho * np.sqrt(ee*mm)
    a, b, c = mixed_coeff(kk, rho, ww)
    ff = a*mm+b+c/mm
    scale = np.maximum(1, np.abs(a*mm)+np.abs(b)+np.abs(c/mm))
    return {'max_scaled_tuning_residual_from_rounded_coordinates': float(np.max(np.abs(ff)/scale)),
            'max_internal_k_relative_error': float(np.max(np.abs(kk-np.array([float(x['internal_k']) for x in table]))/kk))}


def main():
    cases = json.loads((ROOT/'data/high_precision_critical_points.json').read_text())
    cases = {x['classification']: x for x in cases}
    for name in ['fig1_spheroid_convergence','fig2_static_profiles','fig3_regular_tuning_graph','fig4_asymptotic_convergence','fig5_saddle_reconnection','fig6_closed_loop_birth','fig7_critical_scaling','fig8_critical_map','fig8_discriminant_root_curves','fig9_design_interpretation','fig10_weak_loss']:
        rr = rows(name)
        if not rr:
            raise RuntimeError('Empty CSV: '+name)
        for row in rr:
            if None in row:
                raise RuntimeError('Incorrect column count: '+name)
            for v in row.values():
                if v and v.lower() in {'nan','inf','-inf','infinity','indeterminate'}:
                    raise RuntimeError('Nonfinite cell in '+name)

    rr=rows('fig1_spheroid_convergence')
    errors=[abs(float(r['Delta'])-(float(r['B'])**2-4*float(r['A'])*float(r['C']))) for r in rr]
    check('figure_1', {'rows':len(rr),'discriminant_identity_error':max(errors),'last_trace_E':float(rr[-1]['trace_E']),'last_trace_H':float(rr[-1]['trace_H']),'last_current_discrepancy':float(rr[-1]['affine_current_error']),'last_power_discrepancy':float(rr[-1]['quadratic_power_error']),'max_condition':max(float(r['max_condition']) for r in rr)}, {'discriminant_identity_error':1e-14,'last_trace_E':2e-8,'last_trace_H':3e-9,'last_current_discrepancy':2e-9,'last_power_discrepancy':1e-11})

    rr=rows('fig2_static_profiles');errors=[];static_errors=[]
    for pol in ['TM','TE']:
        tt=[r for r in rr if r['polarization']==pol]
        e,m,rh=np.array([[float(r[k]) for r in tt] for k in ['epsilon','mu','rho']])
        kk=rh*np.sqrt(e*m);a,b,c=coefficients(1,pol,kk,rh)
        actual=rh**5*(a*m+b+c/m)
        stored=np.array([float(r['scaled_exact_q']) for r in tt])
        expected=(e+2)/4 if pol=='TM' else -(m+2)*(m+9)/(14*m)
        errors.extend(abs(actual-stored));static_errors.extend(abs(expected-np.array([float(r['scaled_static_q']) for r in tt])))
    check('figure_2',{'rows':len(rr),'max_scaled_power_difference':float(max(errors)),'max_static_difference':float(max(static_errors))},{'max_scaled_power_difference':2e-11,'max_static_difference':1e-14})

    rr=rows('fig3_regular_tuning_graph');metrics={'rows':len(rr)};ee=[];gg=[];ss=[];cc=[]
    for r in rr:
        rho,w,m,e=map(number,[r[k] for k in ['rho','w','mu','epsilon_exact']])
        e0=-2+(2*w/7)*(m+11+18/m)
        ee.append(abs(rho**5*mf(e,m,rho,w)));gg.append(abs(rho**5*mg(e,m,rho,w)[0]-number(r['scaled_F_epsilon_exact'])))
        ss.append(abs(e0-number(r['epsilon_static'])));cc.append(abs(e0-4*rho*rho*p2(e0,m,w)-number(r['epsilon_second_order'])))
    metrics.update(max_scaled_tuning_residual=float(max(ee)),max_gradient_difference=float(max(gg)),max_static_graph_difference=float(max(ss)),max_corrected_graph_difference=float(max(cc)),min_sampled_scaled_derivative=min(float(r['scaled_F_epsilon_exact']) for r in rr),max_sampled_scaled_derivative=max(float(r['scaled_F_epsilon_exact']) for r in rr))
    check('figure_3',metrics,{'max_scaled_tuning_residual':1e-14,'max_gradient_difference':1e-14,'max_static_graph_difference':1e-14,'max_corrected_graph_difference':1e-14})

    rr=rows('fig4_asymptotic_convergence');diff=[]
    e,m,w=mp.mpf(-3),mp.mpf('-2.5'),mp.mpf('1.69')
    for r in rr:
        rho=number(r['rho']);exact=rho**5*mf(e,m,rho,w);lead=p0(e,m,w);cor=lead+rho*rho*p2(e,m,w)
        grad=[rho**5*x for x in mg(e,m,rho,w)]
        g0=[mp.diff(lambda x:p0(x,m,w),e),mp.diff(lambda x:p0(e,x,w),m)]
        g2=[g0[0]+rho*rho*mp.diff(lambda x:p2(x,m,w),e),g0[1]+rho*rho*mp.diff(lambda x:p2(e,x,w),m)]
        vals=[abs(exact-lead),abs(exact-cor),mp.sqrt(sum((a-b)**2 for a,b in zip(grad,g0))),mp.sqrt(sum((a-b)**2 for a,b in zip(grad,g2)))]
        for key,value in zip(['power_leading_error','power_second_order_error','gradient_leading_error','gradient_second_order_error'],vals):
            diff.append(abs(value-number(r[key]))/max(abs(value),mp.mpf('1e-40')))
    x=np.log([float(r['rho']) for r in rr[:7]])
    slopes={key:float(np.polyfit(x,np.log([float(r[key]) for r in rr[:7]]),1)[0]) for key in ['power_leading_error','power_second_order_error','gradient_leading_error','gradient_second_order_error']}
    check('figure_4',{'rows':len(rr),'max_relative_error_in_stored_errors':float(max(diff)),'fit_points':7,'fit_slopes':slopes},{'max_relative_error_in_stored_errors':1e-10})

    rr=rows('fig5_saddle_reconnection');metrics=vector_zero_errors(rr,3);metrics['rows']=len(rr)
    check('figure_5',metrics,{'max_scaled_tuning_residual_from_rounded_coordinates':2e-11,'max_internal_k_relative_error':1e-14})

    rr=rows('fig6_closed_loop_birth');loop=[r for r in rr if r['dataset']=='loop'];dd=[r for r in rr if r['dataset']=='discriminant'];metrics=vector_zero_errors(loop,3);metrics.update(loop_rows=len(loop),discriminant_rows=len(dd))
    points=np.linspace(0,len(dd)-1,min(81,len(dd))).astype(int);err=[]
    for ix in points:
        r=dd[ix];k,w=map(number,[r['internal_k'],r['w']]);a,b,c=abc(k,mp.mpf(3),w)
        err.append(abs(b*b-4*a*c-number(r['Delta'])))
    c=cases['definite'];kc,wc,dkk,dw=map(number,[c[k] for k in ['k','w','D_kk','D_w']]);delta=mp.mpf('1e-7');sx=mp.sqrt(-2*dw*delta/dkk);sy=mp.sqrt(dw*delta)
    errcoords=[]
    for ix in np.linspace(0,len(loop)-1,81).astype(int):
        r=loop[ix];k,w,mu=map(number,[r['internal_k'],r['w'],r['mu']]);a,b,_=abc(k,mp.mpf(3),w)
        errcoords.append(max(abs((k-kc)/sx-number(r['normal_form_x'])),abs((2*a*mu+b)/sy-number(r['normal_form_y']))))
    metrics.update(discriminant_samples=len(points),max_discriminant_difference=float(max(err)),normal_coordinate_samples=81,max_normal_coordinate_difference=float(max(errcoords)))
    check('figure_6',metrics,{'max_scaled_tuning_residual_from_rounded_coordinates':2e-11,'max_internal_k_relative_error':1e-14,'max_discriminant_difference':1e-15,'max_normal_coordinate_difference':1e-8})

    rr=rows('fig7_critical_scaling');metrics={'rows':len(rr),'fits':{},'smallest_increment_ratios':{}};errors={k:[] for k in ['distance','sensitivity','distance_prediction','sensitivity_prediction']}
    for typ in ['saddle','definite']:
        c=cases[typ];e0,m0,w0=map(number,[c[k] for k in ['epsilon','mu','w']]);H=mp.matrix([[number(v) for v in row] for row in c['hessian']]);h=min(mp.eigsy(H,eigvals_only=True));fw=number(c['F_w']);tt=[r for r in rr if r['type']==typ]
        for r in tt:
            e,m,delta=map(number,[r[k] for k in ['epsilon','mu','delta_w']]);distance=mp.sqrt((e-e0)**2+(m-m0)**2);g=mg(e,m,mp.mpf(3),w0+delta);kap=1/mp.sqrt(sum(v*v for v in g));pred=mp.sqrt(-2*fw*delta/h);kpred=1/mp.sqrt(2*abs(fw*h*delta))
            for key,value,col in [('distance',distance,'material_distance_exact'),('sensitivity',kap,'inverse_gradient_exact'),('distance_prediction',pred,'distance_quadratic_prediction'),('sensitivity_prediction',kpred,'inverse_gradient_quadratic_prediction')]:
                errors[key].append(abs(value-number(r[col]))/abs(value))
        x=np.log([float(r['delta_w']) for r in tt[:6]])
        metrics['fits'][typ]={key:float(np.polyfit(x,np.log([float(r[key]) for r in tt[:6]]),1)[0]) for key in ['material_distance_exact','inverse_gradient_exact']}
        first=min(tt,key=lambda r:float(r['delta_w']))
        metrics['smallest_increment_ratios'][typ]={'distance':float(first['material_distance_exact'])/float(first['distance_quadratic_prediction']),'sensitivity':float(first['inverse_gradient_exact'])/float(first['inverse_gradient_quadratic_prediction'])}
    metrics.update({k+'_max_relative_difference':float(max(v)) for k,v in errors.items()})
    check('figure_7',metrics,{'distance_max_relative_difference':2e-8,'sensitivity_max_relative_difference':2e-8,'distance_prediction_max_relative_difference':1e-13,'sensitivity_prediction_max_relative_difference':1e-13})

    rr=rows('fig8_critical_map');counts=Counter(float(r['rho']) for r in rr);stored=rows('fig8_candidate_counts')
    if any(int(r['total'])!=counts.get(float(r['rho']),0) for r in stored):raise RuntimeError('Figure 8 count table mismatch')
    errors=[];derivs=[];signs_ok=True
    for r in rr:
        k,rho,w=map(number,[r[k] for k in ['k','rho','w']]);a,b,c=abc(k,rho,w);dd=lambda x: (lambda q:q[1]**2-4*q[0]*q[2])(abc(x,rho,w));scale=max(mp.mpf(1),abs(b*b)+abs(4*a*c));errors.append(abs(dd(k))/scale);derivs.append(abs(mp.diff(dd,k))/scale)
        e,mu=number(r['epsilon']),number(r['mu'])
        if not (mp.mpf('.1')<=abs(e)<=100 and mp.mpf('.1')<=abs(mu)<=100 and e*mu>0):raise RuntimeError('Figure 8 contrast window violated')
        signs_ok &= (r['type']=='saddle') == (number(r['Delta_kk'])>0)
    if not signs_ok:raise RuntimeError('Figure 8 classification sign mismatch')
    curves=rows('fig8_discriminant_root_curves');sample=np.linspace(0,len(curves)-1,121).astype(int);ce=[]
    for ix in sample:
        r=curves[ix];k,rho,w=map(number,[r[k] for k in ['k','rho','w']]);a,b,c=abc(k,rho,w);ce.append(abs(b*b-4*a*c)/max(mp.mpf(1),abs(b*b)+abs(4*a*c)))
    check('figure_8',{'candidates':len(rr),'saddle':sum(r['type']=='saddle' for r in rr),'definite':sum(r['type']=='definite' for r in rr),'counts':dict(sorted(counts.items())),'first_detected_rho':min(counts),'scan_grid_min':min(float(r['rho']) for r in stored),'curve_rows':len(curves),'curve_samples':len(sample),'max_scaled_Delta':float(max(errors)),'max_scaled_Delta_k':float(max(derivs)),'max_sampled_root_curve_scaled_Delta':float(max(ce))},{'max_scaled_Delta':1e-12,'max_scaled_Delta_k':1e-11,'max_sampled_root_curve_scaled_Delta':1e-11})

    rr=rows('fig9_design_interpretation');c=cases['definite'];e0,m0,w0=map(number,[c[k] for k in ['epsilon','mu','w']]);rho=mp.mpf(3);base=resource(e0,m0,rho,w0)/(1+w0);reserr=[];kaperr=[];coorderr=[];metrics={'rows':len(rr),'loops':{}}
    for delta in sorted(set(r['delta_w'] for r in rr),key=float):
        tt=[r for r in rr if r['delta_w']==delta];metrics['loops'][delta]={'max_abs_epsilon_excursion_percent':max(abs(float(r['relative_epsilon_percent'])) for r in tt),'max_abs_mu_excursion_percent':max(abs(float(r['relative_mu_percent'])) for r in tt),'resource_ratio_min':min(float(r['fixed_power_resource_ratio']) for r in tt),'resource_ratio_max':max(float(r['fixed_power_resource_ratio']) for r in tt)}
    for r in rr:
        e,m,delta=map(number,[r[k] for k in ['epsilon','mu','delta_w']]);w=w0+delta
        coorderr.extend([abs(100*(e/e0-1)-number(r['relative_epsilon_percent'])),abs(100*(m/m0-1)-number(r['relative_mu_percent']))])
    for ix in np.linspace(0,len(rr)-1,82).astype(int):
        r=rr[ix];e,m,delta=map(number,[r[k] for k in ['epsilon','mu','delta_w']]);w=w0+delta;res=resource(e,m,rho,w)/((1+w)*base);g=mg(e,m,rho,w);kap=(1+w)/(2*rho*rho*mp.sqrt((e*g[0])**2+(m*g[1])**2))
        reserr.append(abs(res-number(r['fixed_power_resource_ratio'])));kaperr.append(abs(kap-number(r['relative_metric_inverse_sensitivity']))/kap)
    metrics.update(resource_samples=len(reserr),max_resource_ratio_difference=float(max(reserr)),max_sensitivity_relative_difference=float(max(kaperr)),max_excursion_difference_percent=float(max(coorderr)))
    check('figure_9',metrics,{'max_resource_ratio_difference':1e-13,'max_sensitivity_relative_difference':1e-8,'max_excursion_difference_percent':1e-12})

    rr=rows('fig10_weak_loss');metrics={'rows':len(rr),'cases':[]};errs=[];balances=[];norms=[];grads=[]
    for r in rr:
        e,m,w,etae,etam=map(number,[r[k] for k in ['epsilon_real','mu_real','w','eta_e','eta_m']]);rho=mp.mpf(3)
        tm=lossy_channel(mp.mpc(e,etae),mp.mpc(m,etam),rho,'TM',40,True);te=lossy_channel(mp.mpc(e,etae),mp.mpc(m,etam),rho,'TE',40,True)
        direct=tm['direct_F']+w*te['direct_F'];normal=tm['F']+w*te['F'];ge=mp.diff(lambda x:lossy_mixed(x,m,w,etae,etam,rho),e);gm=mp.diff(lambda x:lossy_mixed(e,x,w,etae,etam,rho),m)
        errs.append(abs(direct));grads.append(max(abs(ge),abs(gm)));balances.extend([tm['power_balance_error'],te['power_balance_error']]);norms.extend([tm['direct_resource'],te['direct_resource']])
        metrics['cases'].append({'type':r['type'],'loss':r['loss'],'eta':float(etae+etam),'w_shift':float(r['w_shift']),'direct_tuning_residual_from_rounded_coordinates':float(abs(direct)),'reactive_identity_error':float(abs(normal-direct)),'max_real_material_derivative':float(max(abs(ge),abs(gm)))})
    metrics.update(max_direct_tuning_residual_from_rounded_coordinates=float(max(errs)),max_real_material_derivative_from_rounded_coordinates=float(max(grads)),max_direct_power_balance_error=float(max(balances)),min_direct_source_resource=float(min(norms)))
    check('figure_10',metrics,{'max_direct_tuning_residual_from_rounded_coordinates':2e-12,'max_real_material_derivative_from_rounded_coordinates':2e-10,'max_direct_power_balance_error':1e-30})

    REPORT['status']='PASS'
    (ROOT/'verification_results').mkdir(exist_ok=True)
    (ROOT/'verification_results/numerical_consistency.json').write_text(json.dumps(REPORT,indent=2))
    print('PASS: supplied figure data are consistent with the checked Maxwell evaluations.',flush=True)


if __name__=='__main__':
    main()
