"""Run the independent reference workflow without overwriting archived data.

The optional reference implementation has its own sampling grids and plots.
It is copied into an isolated output tree before it runs. For the publication
layout, use the Mathematica CSV-only notebook with the archived data instead.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]


def snapshot():
    paths=[ROOT/'MaxwellManifoldFeasibility.m',ROOT/'RunAllFigures.nb',
           ROOT/'Maxwell_Manifold_clean_figures.nb',ROOT/'audit/critical_point_certificates.json']
    paths+=list((ROOT/'data').glob('*'))+list((ROOT/'clean_figures').glob('*.pdf'))
    return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if p.is_file()}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'independent_run',help='Isolated reference output directory')
    parser.add_argument('--checks-only',action='store_true',help='Run exact certificates and symbolic checks; skip the figure studies')
    args=parser.parse_args();out=args.output.resolve()
    protected=[ROOT,ROOT/'analysis',ROOT/'data',ROOT/'audit',ROOT/'clean_figures']
    if out in ROOT.parents or any(out==p or p in out.parents for p in protected[1:]) or out==ROOT:
        parser.error('Choose a separate output tree, not a canonical input directory.')
    before=snapshot()
    (out/'data').mkdir(parents=True,exist_ok=True);(out/'audit').mkdir(exist_ok=True)
    shutil.copytree(ROOT/'analysis',out/'analysis',dirs_exist_ok=True,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    for name in ['data/high_precision_critical_points.json','audit/critical_point_certificates.json']:
        shutil.copy2(ROOT/name,out/name)
    jobs=['verify_certificates.py','small_size_certificate.py','verify_loss_direction.py',
          'static_symbolic.py','geometry_structure_audit.py']
    if not args.checks_only:
        jobs+=['build_study.py','higher_mode_persistence.py','spheroid_check.py','critical_map.py',
               'extended_study.py','weak_loss.py','plot_weak_loss.py']
    print('Independent workflow output:',out,flush=True)
    for name in jobs:
        print('\nRunning',name,flush=True)
        subprocess.run([sys.executable,str(out/'analysis'/name)],cwd=out,check=True)
    if not args.checks_only:
        summary=json.loads((out/'audit/study_summary.json').read_text())
        a=summary['independent_Maxwell_checks']
        if not (a['max_modal_relative_error']<1e-10 and a['max_power_balance_relative_error']<1e-8):
            raise RuntimeError('Independent Maxwell check failed')
        if max(abs(x['direct_F']) for x in a['critical_point_direct_checks'])>=1e-11:
            raise RuntimeError('Critical-point direct volume-work check failed')
        s=summary['asymptotic_slopes']
        for key,low,high in [('power_leading_slope',1.9,2.1),('power_second_order_slope',3.9,4.1),('gradient_leading_slope',1.9,2.1),('gradient_second_order_slope',3.9,4.1)]:
            if not low<s[key]<high:raise RuntimeError('Asymptotic slope check failed: '+key)
        for value in summary['critical_scaling_slopes'].values():
            if not (.48<value['distance_slope']<.52 and -.52<value['kappa_slope']<-.48):
                raise RuntimeError('Critical scaling check failed')
        sph=json.loads((out/'audit/spheroid_check.json').read_text());last=sph['convergence'][-1]
        if not (last['trace_E']<2e-8 and last['trace_H']<3e-9 and last['affine_current_error']<2e-9 and last['quadratic_power_error']<1e-11):
            raise RuntimeError('Spheroidal convergence check failed')
        if not (len(sph['root_checks'])==1 and abs(sph['root_checks'][0]['q'])<1e-11 and sph['w1_discriminant']<0):
            raise RuntimeError('Spheroidal reference root check failed')
        extra=json.loads((out/'audit/asymptotic_extended_cases.json').read_text())
        for row in extra:
            if not all(low<row[key]<high for key,low,high in [('power_leading_slope',1.95,2.05),('power_corrected_slope',3.95,4.05),('gradient_leading_slope',1.95,2.05),('gradient_corrected_slope',3.95,4.05)]):
                raise RuntimeError('Additional asymptotic check failed')
        loss=json.loads((out/'audit/weak_loss_independent_checks.json').read_text())
        if max(abs(float(x['direct_tuning_residual'])) for x in loss['checks'])>=1e-25 or loss['max_direct_power_balance_error']>=1e-25:
            raise RuntimeError('Weak-loss direct Maxwell check failed')
    if snapshot()!=before:
        raise RuntimeError('A canonical archived file changed unexpectedly.')
    print('\nPASS: independent workflow completed; canonical Mathematica data and publication figures are unchanged.',flush=True)


if __name__=='__main__':main()
