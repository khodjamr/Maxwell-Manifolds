"""Read-only release, data and exact-certificate checks (standard library only).

Run from any working directory: python analysis/verify_release.py
No numerical Maxwell engine or third-party Python package is imported.
"""
from pathlib import Path
from collections import Counter
from fractions import Fraction as Q
import csv
import hashlib
import json
import math
import sys

sys.dont_write_bytecode=True

ROOT=Path(__file__).resolve().parents[1]


def require(test,message):
    if not test:raise RuntimeError(message)


def table(stem):
    with (ROOT/'data'/(stem+'.csv')).open(encoding='utf-8-sig',newline='') as f:
        result=list(csv.DictReader(f))
    require(bool(result),'Empty table: '+stem)
    require(all(None not in row for row in result),'Wrong number of CSV columns: '+stem)
    return result


def slope(rows,xkey,ykey,n):
    x=[math.log(float(r[xkey])) for r in rows[:n]]
    y=[math.log(float(r[ykey])) for r in rows[:n]]
    xm=sum(x)/len(x);ym=sum(y)/len(y)
    return sum((a-xm)*(b-ym) for a,b in zip(x,y))/sum((a-xm)**2 for a in x)


def main():
    manifest=ROOT/'MANIFEST_SHA256.csv'
    require(manifest.is_file(),'The release checksum manifest is missing.')
    with manifest.open(newline='') as f:entries=list(csv.DictReader(f))
    for r in entries:
        p=(ROOT/r['path']).resolve()
        require(ROOT in p.parents,'Invalid manifest path: '+r['path'])
        require(p.is_file(),'Missing file: '+r['path'])
        require(p.stat().st_size==int(r['bytes']),'File size differs: '+r['path'])
        require(hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],'Checksum differs: '+r['path'])
    print('PASS: %d archived file checksums.'%len(entries))

    expected={'fig1_spheroid_convergence':11,'fig2_static_profiles':1266,'fig3_regular_tuning_graph':241,
              'fig4_asymptotic_convergence':15,'fig5_saddle_reconnection':9518,'fig6_closed_loop_birth':5519,
              'fig7_critical_scaling':26,'fig8_critical_map':201,'fig8_discriminant_root_curves':17407,
              'fig9_design_interpretation':644,'fig10_weak_loss':12}
    for name,n in expected.items():require(len(table(name))==n,'Unexpected primary CSV row count: '+name)
    print('PASS: eleven primary CSV tables for ten figures.')
    require(len(table('fig4_extended_asymptotics'))==55,'Additional asymptotic table is incomplete.')
    require(len(table('higher_multipole_persistence'))==8,'Quadrupole persistence table is incomplete.')
    points=table('fig8_critical_map');counts=table('fig8_candidate_counts')
    require(len(counts)==21,'The critical-map scan grid must have 21 sizes.')
    for row in counts:
        rho=float(row['rho']);same=[r for r in points if float(r['rho'])==rho]
        require(int(row['total'])==len(same),'Candidate count mismatch at rho='+str(rho))
        for typ in ['saddle','definite']:
            require(int(row[typ])==sum(r['type']==typ for r in same),'Classification count mismatch')
    require(all(float(r['epsilon'])*float(r['mu'])>0 for r in points),'The map must be restricted to equal-sign material quadrants.')
    require(min(float(r['rho']) for r in points)==1.125,'First recorded candidate size differs.')
    print('PASS: auxiliary tables and non-exhaustive map counts (not an exclusion theorem).')
    regular=table('fig3_regular_tuning_graph')
    require(min(float(r['scaled_F_epsilon_exact']) for r in regular)>0,'A sampled regular-graph derivative is not positive.')
    rr=table('fig4_asymptotic_convergence')
    fits={k:slope(rr,'rho',k,7) for k in list(rr[0])[1:]}
    require(abs(fits['power_leading_error']-1.997)<.0005,'Power leading-error fit differs from manuscript.')
    require(abs(fits['power_second_order_error']-4)<.0005,'Power corrected-error fit differs from manuscript.')
    require(abs(fits['gradient_leading_error']-2)<.0005 and abs(fits['gradient_second_order_error']-4)<.0005,'Gradient fits differ from manuscript.')
    rr=table('fig7_critical_scaling')
    for typ,distance,inverse in [('saddle',.49995,-.50001),('definite',.50000,-.50000)]:
        selected=sorted([r for r in rr if r['type']==typ],key=lambda r:float(r['delta_w']))
        require(abs(slope(selected,'delta_w','material_distance_exact',6)-distance)<.000005,'Displacement fit differs from manuscript.')
        require(abs(slope(selected,'delta_w','inverse_gradient_exact',6)-inverse)<.000005,'Inverse-sensitivity fit differs from manuscript.')
    print('PASS: sampled positivity and reported asymptotic/critical fits.')
    rr=table('fig9_design_interpretation');large=[r for r in rr if math.isclose(float(r['delta_w']),1e-5,rel_tol=1e-12)]
    require(len(large)==322,'The larger design loop is incomplete.')
    lo=min(float(r['fixed_power_resource_ratio']) for r in large);hi=max(float(r['fixed_power_resource_ratio']) for r in large)
    require(abs(100*(lo-1)+.422)<.0005 and abs(100*(hi-1)-.468)<.0005,'Fixed-power resource range differs from manuscript.')
    print('PASS: fixed-power resource range and both design-loop datasets.')

    from verify_certificates import interval_namespace,verify
    from small_size_certificate import bound
    from verify_loss_direction import verify_cases
    certificates=json.loads((ROOT/'audit/critical_point_certificates.json').read_text())
    require(len(certificates)==2,'Exactly two critical-point certificates are required.')
    ns=interval_namespace();results=[verify(c,ns) for c in certificates]
    require({x['type'] for x in results}=={'saddle','definite'},'Both critical-transition types must be certified.')
    require(float(results[0]['contraction_bound'])<3.05e-14 and float(results[1]['contraction_bound'])<1.05e-14,'A stated contraction bound failed.')
    print('PASS: two exact rational critical-point contraction certificates.')
    for rho,M,lower in [(Q(1,100),100,Q(113,1000)),(Q(1,10),10,Q(9,1000))]:
        result=bound(rho,M,4,Q(1,4))
        require(Q(result['four_mixed_derivative_lower_bound'])/4>lower,'A quantitative joint-domain exclusion bound failed.')
    require(Q(bound(Q(1,200))['four_mixed_derivative_lower_bound'])/4>Q(198,1000),'The original bounded-contrast exclusion bound failed.')
    print('PASS: all three stated sufficient dipole exclusion bounds.')
    loss=verify_cases();rates={'saddle':Q('16.30703156275055'),'definite':Q('-0.02851444687408144')}
    for result in loss:
        a,b=map(Q,result['magnetic_control_rate_enclosure'])
        require(a<=rates[result['type']]<=b,'The plotted zero-loss normalization falls outside its certified derivative enclosure.')
    print('PASS: full-material loss-gradient signs and plotted magnetic derivative normalizations.')
    print('PASS: release verification completed without altering the archived files.')


if __name__=='__main__':main()
