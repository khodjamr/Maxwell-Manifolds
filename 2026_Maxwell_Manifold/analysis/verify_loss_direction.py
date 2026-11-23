"""Exact bounds on the full-material gradient and magnetic control shifts.

The two supplied rational critical-point boxes are used; there is no new root
search. This program and its dependencies use only the Python standard library.
"""
from pathlib import Path
from fractions import Fraction as Q
import json
from verify_certificates import interval_namespace
ROOT=Path(__file__).resolve().parents[1]


def verify_cases():
    ns=interval_namespace();I=ns['I'];out=[]
    for c in json.loads((ROOT/'audit/critical_point_certificates.json').read_text()):
        cen=[Q(z) for z in c['centre_rational']];rad=[Q(z) for z in c['box_radii_rational']]
        k,w=[I(x-r,x+r) for x,r in zip(cen,rad)];rho=c['rho']
        tm,te=ns['coeff_jets'](k,rho)
        a=tm[0].c[0]+w*te[0].c[0];b=tm[1].c[0]+w*te[1].c[0];mu=-b/(2*a)
        fl=(1+w+2*rho**5*(w*te[1].c[0]-tm[1].c[0]))/(2*rho*rho*mu)
        fw=te[0].c[0]*mu+te[1].c[0]+te[2].c[0]/mu
        if fl.lo<=0<=fl.hi or fw.lo<=0<=fw.hi:
            raise RuntimeError('Required loss/control derivative sign is not certified.')
        rate=-fl/fw
        out.append(dict(type=c['type'],F_mu_imag_enclosure=fl.bounds(),
                        F_w_enclosure=fw.bounds(),magnetic_control_rate_enclosure=rate.bounds(),
                        full_four_real_dimensional_gradient_nonzero=True,
                        meaning='Criticality is sectional; the full complex-material tuning set is regular at this lossless point.'))
    return out


def main():
    out=verify_cases()
    (ROOT/'audit/loss_direction_certificate.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2));print('PASS: full-material regularity and magnetic control-rate signs at both certified points.')


if __name__=='__main__':main()
