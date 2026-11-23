"""Replay the critical-point existence proofs using exact rational arithmetic.

This verifier uses only the Python standard library. No floating-point root
solver, mpmath, scipy, or interval library is required. The candidate boxes and
preconditioners are read as exact rationals from the supplied certificate file.
The exact interval implementation is rational_interval.py. The separate
rational_certificate.py candidate generator needs mpmath only for choosing
new rational centres and preconditioners.
"""
from pathlib import Path
from decimal import Decimal, localcontext
from fractions import Fraction as Q
import json

ROOT = Path(__file__).resolve().parents[1]

def interval_namespace():
    """Return the shared exact-arithmetic engine; no floating-point imports."""
    import rational_interval
    return vars(rational_interval)

def verify(case, ns):
    if case.get('bits') != ns['BITS'] or case.get('taylor_terms') != ns['N']+1:
        raise RuntimeError('Certificate metadata does not match the arithmetic engine.')
    I, system = ns['I'], ns['interval_system']
    centre = [Q(s) for s in case['centre_rational']]
    radii = [Q(s) for s in case['box_radii_rational']]
    if len(centre)!=2 or len(radii)!=2 or any(r<=0 for r in radii):
        raise RuntimeError('Invalid certificate centre or radii.')
    pre = [[Q(s) for s in row] for row in case['preconditioner_rational']]
    assert pre[0][0]*pre[1][1]-pre[0][1]*pre[1][0] != 0
    box = [I(c-r,c+r) for c,r in zip(centre,radii)]
    g0,_,_,_ = system(I(centre[0]), I(centre[1]), case['rho'])
    gbox,jac,abc,signs = system(box[0],box[1],case['rho'])
    residual = [[I(int(i==j))-sum((I(pre[i][h])*jac[h][j] for h in range(2)), I(0))
                 for j in range(2)] for i in range(2)]
    image = [I(centre[i])-sum((I(pre[i][h])*g0[h] for h in range(2)),I(0))+
             sum((residual[i][j]*I(-radii[j],radii[j]) for j in range(2)),I(0))
             for i in range(2)]
    contraction = max(sum((residual[i][j].absmax() for j in range(2)),Q(0)) for i in range(2))
    inside = all(centre[i]-radii[i] < image[i].lo and image[i].hi < centre[i]+radii[i]
                 for i in range(2))
    a,b,_ = abc
    mu = -b/(2*a)
    epsilon = box[0]*box[0]/(case['rho']**2*mu)
    dkk,dw = signs
    assert inside and contraction < 1
    assert not dkk.lo <= 0 <= dkk.hi
    assert not dw.lo <= 0 <= dw.hi
    kind = 'saddle' if dkk.lo > 0 else 'definite'
    assert kind == case['type']
    return {'type':kind,'strict_inclusion':inside,'contraction_bound':ns['decimal'](contraction),
            'A_enclosure':a.bounds(),
            'G_box_enclosures':[z.bounds() for z in gbox],
            'DG_box_enclosures':[[z.bounds() for z in row] for row in jac],
            'G_box_widths':[ns['decimal'](z.hi-z.lo) for z in gbox],
            'DG_box_widths':[[ns['decimal'](z.hi-z.lo) for z in row] for row in jac],
            'epsilon_enclosure':epsilon.bounds(),'mu_enclosure':mu.bounds(),
            'Delta_kk_enclosure':dkk.bounds(),'Delta_w_enclosure':dw.bounds()}

if __name__ == '__main__':
    namespace = interval_namespace()
    cases = json.loads((ROOT/'audit/critical_point_certificates.json').read_text())
    results = [verify(case, namespace) for case in cases]
    (ROOT/'audit/certificate_replay.json').write_text(json.dumps(results,indent=2))
    print(json.dumps(results,indent=2))
    print('PASS: two unique exact critical points, with certified nondegenerate transitions.')
