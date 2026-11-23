"""Verify the dipole series and scaled-power coefficients in Appendix B.

All coefficients are exact rationals. This finite algebra check supports the
printed formulas; the analytic remainder theorem remains in the manuscript.
Requires sympy (only the exact certificate replayers are standard-library-only).
"""
from pathlib import Path
import json
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
t,z,e,m,w=s.symbols('t zeta epsilon mu w',nonzero=True,real=True)
a=[s.Rational((-1)**n*(2*n+2),s.factorial(2*n+3)) for n in range(7)]
f=sum(c*t**n for n,c in enumerate(a))
v=sum(sum(a[j]*a[n-j] for j in range(n+1))*t**n/(2*n+5) for n in range(7))
W=2*(f+t*s.diff(f,t));ntm=t*v+f*W
series=lambda expr,n=4:s.series(expr,t,0,n).removeO().expand()
te=[series(f*f/v),series(f*W/v),series(W*W/v)]
tm=[series(W*W/ntm),series(t*f*W/ntm),series(t*t*f*f/ntm)]
expected_te=[5-2*t/7-2*t*t/735,10-11*t/7+17*t*t/735,20-36*t/7+223*t*t/735]
expected_tm=[2-2*t/5+t*t/350,t-t*t/10,t*t/2]
for actual,expected in zip(te+tm,expected_te+expected_tm):
    if s.expand(series(actual,3)-expected)!=0:
        raise RuntimeError('A printed dipole trace-ratio coefficient is inconsistent.')
if tm[2].coeff(t,3)!=0:
    raise RuntimeError('The TM R22 cubic coefficient must vanish.')
T=series(f*f/ntm);S=2*T+t*s.diff(T,t)
qte=m*(1-z+z*z)*s.diff(te[0],t)+2*s.diff(te[1],t)+(1+z)*s.diff(te[2],t)/m
qtm=m*z*(1+z)*s.diff(tm[0],t)+2*s.diff(tm[1],t)+e*(1-z+z*z)*S
mixed=s.series(((qtm+w*qte)/4).subs(t,z*e*m),z,0,2).removeO().expand()
phi0=(e+2)/4-w*(m+2)*(m+9)/(14*m)
phi2=-m*(e+1)/10-e/4+w*(-2*e*m*m+105*m+34*e*m+223*e-1890/m)/1470
if s.simplify(mixed.coeff(z,0)-phi0)!=0 or s.simplify(mixed.coeff(z,1)-phi2)!=0:
    raise RuntimeError('A scaled-power coefficient is inconsistent.')
out={'status':'PASS','trace_ratios_TE':[str(x) for x in expected_te],
     'trace_ratios_TM':[str(x) for x in expected_tm],
     'TM_R22_cubic_coefficient':0,'Phi0':str(s.factor(phi0)),'Phi2':str(s.factor(phi2)),
     'scope':'Exact finite series algebra; uniform analytic remainder proved in manuscript.'}
(ROOT/'audit').mkdir(exist_ok=True)
(ROOT/'audit/symbolic_series_check.json').write_text(json.dumps(out,indent=2))
print('PASS: six dipole trace-ratio series, vanishing TM cubic coefficient, Phi0 and Phi2.')
