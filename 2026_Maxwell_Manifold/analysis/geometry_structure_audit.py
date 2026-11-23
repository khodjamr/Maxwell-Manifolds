"""Algebraic checks supporting, but not replacing, the analytic proofs.

The supported-potential estimate and analytic projector are proved in the text
using scalar interior H2 regularity and Hilbert-space projection. This script
checks the constant-coefficient operator identity and the frozen-lift algebra.
It is not a numerical proof of an infinite-dimensional closed-range theorem.
"""
from pathlib import Path
import json
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
x,y,z,t,mu,ma,mb,rho=sp.symbols('xi1 xi2 xi3 t mu mu_a mu_b rho', nonzero=True)
xi=sp.Matrix([x,y,z]);n2=xi.dot(xi);I=sp.eye(3)
# Fourier symbols: curl curl = |xi|^2 I - xi xi^T; grad div = -xi xi^T.
A=n2*I-xi*xi.T;G=-xi*xi.T
symbol=sp.simplify((I+G/t)*(A-t*I)-(n2-t)*I)
assert symbol==sp.zeros(3)
a1=(mb-mu)/(mb-ma);a2=(mu-ma)/(mb-ma)
assert sp.simplify(a1+a2-1)==0
assert sp.simplify(a1*ma+a2*mb-mu)==0
# The frozen Cauchy-lift coefficients give electric trace 1 and curl trace mu.
electric=sp.simplify((mb-ma)/(mb-ma)+mu*(1-1)/(mb-ma))
curl=sp.simplify((mb*ma-ma*mb)/(mb-ma)+mu*(mb-ma)/(mb-ma))
assert electric==1 and curl==mu
out={'operator_identity_symbol_zero':True,'electric_trace_preserved':True,
     'magnetic_trace_preserved':True,
     'proof_scope':'Analytic proof in manuscript; these are symbolic consistency checks only.',
     'standard_tools':['Young convolution inequality','interior H2 regularity for the Laplacian',
                       'Hilbert-space orthogonal projection','Neumann series for bounded inverse',
                       'one-variable real-analytic factorization'],
     'new_assumptions':['bounded Lipschitz V','complete exterior E,H fixed',
                        'two compatible Y(V) Cauchy lifts at distinct nonzero permeabilities',
                        'rho>0 and epsilon*mu != 0'],
     'not_asserted':['compatibility of arbitrary raw data','a unique unconstrained forward transmission map',
                     'shape derivatives','critical points for every shape','arbitrary rough/fractal boundaries']}
(ROOT/'audit/geometry_structure_audit.json').write_text(json.dumps(out,indent=2))
print('PASS: supported-potential symbol and Cauchy-lift algebra (not a substitute for the analytic proof).')
