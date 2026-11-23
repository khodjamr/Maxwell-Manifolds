"""Exact rational bound for dipole exclusion on a specified material window.

The result is conservative and is not the largest admissible electrical size.
Only Python's standard library is used, including in the interval arithmetic.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial,comb
import json
from verify_certificates import interval_namespace
ROOT=Path(__file__).resolve().parents[1]

def bound(rho=Q(1,200),M=100,wmax=4,tmax=None):
    """Bound with independent external size, internal parameter and permeability.

    If tmax is omitted, |epsilon|, |mu| <= M implies |t| <= rho^2 M^2.
    Giving tmax instead certifies the larger joint domain |t| <= tmax,
    |mu| <= M.  No separate bound on epsilon is then required.
    """
    ns=interval_namespace();I,J=ns['I'],ns['J'];s=rho*rho;r=s*M*M if tmax is None else Q(tmax)
    assert 0<r<=1
    nmax=30;aa=[Q((-1)**n*(2*n+2),factorial(2*n+3)) for n in range(nmax+1)]
    bb=[sum((aa[i]*aa[n-i] for i in range(n+1)),Q(0))/(2*n+5) for n in range(nmax+1)]
    tb=I(-r,r)
    def jet(coeff,majorant):
        vals=[]
        for m in range(4):
            val=I(0)
            for n in range(nmax,m-1,-1):val=val*tb+I(comb(n,m)*coeff[n])
            # |a_n| <= (1/3) 10^-n, |b_n| <= (1/18) 10^-n.
            # Successive tail majorants have ratio <= (1/10)*r*(32/29)<1/2.
            rem=2*majorant*comb(nmax+1,m)*Q(1,10)**(nmax+1)*r**(nmax+1-m)
            vals.append(val+I(-rem,rem))
        return J.make(vals)
    f=jet(aa,Q(1,3));v=jet(bb,Q(1,18));t=J.make([tb,1,0,0]);W=2*(f+t*f.deriv());vn=t*v+f*W
    rte=[f*f/v,f*W/v,W*W/v];rtm=[W*W/vn,t*f*W/vn];T=f*f/vn
    mt=[(2*z.c[2]).absmax() for z in rte];me=[(2*z.c[2]).absmax() for z in rtm]
    ms=(3*T.c[1]+tb*2*T.c[2]).absmax()
    tm_error=M*M*s*s*(1+s)*me[0]+2*s*M*me[1]+s*(1+s)+(1+s+s*s)*2*r*ms
    te_bound=s*(M*M*(1+s+s*s)*mt[0]+2*M*mt[1]+(1+s)*mt[2])
    margin=1-tm_error-wmax*te_bound
    return dict(rho_rational=str(rho),M=M,wmax=wmax,t_abs_bound=str(r),
                domain='contrast_box' if tmax is None else 'independent_internal_parameter',
                nu_TE_enclosure=v.c[0].bounds(),nu_TM_enclosure=vn.c[0].bounds(),
                ratio_second_derivative_bounds_TE=[ns['decimal'](x) for x in mt],
                ratio_second_derivative_bounds_TM=[ns['decimal'](x) for x in me],
                S_prime_abs_bound=ns['decimal'](ms),
                four_TM_derivative_error_bound=ns['decimal'](tm_error),
                four_TE_derivative_abs_bound=ns['decimal'](te_bound),
                four_mixed_derivative_lower_bound=ns['decimal'](margin),
                scaled_permittivity_derivative_lower_bound=ns['decimal'](margin/4),
                certified=margin>0)

if __name__=='__main__':
    out=[]
    for rr in [Q(1,100),Q(1,200),Q(1,500),Q(1,1000)]:
        try:res=bound(rr)
        except AssertionError:continue
        out.append(res);print(json.dumps(res,indent=2))
    for rr,mm in [(Q(1,100),100),(Q(1,10),10)]:
        res=bound(rr,mm,4,Q(1,4))
        out.append(res);print(json.dumps(res,indent=2))
    required=[(Q(1,100),100,Q(113,1000)),(Q(1,10),10,Q(9,1000))]
    for rr,mm,threshold in required:
        result=bound(rr,mm,4,Q(1,4))
        if Q(result['four_mixed_derivative_lower_bound'])/4<=threshold:
            raise RuntimeError('A stated joint size-contrast lower bound failed.')
    if Q(bound(Q(1,200))['four_mixed_derivative_lower_bound'])/4<=Q(198,1000):
        raise RuntimeError('The original bounded-contrast lower bound failed.')
    (ROOT/'audit/small_size_certificate.json').write_text(json.dumps(out,indent=2))
