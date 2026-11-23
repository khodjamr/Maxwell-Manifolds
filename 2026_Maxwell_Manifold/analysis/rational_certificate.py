"""Optional candidate-certificate generator; supplied proof inputs are unchanged.

New candidate files go to verification_results/generated_certificates/.
Replay the archived proofs with verify_certificates.py or verify_release.py.
"""
from pathlib import Path
from fractions import Fraction as Q
import json
import mpmath as mp
from refine_critical import abc,disc,refine
from rational_interval import I,J,BITS,GRID,N,decimal,interval_system

def certificate(case):
 mp.mp.dps=85
 rho=int(Q(case['rho']));kk=mp.mpf(case['k']);ww=mp.mpf(case['w'])
 centre=[Q(mp.nstr(kk,52)),Q(mp.nstr(ww,52))];radius=[Q(1,10**18),Q(1,10**18)]
 X=[I(c-r,c+r) for c,r in zip(centre,radius)]
 g0,_,_,_=interval_system(I(centre[0]),I(centre[1]),rho)
 _,Jac,ABC,signs=interval_system(X[0],X[1],rho)
 jf=mp.matrix([[mp.diff(lambda x:disc(x,mp.mpf(rho),ww),kk),mp.diff(lambda x:disc(kk,mp.mpf(rho),x),ww)],
 [mp.diff(lambda x:disc(x,mp.mpf(rho),ww),kk,2),mp.diff(lambda x:mp.diff(lambda y:disc(x,mp.mpf(rho),y),ww),kk)]])
 pre=jf**-1;C=[[Q(mp.nstr(pre[i,j],65)) for j in range(2)] for i in range(2)]
 assert C[0][0]*C[1][1]-C[0][1]*C[1][0]!=0
 R=[[I(int(i==j))-sum((I(C[i][h])*Jac[h][j] for h in range(2)),I(0)) for j in range(2)] for i in range(2)]
 dx=[I(-r,r) for r in radius]
 K=[I(centre[i])-sum((I(C[i][h])*g0[h] for h in range(2)),I(0))+sum((R[i][j]*dx[j] for j in range(2)),I(0)) for i in range(2)]
 q=max(sum((R[i][j].absmax() for j in range(2)),Q(0)) for i in range(2))
 # Banach's theorem is applied to the EXACT rational box, not its
 # slightly larger dyadic enclosure used to evaluate the Jacobian.
 inclusion=all(centre[i]-radius[i]<K[i].lo and K[i].hi<centre[i]+radius[i] for i in range(2));assert inclusion and q<1
 A,B,_=ABC;mu=-B/(2*A);eps=X[0]*X[0]/(rho*rho*mu)
 dkk,dw=signs;assert not dkk.lo<=0<=dkk.hi;assert not dw.lo<=0<=dw.hi
 out={'rho':rho,'bits':BITS,'taylor_terms':N+1,'centre_rational':[str(x) for x in centre],'box_radii_rational':[str(x) for x in radius],
 'preconditioner_rational':[[str(x) for x in row] for row in C],
 'k_box':X[0].bounds(),'w_box':X[1].bounds(),'krawczyk_image':[x.bounds() for x in K],
 'contraction_bound':decimal(q),'strict_inclusion':inclusion,'epsilon_box':eps.bounds(),'mu_box':mu.bounds(),
 'A_box':A.bounds(),'D_kk_box':dkk.bounds(),'D_w_box':dw.bounds(),
 'type':'saddle' if dkk.lo>0 else 'definite'}
 print(json.dumps(out,indent=2),flush=True);return out

if __name__=='__main__':
 mp.mp.dps=80
 cases=[refine(3,6.003022377754758,3.458302594955618),refine(3,8.449848767106866,1.0658307471841615)]
 root=Path(__file__).resolve().parents[1]/'verification_results/generated_certificates'
 (root/'data').mkdir(parents=True,exist_ok=True)
 (root/'audit').mkdir(exist_ok=True)
 (root/'data/high_precision_critical_points.json').write_text(json.dumps(cases,indent=2))
 certs=[certificate(c) for c in cases]
 (root/'audit/critical_point_certificates.json').write_text(json.dumps(certs,indent=2))
