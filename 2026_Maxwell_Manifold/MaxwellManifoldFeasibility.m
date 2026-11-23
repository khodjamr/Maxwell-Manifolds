(* ::Package:: *)
(* Mathematica 9 companion, repair v3.2.5 (21 November 2026).
   Precision-aware validated quadrature and explicit function exits.
   Return inside Do can exit the loop rather than the enclosing Module.
   Local Catch/Throw tags below make all function exits unambiguous.
   See https://reference.wolfram.com/language/ref/message/Break/nofunc.html
   Maxwell formulas and study parameters are unchanged.
   Native Mathematica 9 is not available in the repair environment. *)
BeginPackage["MaxwellManifoldFeasibility`"];
SelfTest::usage="SelfTest[] checks quadrature, root solvers, independent Maxwell powers, critical points, source-resource variation, weak loss and the spheroid solver sphere benchmark.";
RunStudy::usage="RunStudy[] recomputes ten figure PDFs and labelled CSV files in figures and data beside this package.";
ChannelCoefficients::usage="ChannelCoefficients[ell,polarization,k,rho] returns {A,B,C}; k is the internal wavenumber, polarization is TE or TM.";
TuningPower::usage="TuningPower[epsilon,mu,rho,w] gives the exact dipole power for TM amplitude 1 and TE amplitude Sqrt[w].";
MaterialGradient::usage="MaterialGradient[epsilon,mu,rho,w] gives the two actual material partial derivatives.";
RefineCriticalPoint::usage="RefineCriticalPoint[rho,kSeed,wSeed] refines Delta=Delta_k=0. This is numerical refinement, not an interval proof.";
RunSphericalStudy::usage="RunSphericalStudy[] recomputes the original six spherical figure/data pairs.";
RunExtendedStudy::usage="RunExtendedStudy[] recomputes spheroid, critical map, design, and weak-loss Figures 1 and 8--10. This can take substantial time.";
RunSpheroidStudy::usage="RunSpheroidStudy[] independently matches both Cauchy traces and generates Figure 1 and its CSV.";
RunCriticalMap::usage="RunCriticalMap[] produces the non-exhaustive numerical critical map in Figure 8.";
RunDesignStudy::usage="RunDesignStudy[] produces fixed-power resource and relative-sensitivity Figure 9.";
SpheroidSelfTest::usage="SpheroidSelfTest[] benchmarks the spheroidal implementation on a sphere.";
SourceResource::usage="SourceResource[epsilon,mu,rho,w] returns the current resource for amplitudes {1,Sqrt[w]}; divide by (1+w) to fix total radiated power.";
RelativeSensitivity::usage="RelativeSensitivity[epsilon,mu,rho,w] uses F/P_rad and the explicitly stated logarithmic-material metric.";
LossyTuningPower::usage="LossyTuningPower[epsilon,mu,rho,w] uses the Hermitian-norm spherical dipole source for complex constitutive values.";
LossySourceResource::usage="LossySourceResource[epsilon,mu,rho,w] gives the complex-material spherical dipole current resource.";
WeakLossSelfTest::usage="WeakLossSelfTest[] checks resource-derived power against independently matched lossy fields and volume work.";
RunWeakLossStudy::usage="RunWeakLossStudy[] computes Figure 10 and its CSV: critical sections with fixed loss and retuned material/composition controls.";
QuadratureSelfTest::usage="QuadratureSelfTest[] checks quadrature nodes, positive weights, polynomial moments, and precision-specific cache entries before field calculations.";
WriteDiagnosticReport::usage="WriteDiagnosticReport[] writes kernel, package, quadrature and most recent SelfTest results to data/Maxwell_diagnostic.txt beside this package.";
PackageVersion::usage="PackageVersion[] returns the installed numerical-package revision.";
RootSolverSelfTest::usage="RootSolverSelfTest[] checks successful bisection and Newton exits and a Maxwell material root before figure generation.";
RunFigure7::usage="RunFigure7[] generates only Figure 7 and its CSV, using guarded precision and a bracketed normalized ray coordinate.";
Figure7SelfTest::usage="Figure7SelfTest[] checks the most sensitive sample on each of the two Figure 7 rays.";
Begin["`Private`"];
$Directory=DirectoryName[$InputFileName];$Precision=50;
$PackageVersion="3.2.5";$SelfTestPassed=False;
$ComputationActive=False;$CurrentStage="Package loaded";
$LastRootChecks={};$LastFailure={};$LastRootState={};$LastNewtonState={};
$Figure7Context={};$LastFigure7Checks={};$LastRefinement={};$DebugCriticalRefinement=False;
$Figure7Generated=False;
$LastQuadratureChecks={};$LastCoreChecks={};$LastCriticalChecks={};
$LastResourceChecks={};$LastWeakLossChecks={};$LastSpheroidChecks={};
$LastWeakLossNormContext={};$LastWeakLossSolverState={};$LastWeakLossCase={};
PackageVersion[]:=$PackageVersion;
val[record_,key_]:=key/.record;
hp[x_]:=N[If[MachineNumberQ[x],Rationalize[x,0],x],$Precision];
space[a_,b_,n_Integer]:=Table[a+(b-a)j/(n-1),{j,0,n-1}];
logspace[a_,b_,n_Integer]:=Table[10^(a+(b-a)j/(n-1)),{j,0,n-1}];
finite[x_]:=NumberQ[x]&&FreeQ[x,Indeterminate|ComplexInfinity|DirectedInfinity];
realQ[x_]:=NumberQ[x]&&TrueQ[Im[x]==0];
materialQ[e_,m_,rho_,w_]:=And@@(realQ/@{e,m,rho,w})&&TrueQ[e m!=0&&rho>0&&w>=0];
(* A failed study exits its outer checked computation immediately. *)
SetAttributes[checkedComputation,HoldAll];
checkedComputation[expr_]:=If[TrueQ[$ComputationActive],expr,
 Block[{$IterationLimit=100000000,$ComputationActive=True,$ComputationFailureTag=Unique["MaxwellFailure$"]},
  Module[{result},result=Catch[expr,$ComputationFailureTag];
   If[result===$Failed,WriteDiagnosticReport[]];result]]];
fail[s_]:=($LastFailure={"Stage"->$CurrentStage,"Message"->s,
 "WorkingDigits"->$Precision,"Figure7Context"->$Figure7Context};
 Print["MANIFOLD ERROR [",$CurrentStage,"]: ",s];
 If[TrueQ[$ComputationActive],Throw[$Failed,$ComputationFailureTag],$Failed]);
internalK[e_,m_,rho_]:=If[TrueQ[e m>0],rho Sqrt[e m],I rho Sqrt[-e m]];

(* ================================================================== *)
(* EXACT RADIAL RATIOS, POWERS, AND ACTUAL MATERIAL GRADIENTS             *)
(* ================================================================== *)
sj[0,x_]:=Sin[x]/x;
sj[1,x_]:=Sin[x]/x^2-Cos[x]/x;
sj[l_Integer,x_]/;l>=2:=sj[l,x]=Expand[(2l-1)sj[l-1,x]/x-sj[l-2,x]];
sy[0,x_]:=-Cos[x]/x;
sy[1,x_]:=-Cos[x]/x^2-Sin[x]/x;
sy[l_Integer,x_]/;l>=2:=sy[l,x]=Expand[(2l-1)sy[l-1,x]/x-sy[l-2,x]];
radial[l_Integer]:=radial[l]=Module[{j=sj[l,z],d,vm,vn},
 d=D[j,z]+j/z;vm=(j^2-sj[l-1,z]sj[l+1,z])/2;vn=vm+j d/z;{j,d,vm,vn}];
ratio[l_Integer,p_String]:=ratio[l,p]=Module[{j,d,vm,vn,u,v,n},
 {j,d,vm,vn}=radial[l];{u,v,n}=If[p=="TE",{j,d,vm},{d,j,vn}];
 {u^2/n,z u v/n,z^2 v^2/n}];
abcExpr[l_Integer,p_String]:=abcExpr[l,p]=Module[{j,y,d,dy,aa,bb,cross,q},
 j=sj[l,r];y=sy[l,r];d=D[j,r]+j/r;dy=D[y,r]+y/r;
 {aa,bb}=If[p=="TE",{j^2+y^2,d^2+dy^2},{d^2+dy^2,j^2+y^2}];
 cross=j d+y dy;q=ratio[l,p];
 {r bb D[q[[1]],z]/(8z),-cross D[q[[2]],z]/(4z),aa D[q[[3]],z]/(8r z)}];
abcDerivative[l_Integer,p_String,n_Integer]:=abcDerivative[l,p,n]=D[abcExpr[l,p],{z,n}];
abc[l_,p_,k_,rho_,n_:0]:=N[abcDerivative[l,p,n]/.{z->hp[k],r->hp[rho]},$Precision];
ChannelCoefficients[l_Integer,p_String,k_?NumericQ,rho_?NumericQ]:=
 If[l<1||!MemberQ[{"TE","TM"},p]||!realQ[rho]||!TrueQ[rho>0]||k==0||!realQ[k^2],
 fail["Require ell>=1, TE or TM, positive real rho, and nonzero k with real k^2."],Re[abc[l,p,k,rho]]];
mix[k_,rho_,w_,n_:0]:=abc[1,"TM",k,rho,n]+w abc[1,"TE",k,rho,n];
f[e_,m_,rho_,w_]:=Module[{k=internalK[e,m,rho],q},q=mix[k,rho,w];Re[q.{m,1,1/m}]];
TuningPower[e_?NumericQ,m_?NumericQ,rho_?NumericQ,w_?NumericQ]:=
 If[!materialQ[e,m,rho,w],fail["Require nonzero real material parameters, rho>0, w>=0."],f[hp[e],hp[m],hp[rho],hp[w]]];
gradient[e_,m_,rho_,w_]:=Module[{k=internalK[e,m,rho],q,qk,fk},
 q=mix[k,rho,w];qk=mix[k,rho,w,1];fk=qk.{m,1,1/m};
 Re[{rho^2 m fk/(2k),q[[1]]-q[[3]]/m^2+rho^2 e fk/(2k)}]];
MaterialGradient[e_?NumericQ,m_?NumericQ,rho_?NumericQ,w_?NumericQ]:=
 If[!materialQ[e,m,rho,w],fail["Require nonzero real material parameters, rho>0, w>=0."],
 gradient[hp[e],hp[m],hp[rho],hp[w]]];
hessian[e_,m_,rho_,w_]:=Module[{k=internalK[e,m,rho],q,qk,qkk,fk,fkk,fkm,ke,km,kee,kmm,kem},
 q=mix[k,rho,w];qk=mix[k,rho,w,1];qkk=mix[k,rho,w,2];
 fk=qk.{m,1,1/m};fkk=qkk.{m,1,1/m};fkm=qk[[1]]-qk[[3]]/m^2;
 ke=rho^2 m/(2k);km=rho^2 e/(2k);kee=-rho^4 m^2/(4k^3);kmm=-rho^4 e^2/(4k^3);kem=rho^2/(4k);
 Re[{{fkk ke^2+fk kee,fkk ke km+fkm ke+fk kem},
 {fkk ke km+fkm ke+fk kem,fkk km^2+2fkm km+2q[[3]]/m^3+fk kmm}}]];
delta[k_,rho_,w_,n_:0,j_:0]:=Module[{q,qk,qkk,qw,qkw,a,b,c,ak,bk,ck,aw,bw,cw,moduleExitTag},Catch[
 q=mix[k,rho,w];{a,b,c}=q;
 If[n==0&&j==0,Throw[Re[b^2-4a c],moduleExitTag]];
 qk=mix[k,rho,w,1];{ak,bk,ck}=qk;
 If[n==1&&j==0,Throw[Re[2b bk-4(ak c+a ck)],moduleExitTag]];
 If[n==2&&j==0,qkk=mix[k,rho,w,2];
 Throw[Re[2(bk^2+b qkk[[2]])-4(qkk[[1]]c+2ak ck+a qkk[[3]])],moduleExitTag]];
 qw=abc[1,"TE",k,rho];{aw,bw,cw}=qw;
 If[n==0&&j==1,Throw[Re[2b bw-4(aw c+a cw)],moduleExitTag]];
 If[n==1&&j==1,qkw=abc[1,"TE",k,rho,1];
 Throw[Re[2(bk bw+b qkw[[2]])-4(qkw[[1]]c+ak cw+aw ck+a qkw[[3]])],moduleExitTag]];
 fail["Unsupported discriminant derivative."],moduleExitTag]];
RefineCriticalPoint[rho0_?NumericQ,k0_?NumericQ,w0_?NumericQ]:=
 checkedComputation[Block[{$Precision=Max[80,$Precision+30],$IterationLimit=10000000},Module[
 {rho=hp[rho0],k=N[If[Precision[k0]===Infinity,k0,Round[k0 10^24]/10^24],$Precision],
  w=N[If[Precision[w0]===Infinity,w0,Round[w0 10^24]/10^24],$Precision],
  g,jac,step,q,m,e,it,converged=False,moduleExitTag},Catch[
 (* k0 and w0 are only starting guesses. Exact decimal representatives
    permit fresh evaluations; no extra accuracy is assigned to physical rho. *)
 Do[g={delta[k,rho,w],delta[k,rho,w,1]};
  jac={{delta[k,rho,w,1],delta[k,rho,w,0,1]},
   {delta[k,rho,w,2],delta[k,rho,w,1,1]}};
  If[!VectorQ[g,finite]||!MatrixQ[jac,finite],
   Throw[fail["Invalid critical-point refinement system."],moduleExitTag]];
  step=LinearSolve[jac,g];
  If[!VectorQ[step,finite],Throw[fail["Critical-point refinement solve failed."],moduleExitTag]];
  {k,w}={k,w}-step;
  If[TrueQ[$DebugCriticalRefinement],Print["Refinement ",it,"; step ",N[Max[Abs[step]],6]]];
  $LastRefinement={"WorkingDigits"->$Precision,"Iteration"->it,"rho"->rho,
   "k"->k,"w"->w,"Residual"->g,"Step"->step};
  If[TrueQ[Max[Abs[step]]<10^-32],converged=True;Break[]],{it,1,30}];
 If[!TrueQ[converged],Throw[fail["Critical-point refinement did not converge."],moduleExitTag]];
 q=mix[k,rho,w];m=-q[[2]]/(2q[[1]]);e=k^2/(rho^2 m);
 {"rho"->rho,"k"->k,"w"->w,"epsilon"->e,"mu"->m,
  "F"->f[e,m,rho,w],"Gradient"->gradient[e,m,rho,w],
  "Hessian"->hessian[e,m,rho,w],"DeltaKK"->delta[k,rho,w,2],
  "DeltaW"->delta[k,rho,w,0,1],"Fw"->(abc[1,"TE",k,rho].{m,1,1/m}),
  "RefinementWorkingDigits"->$Precision},moduleExitTag]]]];
Clear[criticalCaseCache];
criticalCaseCache[digits_Integer]:=criticalCaseCache[digits]=Block[{$Precision=digits},
 (* Exact decimal representatives below are numerical starting guesses,
    independently refined at 100 digits. Maxwell formulas recompute all fields. *)
 {RefineCriticalPoint[3,1200604475550781732433718828631948952213302977087253411242014659353897148801236367637/200000000000000000000000000000000000000000000000000000000000000000000000000000000000,3458302594917308355158240717657493498583797591805610750839394301939515695922556835309/1000000000000000000000000000000000000000000000000000000000000000000000000000000000000],
  RefineCriticalPoint[3,4224924383550896825011620906675011561174325310373944312324221295578687193654808485479/500000000000000000000000000000000000000000000000000000000000000000000000000000000000,1065830747184161970372160348743671439084389855372429069393777596031213347491685940231/1000000000000000000000000000000000000000000000000000000000000000000000000000000000000]}];
cases[]:=criticalCaseCache[$Precision];

(* ================================================================== *)
(* INDEPENDENT MAXWELL BOUNDARY MATCHING AND VOLUME-WORK SELF-TEST        *)
(* ================================================================== *)
(* Quadrature is cached by BOTH order and precision.  Clear the old
   one-argument memoized rules so Get[] also repairs an existing session. *)
Clear[gauss,gl,legendreRule,radialRule,legendrePair,guardedLegendreRule];
$NativeQuadratureAvailable=False;
If[!StringMatchQ[$Version,"*Mathics*"],
 Quiet[Check[Needs["NumericalDifferentialEquationAnalysis`"];
  $NativeQuadratureAvailable=Length[DownValues[
   NumericalDifferentialEquationAnalysis`GaussianQuadratureWeights]]>0,
  $NativeQuadratureAvailable=False]]];
$QuadratureFallbackUsed=False;
$QuadratureBackend=If[TrueQ[$NativeQuadratureAvailable],
 "Wolfram GaussianQuadratureWeights","guarded Legendre recurrence"];

legendrePair[n_Integer,x_]:=Module[{p0=1,p1=x,p2,j,moduleExitTag},Catch[
 If[n==0,Throw[{1,0},moduleExitTag]];
 Do[p2=((2j-1)x p1-(j-1)p0)/j;p0=p1;p1=p2,{j,2,n}];
 {p1,n(x p1-p0)/(x^2-1)},moduleExitTag]];

guardedLegendreRule[n_Integer,prec_Integer]:=Module[
 {work=prec+25+Ceiling[2Log[10,n+1]],tol=10^(-prec-5),pairs={},x,p,dp,dx,
  converged,i,it,wt,moduleExitTag},Catch[
 (* Only positive roots are iterated.  Reflection gives exact symmetry;
    the odd-order central root is inserted as exact zero. *)
 Do[x=N[Cos[Pi(i-1/4)/(n+1/2)],work];converged=False;
  Do[{p,dp}=legendrePair[n,x];
   If[!finite[p]||!finite[dp]||!TrueQ[Abs[dp]>0],Throw[$Failed,moduleExitTag]];
   dx=p/dp;x=x-dx;
   If[!finite[x]||!TrueQ[0<x<1],Throw[$Failed,moduleExitTag]];
   If[TrueQ[Abs[dx]<tol],converged=True;Break[]],{it,1,20}];
  If[!TrueQ[converged],Throw[$Failed,moduleExitTag]];
  {p,dp}=legendrePair[n,x];wt=2/((1-x^2)dp^2);
  pairs=Join[pairs,{{-x,wt},{x,wt}}],{i,1,Floor[n/2]}];
 If[OddQ[n],{p,dp}=legendrePair[n,0];AppendTo[pairs,{0,2/dp^2}]];
 N[SortBy[pairs,First],prec],moduleExitTag]];

validLegendreRuleQ[pairs_,n_,prec_]:=Module[{x,w,err,tol=10^(-Max[12,prec-10]),mom,moduleExitTag},Catch[
 If[!MatrixQ[pairs,finite]||Dimensions[pairs]!={n,2},Throw[False,moduleExitTag]];
 {x,w}=Transpose[pairs];
 If[!And@@(TrueQ[-1<#<1]&/@x)||!And@@(TrueQ[#>0]&/@w),Throw[False,moduleExitTag]];
 If[n>1&&!And@@(TrueQ[#>0]&/@(Rest[x]-Most[x])),Throw[False,moduleExitTag]];
 mom=Table[Abs[If[j==0,Total[w],Total[w x^j]]-If[EvenQ[j],2/(j+1),0]],{j,0,Min[10,2n-1]}];
 err=Max[mom];TrueQ[finite[err]&&err<tol],moduleExitTag]];

legendreRule[n_Integer,prec_Integer]:=Module[{pairs=$Failed,result,moduleExitTag},Catch[
 If[n<1||prec<20,Throw[fail["Require positive quadrature order and at least 20 digits."],moduleExitTag]];
 If[TrueQ[$NativeQuadratureAvailable],
  pairs=Quiet[Check[NumericalDifferentialEquationAnalysis`GaussianQuadratureWeights[
   n,-1,1,prec+20],$Failed]];
  If[MatrixQ[pairs,finite],pairs=N[SortBy[pairs,First],prec]]];
 If[!validLegendreRuleQ[pairs,n,prec],
  $QuadratureFallbackUsed=True;
  If[TrueQ[$NativeQuadratureAvailable],$QuadratureBackend="Wolfram quadrature with guarded fallback"];
  pairs=guardedLegendreRule[n,prec]];
 If[!validLegendreRuleQ[pairs,n,prec],
  Throw[fail["Invalid quadrature rule: order "<>ToString[n]<>" at "<>ToString[prec]<>" digits."],moduleExitTag]];
 result=Transpose[pairs];legendreRule[n,prec]=result;result,moduleExitTag]];

radialRule[n_Integer,prec_Integer]:=Module[{rule,x,result,moduleExitTag},Catch[
 rule=legendreRule[n,prec];If[rule===$Failed,Throw[$Failed,moduleExitTag]];
 x=(1+rule[[1]])/2;result={x,rule[[2]]x^2/2};
 If[!And@@(TrueQ[0<#<1]&/@x)||!And@@(TrueQ[#>0]&/@result[[2]]),
  Throw[fail["Invalid radial quadrature nodes or weights."],moduleExitTag]];
 radialRule[n,prec]=result;result,moduleExitTag]];
gauss[n_Integer]:=radialRule[n,$Precision];
gl[n_Integer]:=legendreRule[n,$Precision];

QuadratureSelfTest[]:=checkedComputation[Block[{$IterationLimit=10000000},Module[
 {rows={},rule,rr,n,prec,momentError,radialError,ok=True,pt},
 Do[{n,prec}=pt;rule=legendreRule[n,prec];rr=radialRule[n,prec];
  If[rule===$Failed||rr===$Failed,ok=False;AppendTo[rows,{n,prec,"FAILED","FAILED"}],
   momentError=Max[Table[Abs[If[j==0,Total[rule[[2]]],Total[rule[[2]]rule[[1]]^j]]-
    If[EvenQ[j],2/(j+1),0]],{j,0,Min[10,2n-1]}]];
   radialError=Abs[Total[rr[[2]]]-1/3];
   ok=ok&&TrueQ[momentError<10^(-prec+10)&&radialError<10^(-prec+10)];
   AppendTo[rows,{n,prec,momentError,radialError}]],
  {pt,{{28,50},{28,35},{37,35}}}];
 $LastQuadratureChecks=rows;
 Print["Quadrature backend: ",$QuadratureBackend];
 Print[Grid[Prepend[rows,{"Order","Digits","Polynomial moment error","Radial mass error"}],Frame->All]];
 Print[If[TrueQ[ok],"QuadratureSelfTest PASS","QuadratureSelfTest FAIL"]];TrueQ[ok]]]];

waves[l_,k_,rr_,out_:False]:=Module[{x=k rr,j,jp,jpp,d,dp},
 j=N[sj[l,z]/.z->x,$Precision];jp=N[D[sj[l,z],z]/.z->x,$Precision];
 If[out,j+=I N[sy[l,z]/.z->x,$Precision];jp+=I N[D[sy[l,z],z]/.z->x,$Precision]];
 jpp=-2jp/x-(1-l(l+1)/x^2)j;d=jp+j/x;dp=jpp+jp/x-j/x^2;
 {{j},{Sqrt[l(l+1)]j/x,d},{rr jp},{Sqrt[l(l+1)]rr(jp/x-j/x^2),rr dp}}];
direct[e_,m_,rho_,l_,p_]:=Module[{k=internalK[e,m,rho],i,j,b,o,U,V,Up,Vp,mat,h,bo,nodes,weights,v,u,uk,ef,nu=0,pow=0,n,rule,moduleExitTag},Catch[
 rule=gauss[28];If[rule===$Failed,Throw[$Failed,moduleExitTag]];
 {i,j}=If[p=="TE",{1,2},{2,1}];b=waves[l,k,hp[1]];o=waves[l,rho,hp[1],True];
 {U,V,Up,Vp}={Last[b[[i]]],Last[b[[j]]],Last[b[[i+2]]],Last[b[[j+2]]]};
 mat={{U,-Last[o[[i]]]},{k V/(I rho m),I Last[o[[j]]]}};
 {h,bo}=LinearSolve[mat,{-I rho m Up/(2k),-(V+k Vp)/(2k)}];{nodes,weights}=rule;
 Do[v=waves[l,k,nodes[[n]]];u=v[[i]];uk=v[[i+2]];ef=I rho m uk/(2k)+h u;
 nu+=weights[[n]]Total[Abs[u]^2];pow-=weights[[n]]Total[ef Conjugate[u]]/2,{n,1,Length[nodes]}];
 {"Reactive"->Im[pow]/Abs[bo]^2,"PowerError"->Abs[Re[pow]/(Abs[bo]^2/(2rho^2))-1]},moduleExitTag]];
CoreSelfTest[]:=checkedComputation[Block[{$IterationLimit=1000000},Module[{tests={},er,mr,rho,q,ref,rel,c,ok,moduleExitTag},Catch[
 Do[{er,mr,rho}=hp/@point;Do[ref=direct[er,mr,rho,1,pol];If[ref===$Failed,Throw[False,moduleExitTag]];q=abc[1,pol,internalK[er,mr,rho],rho].{mr,1,1/mr};
 rel=Abs[Re[q]-val[ref,"Reactive"]]/Max[1,Abs[q]];
 AppendTo[tests,{ToString[{point,pol},InputForm],rel,val[ref,"PowerError"]}],{pol,{"TM","TE"}}],{point,{{-3,-17/10,2/5},{-4,2,3/10}}}];
 (* The volume-work tolerance is consistent with the independent reference
    tests and permits alternate interpreters' complex-arithmetic roundoff.
    Critical-point residuals are tested separately at much tighter tolerance. *)
 c=cases[];ok=Max[Flatten[tests[[All,{2,3}]]]]<10^-11&&
 And@@((Abs[val[#,"F"]]<10^-24&&Norm[val[#,"Gradient"]]<10^-22)&/@c)&&val[c[[1]],"DeltaKK"]>0&&val[c[[2]],"DeltaKK"]<0;
 $LastCoreChecks=tests;$LastCriticalChecks=({val[#,"F"],val[#,"Gradient"]}&/@c);
 Print[Grid[Prepend[tests,{"Case","Relative power error","Real-power balance error"}],Frame->All]];
 Print["Critical residuals: ",({val[#,"F"],val[#,"Gradient"]}&/@c)];
 Print[If[TrueQ[ok],"SelfTest PASS","SelfTest FAIL"]];TrueQ[ok],moduleExitTag]]]];

(* ================================================================== *)
(* ROOTS, EXPORTS, AND SMALL-SIZE FORMULAS                              *)
(* ================================================================== *)
root[fn_,aa_,bb_]:=Module[{a=hp[aa],b=hp[bb],fa,fb,x,fx,it,moduleExitTag},Catch[
 fa=fn[a];fb=fn[b];
 $LastRootState={"Iteration"->0,"Bracket"->{a,b},"EndpointValues"->{fa,fb}};
 If[!finite[fa]||!finite[fb]||!realQ[fa]||!realQ[fb],
  Throw[fail["Nonfinite or nonreal root-bracket value."],moduleExitTag]];
 If[TrueQ[Abs[fa]<10^-32],Throw[a,moduleExitTag]];
 If[TrueQ[Abs[fb]<10^-32],Throw[b,moduleExitTag]];
 If[!TrueQ[fa fb<0],Throw[fail["Root not bracketed."],moduleExitTag]];
 Do[x=(a+b)/2;fx=fn[x];
  $LastRootState={"Iteration"->it,"x"->x,"F"->fx,"BracketWidth"->Abs[b-a],
   "FunctionAccuracy"->Accuracy[fx],"WorkingDigits"->$Precision};
  If[!finite[fx]||!realQ[fx],Throw[fail["Invalid function value during bisection."],moduleExitTag]];
  If[TrueQ[Abs[fx]<10^-32||Abs[b-a]<10^-30],Throw[x,moduleExitTag]];
  If[TrueQ[fa fx<=0],b=x;fb=fx,a=x;fa=fx],{it,1,160}];
 fail["Root failed to converge."],moduleExitTag]];
newton[fn_,df_,x0_]:=Module[{x=hp[x0],fx,dfx,dx,xnext,it,moduleExitTag},Catch[
 Do[fx=fn[x];dfx=df[x];
  $LastNewtonState={"Iteration"->it,"x"->x,"F"->fx,"Derivative"->dfx,
   "FunctionAccuracy"->Accuracy[fx],"WorkingDigits"->$Precision};
  If[!finite[fx]||!finite[dfx]||!TrueQ[Abs[dfx]>0],
   Throw[fail["Invalid function value or zero derivative in Newton iteration."],moduleExitTag]];
  dx=fx/dfx;xnext=x-dx;
  $LastNewtonState=Join[$LastNewtonState,{"Correction"->dx,
   "CorrectionAccuracy"->Accuracy[dx],"StepTolerance"->10^-30}];
  If[!finite[dx]||!finite[xnext],Throw[fail["Invalid Newton step."],moduleExitTag]];
  x=xnext;If[TrueQ[Abs[dx]<10^-30],Throw[x,moduleExitTag]],{it,1,50}];
 fail["Newton iteration failed."],moduleExitTag]];

RootSolverSelfTest[]:=checkedComputation[Block[{$IterationLimit=10000000},Module[
 {b,n,e,m=hp[-12],rho=hp[4/25],w=hp[169/100],e0,rows,ok,moduleExitTag},Catch[
 $CurrentStage="Root-solver checks";
 b=root[Function[xx,xx^2-2],1,2];
 n=newton[Function[xx,xx^2-2],Function[xx,2xx],1];
 If[b===$Failed||n===$Failed,Throw[False,moduleExitTag]];
 e0=estatic[m,w];e=root[Function[xx,f[xx,m,rho,w]],e0-4/5,e0+4/5];
 If[e===$Failed,Throw[False,moduleExitTag]];
 rows={{"Bisection",Abs[b^2-2]},{"Newton",Abs[n^2-2]},
  {"Maxwell material root",Abs[rho^5 f[e,m,rho,w]]}};
 $LastRootChecks=rows;ok=And@@(TrueQ[#[[2]]<10^-25]&/@rows);
 Print[Grid[Prepend[rows,{"Root check","Absolute residual"}],Frame->All]];
 Print[If[TrueQ[ok],"RootSolverSelfTest PASS","RootSolverSelfTest FAIL"]];TrueQ[ok],moduleExitTag]]]];

phi0[e_,m_,w_]:=(e+2)/4-w(m+2)(m+9)/(14m);
phi2[e_,m_,w_]:=-m(e+1)/10-e/4+w(-2e m^2+105m+34e m+223e-1890/m)/1470;
gphi0[e_,m_,w_]:={1/4,-w(1-18/m^2)/14};
gphi2[e_,m_,w_]:={-m/10-1/4+w(-2m^2+34m+223)/1470,-(e+1)/10+w(-4e m+105+34e+1890/m^2)/1470};
estatic[m_,w_]:=-2+(2w/7)(m+11+18/m);
ensureOutputDirectories[]:=Module[{dir},Do[dir=FileNameJoin[{$Directory,folder}];
 If[!DirectoryQ[dir],CreateDirectory[dir]],{folder,{"figures","data"}}]];
write[name_,header_,rows_]:=Module[{path,result,moduleExitTag},Catch[
 ensureOutputDirectories[];path=FileNameJoin[{$Directory,"data",name<>".csv"}];
 result=Export[path,Prepend[N[rows,17],header],"CSV"];
 If[result===$Failed,Throw[fail["CSV export failed: "<>path],moduleExitTag]];result,moduleExitTag]];
save[name_,g_]:=Module[{path,result,moduleExitTag},Catch[
 ensureOutputDirectories[];path=FileNameJoin[{$Directory,"figures",name<>".pdf"}];
 result=Export[path,g,"PDF"];
 If[result===$Failed,Throw[fail["PDF export failed: "<>path],moduleExitTag]];g,moduleExitTag]];
plot[series_,labels_,axes_]:=ListLinePlot[N[series,17],Frame->True,FrameLabel->axes,PlotLegends->labels,
 PlotRange->All,ImageSize->650,BaseStyle->{FontFamily->"Times",FontSize->13}];
branchRows[c_,w_,half_,n_]:=Module[{rho=val[c,"rho"],kc=val[c,"k"],wc=val[c,"w"],rows={},q,d,m,e,k,sgn,i=0,ks,ds,folds={},j},
 ks=space[kc-half,kc+half,n];ds=delta[#,rho,w]&/@ks;
 Do[If[ds[[j]]ds[[j+1]]<0,AppendTo[folds,root[Function[xx,delta[xx,rho,w]],ks[[j]],ks[[j+1]]]]],{j,1,n-1}];
 ks=Sort[Join[ks,folds,{kc}]];
 Do[i++;k=kk;q=mix[k,rho,w];d=q[[2]]^2-4q[[1]]q[[3]];
 If[AnyTrueLocal[folds,Abs[#-k]<10^-25&]||(Abs[w-wc]<10^-30&&Abs[k-kc]<10^-25),d=0];
 Do[If[d>=0,m=(-q[[2]]+sgn Sqrt[d])/(2q[[1]]);e=k^2/(rho^2 m);AppendTo[rows,{sgn,i,k,e,m}]],{sgn,{1,-1}}],
 {kk,ks}];rows];
AnyTrueLocal[x_,fn_]:=Or@@(fn/@x);
branchSegments[rows_,sgn_]:=Module[{s=Select[rows,First[#]==sgn&]},(#[[All,{4,5}]]&/@Split[s,#2[[2]]-#1[[2]]==1&])];

RunSphericalStudy[]:=checkedComputation[Block[{$IterationLimit=10000000},Module[
 {out={},rows,fig,ee,mm,rho,w=hp[169/100],q,v,st,e0,e2,ex,grad,c,kc,wc,ec,mc,dw,k,br,series,curves,
 d,x,y,kscale,dscale,ev,vec,H,fw,h,guess,s,ep,mpa,kap,pr,typ,crit,rr,idx,moduleExitTag},Catch[
 If[!DirectoryQ[FileNameJoin[{$Directory,"figures"}]],CreateDirectory[FileNameJoin[{$Directory,"figures"}]]];
 If[!DirectoryQ[FileNameJoin[{$Directory,"data"}]],CreateDirectory[FileNameJoin[{$Directory,"data"}]]];
 If[!TrueQ[ensureSelfTest[]],Throw[$Failed,moduleExitTag]];crit=cases[];
 (* FIGURE 2 AND CSV: electric and magnetic static dipole profiles. *)
 $CurrentStage="Figure 2";Print["Computing Figure 2."];
 rows={};curves={};Do[rho=hp[rr];series={};
 Do[q=abc[1,"TM",rho Sqrt[ee(-3)],rho].{-3,1,-1/3};v=rho^5 Re[q];st=(ee+2)/4;
 AppendTo[rows,{"TM",rho,ee,-3,v,st}];AppendTo[series,{ee,v}],{ee,hp/@space[-11/2,-11/10,181]}];AppendTo[curves,series],{rr,{1/20,3/20,3/10}}];
 fig=plot[Append[curves,Table[{ee,(ee+2)/4},{ee,space[-11/2,-11/10,181]}]],{"rho=.05","rho=.15","rho=.3","Static"},{"epsilon","rho^5 q (TM)"}];
 curves={};Do[rho=hp[rr];series={};Do[q=abc[1,"TE",rho Sqrt[(-3)mm],rho].{mm,1,1/mm};v=rho^5 Re[q];st=-(mm+2)(mm+9)/(14mm);
 AppendTo[rows,{"TE",rho,-3,mm,v,st}];AppendTo[series,{mm,v}],{mm,hp/@space[-12,-11/10,241]}];AppendTo[curves,series],{rr,{1/20,3/20,3/10}}];
 fig=GraphicsRow[{fig,plot[Append[curves,Table[{mm,-(mm+2)(mm+9)/(14mm)},{mm,space[-12,-11/10,241]}]],{"rho=.05","rho=.15","rho=.3","Static"},{"mu","rho^5 q (TE)"}]},ImageSize->1100];
 write["fig2_static_profiles",{"polarization","rho","epsilon","mu","scaled_exact_q","scaled_static_q"},rows];AppendTo[out,save["fig2_static_profiles",fig]];
 (* FIGURE 3 AND CSV: unique permittivity graph and coordinate fold. *)
 $CurrentStage="Figure 3";Print["Computing Figure 3."];
 rho=hp[4/25];rows={};Do[e0=estatic[mm,w];ex=root[Function[xx,f[xx,mm,rho,w]],e0-4/5,e0+4/5];e2=e0-4rho^2 phi2[e0,mm,w];
 grad=rho^5 gradient[ex,mm,rho,w][[1]];AppendTo[rows,{rho,w,mm,ex,e0,e2,grad}],{mm,hp/@space[-12,-7/5,241]}];
 fig=GraphicsRow[{plot[rows[[All,{3,#}]]&/@{4,5,6},{"Exact","Leading","Second order"},{"mu","Tuned epsilon"}],
 plot[{rows[[All,{3,7}]],({#[[3]],1/4}&/@rows)},{"Exact derivative","Static 1/4"},{"mu","d(rho^5 F)/d epsilon"}]},ImageSize->1100];
 write["fig3_regular_tuning_graph",{"rho","w","mu","epsilon_exact","epsilon_static","epsilon_second_order","scaled_F_epsilon_exact"},rows];AppendTo[out,save["fig3_regular_tuning_graph",fig]];
 (* FIGURE 4 AND CSV: value and gradient errors, orders two and four. *)
 $CurrentStage="Figure 4";Print["Computing Figure 4."];
 rows={};Do[rho=hp[rr];v=rho^5 f[-3,-5/2,rho,w];grad=rho^5 gradient[-3,-5/2,rho,w];
 AppendTo[rows,{rho,Abs[v-phi0[-3,-5/2,w]],Abs[v-phi0[-3,-5/2,w]-rho^2 phi2[-3,-5/2,w]],
 Norm[grad-gphi0[-3,-5/2,w]],Norm[grad-gphi0[-3,-5/2,w]-rho^2 gphi2[-3,-5/2,w]]}],{rr,logspace[-2,-3/4,15]}];
 fig=GraphicsRow[Table[ListLogLogPlot[N[rows[[All,{1,#}]]&/@pair,17],Joined->True,Frame->True,FrameLabel->{"rho","Absolute error"},
 PlotLegends->{"Leading","Through rho^2"},ImageSize->550],{pair,{{2,3},{4,5}}}],ImageSize->1100];
 write["fig4_asymptotic_convergence",{"rho","power_leading_error","power_second_order_error","gradient_leading_error","gradient_second_order_error"},rows];AppendTo[out,save["fig4_asymptotic_convergence",fig]];
 (* FIGURE 5 AND CSV: true saddle crossing and branch reconnection. *)
 $CurrentStage="Figure 5";Print["Computing Figure 5."];
 c=crit[[1]];wc=val[c,"w"];ec=val[c,"epsilon"];mc=val[c,"mu"];rows={};curves={};idx=0;
 Do[br=branchRows[c,wc+dw,hp[3/25],1601];Do[AppendTo[rows,Join[{idx,wc+dw},bb]],{bb,br}];
 series=Join[branchSegments[br,1],branchSegments[br,-1]];
 AppendTo[curves,ListLinePlot[N[series,17],Frame->True,FrameLabel->{"epsilon","mu"},
 PlotRange->{{ec-9/200,ec+9/200},{mc-3/20,mc+3/20}},PlotLabel->Row[{"w-w* = ",N[dw,4]}],ImageSize->400]];idx++,{dw,{-1/50,0,1/50}}];
 fig=GraphicsRow[curves,ImageSize->1100];write["fig5_saddle_reconnection",{"panel","w","branch","sample","internal_k","epsilon","mu"},rows];AppendTo[out,save["fig5_saddle_reconnection",fig]];
 (* FIGURE 6 AND CSV: discriminant maximum and closed tuning loop. *)
 $CurrentStage="Figure 6";Print["Computing Figure 6."];
 c=crit[[2]];wc=val[c,"w"];kc=val[c,"k"];rho=val[c,"rho"];dw=hp[1/10^7];rows={};curves={};
 Do[series={};Do[d=delta[k,rho,wc+sgn dw];AppendTo[series,{1000(k-kc),10^11 d}];AppendTo[rows,{"discriminant",wc+sgn dw,0,0,k,"","",d,"",""}],
 {k,space[kc-1/1250,kc+1/1250,801]}];AppendTo[curves,series],{sgn,{-1,0,1}}];
 fig=plot[curves,{"Below","Critical","Above"},{"1000(k-k*)","10^11 Delta"}];kscale=Sqrt[-2val[c,"DeltaW"]dw/val[c,"DeltaKK"]];dscale=Sqrt[val[c,"DeltaW"]dw];
 br=branchRows[c,wc+dw,hp[1/1250],2401];curves={};Do[series={};
 Do[{idx,k,ep,mpa}=bb[[{2,3,4,5}]];q=mix[k,rho,wc+dw];x=(k-kc)/kscale;y=(2q[[1]]mpa+q[[2]])/dscale;AppendTo[series,{x,y}];
 AppendTo[rows,{"loop",wc+dw,sgn,idx,k,ep,mpa,"",x,y}],{bb,Select[br,First[#]==sgn&]}];AppendTo[curves,series],{sgn,{1,-1}}];
 AppendTo[curves,Table[{Cos[rr],Sin[rr]},{rr,space[0,2Pi,401]}]];
 fig=GraphicsRow[{fig,ListLinePlot[N[curves,17],Frame->True,FrameLabel->{"Normal-form x","Normal-form y"},
 PlotStyle->{Red,Red,Directive[Black,Dashed]},PlotLegends->{"Exact +","Exact -","Quadratic circle"},
 Epilog->{Black,PointSize[.015],Point[{0,0}]},AspectRatio->1,PlotRange->All,ImageSize->650]},ImageSize->1100];
 write["fig6_closed_loop_birth",{"dataset","w","branch","sample","internal_k","epsilon","mu","Delta","normal_form_x","normal_form_y"},rows];AppendTo[out,save["fig6_closed_loop_birth",fig]];
 (* FIGURE 7 AND CSV: independently callable critical-scaling block. *)
 AppendTo[out,RunFigure7[]];
 Print["Six figure PDFs and their CSV tables have been recomputed."];Column[out]
,moduleExitTag]]]];

(* ================================================================== *)
(* FIGURE 7 AND CSV: NORMALIZED, BRACKETED CRITICAL-RAY CALCULATION       *)
(* ================================================================== *)
(* The smaller eigenvalue of a real symmetric 2-by-2 matrix is explicit.
   Scalar sums/square roots retain arbitrary precision on all interpreters.
   The first nonzero component fixes the orientation of the selected ray. *)
rayRoot[fn_,df_,aa_,bb_]:=Module[
 {a=hp[aa],b=hp[bb],fa,fb,x,fx,dfx,xnext,it,moduleExitTag},Catch[
 fa=fn[a];fb=fn[b];
 If[!finite[fa]||!finite[fb]||!TrueQ[fa fb<0],
  Throw[fail["Figure 7 normalized root is not bracketed."],moduleExitTag]];
 x=(a+b)/2;
 Do[fx=fn[x];
  $LastRootState={"Solver"->"Safeguarded normalized ray root","Iteration"->it,
   "x"->x,"F"->fx,"BracketWidth"->Abs[b-a],"FunctionAccuracy"->Accuracy[fx],
   "WorkingDigits"->$Precision};
  If[!finite[fx],Throw[fail["Nonfinite Figure 7 normalized power."],moduleExitTag]];
  If[TrueQ[Abs[fx]<10^-32||Abs[b-a]<10^-30],Throw[x,moduleExitTag]];
  If[TrueQ[fa fx<=0],b=x;fb=fx,a=x;fa=fx];
  dfx=df[x];xnext=If[finite[dfx]&&TrueQ[Abs[dfx]>0],x-fx/dfx,(a+b)/2];
  If[!finite[xnext]||!TrueQ[xnext-a>0]||!TrueQ[b-xnext>0],xnext=(a+b)/2];
  $LastRootState=Join[$LastRootState,{"Derivative"->dfx,"NextX"->xnext}];
  x=xnext,{it,160}];
 fail["Figure 7 normalized root did not converge."],moduleExitTag]];

negativePair[hh_]:=Module[{aa,bb,cc,gap,h,v1,v2,v,norm,err,moduleExitTag},Catch[
 If[Dimensions[hh]!={2,2}||!MatrixQ[hh,realQ],
  Throw[fail["Invalid Figure 7 Hessian."],moduleExitTag]];
 aa=hh[[1,1]];bb=(hh[[1,2]]+hh[[2,1]])/2;cc=hh[[2,2]];
 gap=Sqrt[(aa-cc)^2+4bb^2];h=(aa+cc-gap)/2;
 v1={bb,h-aa};v2={h-cc,bb};
 v=If[TrueQ[Total[v1^2]>=Total[v2^2]],v1,v2];norm=Sqrt[Total[v^2]];
 If[!finite[h]||!TrueQ[h<0]||!TrueQ[norm>0],
  Throw[fail["No resolved negative-curvature ray for Figure 7."],moduleExitTag]];
 v=v/norm;If[TrueQ[v[[1]]<0]||(TrueQ[v[[1]]==0]&&TrueQ[v[[2]]<0]),v=-v];
 err=Sqrt[Total[(hh.v-h v)^2]]/Max[1,Max[Abs[Flatten[hh]]]];
 If[!finite[err]||!TrueQ[err<10^-35],
  Throw[fail["Figure 7 eigenpair precision check failed."],moduleExitTag]];
 {h,v,err},moduleExitTag]];

figure7Rows[smoke_:False]:=Block[{$Precision=80,$IterationLimit=100000000},
 $CurrentStage=If[TrueQ[smoke],"Figure 7 self-test","Figure 7"];
 If[!TrueQ[smoke],Print["Figure 7: preparing guarded-precision critical points."]];Module[
 {crit=cases[],rows={},checks={},c,ec,mc,wc,rho,h,vec,eigenError,fw,typ,
  grid,dw,guess,normalizer,u,dist,ep,mm,g,gn,kap,pred,res,normalizedResidual,idx,j,moduleExitTag},Catch[
 $CurrentStage=If[TrueQ[smoke],"Figure 7 self-test","Figure 7"];
 Do[c=crit[[idx]];{ec,mc,wc,rho}=(val[c,#]&/@{"epsilon","mu","w","rho"});
  {h,vec,eigenError}=negativePair[val[c,"Hessian"]];fw=val[c,"Fw"];
  typ=If[idx==1,"saddle","definite"];
  grid=If[TrueQ[smoke],{If[idx==1,10^-8,10^-10]},
   If[idx==1,logspace[-8,-4,13],logspace[-10,-6,13]]];
  Do[dw=hp[grid[[j]]];normalizer=fw dw;guess=Sqrt[-2normalizer/h];
   $Figure7Context={"Type"->typ,"Sample"->j,"DeltaW"->dw,
    "WorkingDigits"->$Precision,"BaseMaterial"->{ec,mc},"BaseComposition"->wc,
    "PredictedDistance"->guess,"Ray"->vec,"PowerNormalizer"->normalizer};
   If[!TrueQ[normalizer>0]||!TrueQ[guess>0],
    Throw[fail["Invalid Figure 7 normalization."],moduleExitTag]];
   (* u=s/guess is of order one. The endpoints are checked by rayRoot[].
      This solves the same F=0 equation on the same physical ray. *)
   u=rayRoot[Function[uu,f[ec+guess uu vec[[1]],mc+guess uu vec[[2]],rho,wc+dw]/normalizer],
    Function[uu,guess(gradient[ec+guess uu vec[[1]],mc+guess uu vec[[2]],rho,wc+dw].vec)/normalizer],0,2];
   If[!finite[u]||!TrueQ[0<u<2],Throw[fail["Invalid Figure 7 ray root."],moduleExitTag]];
   dist=guess u;ep=ec+dist vec[[1]];mm=mc+dist vec[[2]];
   g=gradient[ep,mm,rho,wc+dw];gn=Sqrt[Total[g^2]];
   If[!VectorQ[g,finite]||!TrueQ[gn>0],
    Throw[fail["Invalid Figure 7 material gradient."],moduleExitTag]];
   res=f[ep,mm,rho,wc+dw];normalizedResidual=Abs[res/normalizer];
   If[!finite[normalizedResidual]||!TrueQ[normalizedResidual<10^-24],
    Throw[fail["Figure 7 normalized residual check failed."],moduleExitTag]];
   kap=1/gn;pred=1/Sqrt[2fw Abs[h]dw];
   AppendTo[rows,{typ,dw,Abs[dist],guess,kap,pred,ep,mm,res}];
   AppendTo[checks,{typ,dw,u,normalizedResidual,eigenError}];$LastFigure7Checks=checks;
   If[!TrueQ[smoke]&&(j==1||Mod[j,4]==0||j==Length[grid]),
    Print["Figure 7 ",typ,": sample ",j," of ",Length[grid]," passed."]],
   {j,Length[grid]}],{idx,2}];rows,moduleExitTag]]];

Figure7SelfTest[]:=checkedComputation[Module[{rows},
 rows=figure7Rows[True];
 Print["Figure7SelfTest ",If[ListQ[rows]&&Length[rows]==2,"PASS","FAIL"]];
 ListQ[rows]&&Length[rows]==2]];

RunFigure7[]:=checkedComputation[Block[{$IterationLimit=100000000},Module[
 {rows,curves,fig,typ,moduleExitTag},Catch[
 $Figure7Generated=False;
 If[!TrueQ[ensureSelfTest[]],Throw[$Failed,moduleExitTag]];
 rows=figure7Rows[False];
 curves=Table[Cases[rows,{typ,d_,s_,g_,kk_,pp_,___}:>{d,s}],{typ,{"saddle","definite"}}];
 curves=Join[curves,Table[Cases[rows,{typ,d_,s_,g_,kk_,pp_,___}:>{d,g}],{typ,{"saddle","definite"}}]];
 fig=ListLogLogPlot[N[curves,17],Joined->True,Frame->True,PlotLegends->{"Saddle","Definite","Saddle prediction","Definite prediction"},
 PlotStyle->{Blue,Orange,Directive[Black,Dashed],Directive[Black,Dashed]},FrameLabel->{"w-w*","Material distance"},ImageSize->550];
 curves=Table[Cases[rows,{typ,d_,s_,g_,kk_,pp_,___}:>{d,kk}],{typ,{"saddle","definite"}}];
 curves=Join[curves,Table[Cases[rows,{typ,d_,s_,g_,kk_,pp_,___}:>{d,pp}],{typ,{"saddle","definite"}}]];
 fig=GraphicsRow[{fig,ListLogLogPlot[N[curves,17],Joined->True,Frame->True,PlotLegends->{"Saddle","Definite","Saddle prediction","Definite prediction"},
 PlotStyle->{Blue,Orange,Directive[Black,Dashed],Directive[Black,Dashed]},FrameLabel->{"w-w*","1/Norm[gradient F]"},ImageSize->550]},ImageSize->1100];
 write["fig7_critical_scaling",{"type","delta_w","material_distance_exact","distance_quadratic_prediction","inverse_gradient_exact","inverse_gradient_quadratic_prediction","epsilon","mu","F_residual"},rows];save["fig7_critical_scaling",fig];$Figure7Generated=True;WriteDiagnosticReport[];fig
,moduleExitTag]]]];

(* ================================================================== *)
(* RESOURCE AND RELATIVE-METRIC SENSITIVITY                             *)
(* ================================================================== *)
resourceRaw[e_,m_,rho_,w_]:=Module[{k,j,d,vm,vn,i0,i2,profiles,U,W,nu,aa,bb,cross=-rho^-5,ans=0,weight},
 (* Real normalized profiles are used in both signs of epsilon*mu.
    Exact outgoing products avoid complex-absolute-value roundoff. *)
 If[e m>0,k=rho Sqrt[e m];{j,d,vm,vn}=N[radial[1]/.z->k,$Precision];
 profiles={{d,k j,vn,1},{j/k,d,vm/k^2,w}},
 k=rho Sqrt[-e m];i0=Sinh[k]/k;j=(k Cosh[k]-Sinh[k])/k^2;d=i0-j/k;i2=i0-3j/k;
 vm=(j^2-i0 i2)/2;vn=j d/k-vm;profiles={{d,-k j,vn,1},{j/k,d,vm/k^2,w}}];
 Do[{U,W,nu,weight}=profiles[[it]];
 {aa,bb}=If[it==1,{rho^-2-rho^-4+rho^-6,rho^-2+rho^-4},{rho^-2+rho^-4,rho^-2-rho^-4+rho^-6}];
 ans+=weight(U^2 bb-2U W cross/(rho m)+W^2 aa/(rho^2 m^2))/(2nu),{it,2}];ans];
SourceResource[e_?NumericQ,m_?NumericQ,rho_?NumericQ,w_?NumericQ]:=
 If[!materialQ[e,m,rho,w],fail["Invalid material or pattern."],resourceRaw[hp[e],hp[m],hp[rho],hp[w]]];
RelativeSensitivity[e_?NumericQ,m_?NumericQ,rho_?NumericQ,w_?NumericQ]:=Module[{g,moduleExitTag},Catch[
 If[!materialQ[e,m,rho,w],Throw[fail["Invalid material or pattern."],moduleExitTag]];
 g=MaterialGradient[e,m,rho,w];(1+w)/(2rho^2 Norm[{e,m}g]),moduleExitTag]];

(* ================================================================== *)
(* FIGURE 1 AND CSV: INDEPENDENT SPHEROID CAUCHY-BOUNDARY MATCHING         *)
(* ================================================================== *)
sGeometry[x_,aspect_]:=Module[{rr,q,nr},rr=(1-x^2+x^2/aspect^2)^(-1/2);q=x Sqrt[1-x^2](1/aspect^2-1)rr^2;
 nr=1/Sqrt[1+q^2];{rr,nr,-q nr}];
sWave[l_,k_,rr_List,xx_List,pol_,out_:False,normalise_:True]:=Module[
 {L=l(l+1),yy,pp,j,jp,jpp,d,dp,zz=k rr,zero,mm,mk,nn,nk,u,uk,v,vk,ss,fac},
 yy=Sqrt[(2l+1)/(4Pi)]LegendreP[l,xx];
 pp=Sqrt[(2l+1)/(4Pi)]l(xx LegendreP[l,xx]-LegendreP[l-1,xx])/(Sqrt[L]Sqrt[1-xx^2]);
 (* Built-in arbitrary-precision spherical Bessel functions avoid the
    severe high-order cancellation of explicit sine/cosine recurrences. *)
 j=N[SphericalBesselJ[l,zz],$Precision];jp=N[SphericalBesselJ[l-1,zz]-(l+1)j/zz,$Precision];
 If[out,j+=I N[SphericalBesselY[l,zz],$Precision];
 jp+=I N[SphericalBesselY[l-1,zz]-(l+1)SphericalBesselY[l,zz]/zz,$Precision]];
 jpp=-2jp/zz-(1-L/zz^2)j;d=jp+j/zz;dp=jpp+jp/zz-j/zz^2;zero=ConstantArray[0,Length[rr]];
 mm={zero,zero,-pp j};mk={zero,zero,-pp rr jp};nn={Sqrt[L]yy j/zz,pp d,zero};
 nk={Sqrt[L]yy rr(jp/zz-j/zz^2),pp rr dp,zero};
 {u,uk,v,vk}=If[pol=="TE",{mm,mk,nn,nk},{nn,nk,mm,mk}];
 ss=If[normalise,If[pol=="TE",l,l-1],0];fac=If[normalise,(2l+1)!!/k^ss,1];
 {fac u,fac(uk-ss u/k)/(2k),fac k v,fac((1-ss)v+k vk)/(2k)}];
sTang[v_,nr_,nt_,pol_]:=If[pol=="TE",v[[3]],nr v[[2]]-nt v[[1]]];
sInner[a_,b_,weights_]:=Total[weights Total[a Conjugate[b]]];
sSolve[rho_,k_,mu_,aspect_,pol_,lmax_,ntrain_:100,ncheck_:173,nrad_:60,nangle_:110]:=Module[
 {orders=Range[1,lmax,2],x,wx,rr,nr,nt,sw,oo,et,ht,other,scales,cols={},u,ut,cu,cut,mat,rhs,cs,scaled,sv,sol,c,h,
 xc,wc,rc,nrc,ntc,swc,eout,hout,ee,hh,er,hr,xa,wa,ur,wr,ra,rv,xv,vw,jv,ef,power,res,l,it,rule,moduleExitTag},Catch[
 rule=gl[ntrain];If[rule===$Failed,Throw[$Failed,moduleExitTag]];{x,wx}=rule;{rr,nr,nt}=sGeometry[x,aspect];sw=2Pi wx rr^2/nr;other=If[pol=="TE","TM","TE"];
 oo=sWave[1,rho,rr,x,pol,True,False];et=sTang[oo[[1]],nr,nt,pol];ht=sTang[oo[[3]]/(I rho),nr,nt,other];
 scales=Sqrt[Total[sw Abs[#]^2]]&/@{et,ht};
 If[!VectorQ[scales,finite]||!And@@(TrueQ[#>0]&/@scales),
  Throw[fail["Invalid spheroid trace normalization."],moduleExitTag]];
 Do[{u,ut,cu,cut}=sWave[l,k,rr,x,pol];AppendTo[cols,{I rho mu sTang[ut,nr,nt,pol],sTang[cut,nr,nt,other]}],{l,orders}];
 Do[{u,ut,cu,cut}=sWave[l,k,rr,x,pol];AppendTo[cols,{sTang[u,nr,nt,pol],sTang[cu/(I rho mu),nr,nt,other]}],{l,orders}];
 mat=Transpose[(Join[Sqrt[sw]#[[1]]/scales[[1]],Sqrt[sw]#[[2]]/scales[[2]]]&/@cols)];
 rhs=Join[Sqrt[sw]et/scales[[1]],Sqrt[sw]ht/scales[[2]]];
 If[!MatrixQ[mat,finite]||!VectorQ[rhs,finite],Throw[fail["Nonfinite spheroid matrix or right-hand side."],moduleExitTag]];
 cs=Norm/@Transpose[mat];
 If[!And@@(TrueQ[#>0]&/@cs),Throw[fail["Zero spheroid matrix column."],moduleExitTag]];
 scaled=Transpose[MapThread[#1/#2&,{Transpose[mat],cs}]];sv=SingularValueList[scaled];
 If[!VectorQ[sv,finite]||Length[sv]!=Length[cs]||!TrueQ[Last[sv]>0],Throw[fail["Spheroid matrix has no resolved full column rank."],moduleExitTag]];
 sol=LeastSquares[scaled,rhs]/cs;
 If[!VectorQ[sol,finite],Throw[fail["Spheroid least-squares solve failed."],moduleExitTag]];
 c=Take[sol,Length[orders]];h=Drop[sol,Length[orders]];
 rule=gl[ncheck];If[rule===$Failed,Throw[$Failed,moduleExitTag]];{xc,wc}=rule;{rc,nrc,ntc}=sGeometry[xc,aspect];swc=2Pi wc rc^2/nrc;
 oo=sWave[1,rho,rc,xc,pol,True,False];eout=sTang[oo[[1]],nrc,ntc,pol];hout=sTang[oo[[3]]/(I rho),nrc,ntc,other];
 ee=ConstantArray[0,Length[xc]];hh=ee;
 Do[{u,ut,cu,cut}=sWave[orders[[it]],k,rc,xc,pol];ee+=sTang[I rho mu c[[it]]ut+h[[it]]u,nrc,ntc,pol];
 hh+=sTang[c[[it]]cut+h[[it]]cu/(I rho mu),nrc,ntc,other],{it,Length[orders]}];
 er=Sqrt[Total[swc Abs[ee-eout]^2]/Total[swc Abs[eout]^2]];hr=Sqrt[Total[swc Abs[hh-hout]^2]/Total[swc Abs[hout]^2]];
 rule=gl[nangle];If[rule===$Failed,Throw[$Failed,moduleExitTag]];{xa,wa}=rule;
 rule=gauss[nrad];If[rule===$Failed,Throw[$Failed,moduleExitTag]];{ur,wr}=rule;ra=First[sGeometry[xa,aspect]];
 rv=Flatten[Outer[Times,ra,ur]];xv=Flatten[Table[ConstantArray[xx,nrad],{xx,xa}]];
 vw=Flatten[Outer[Times,2Pi wa ra^3,wr]];jv=ConstantArray[0,{3,Length[rv]}];ef=jv;
 Do[{u,ut,cu,cut}=sWave[orders[[it]],k,rv,xv,pol];jv+=c[[it]]u;ef+=I rho mu c[[it]]ut+h[[it]]u,{it,Length[orders]}];
 power=-sInner[ef,jv,vw]/2;res=Re[sInner[jv,jv,vw]]/2;
 {"q"->Im[power],"Resource"->res,"J"->jv,"E"->ef,"Weights"->vw,"ETrace"->er,"HTrace"->hr,
 "PowerError"->Abs[Re[power]2rho^2-1],"Condition"->First[sv]/Last[sv]},moduleExitTag]];
SpheroidSelfTest[]:=checkedComputation[Block[{$Precision=35,$IterationLimit=1000000},Module[{ss,ref,errors,moduleExitTag},Catch[
 errors=Table[ss=sSolve[hp[4/5],hp[6/5],hp[8/5],hp[1],pol,3,24,37,18,24];
 If[ss===$Failed,Throw[False,moduleExitTag]];ref=ChannelCoefficients[1,pol,hp[6/5],hp[4/5]].{hp[8/5],1,hp[5/8]};
 Max[Abs[val[ss,"q"]-ref]/Max[1,Abs[ref]],val[ss,"ETrace"],val[ss,"HTrace"]],{pol,{"TM","TE"}}];
 $LastSpheroidChecks=errors;Print["Sphere benchmark for spheroid solver: ",errors];TrueQ[Max[errors]<10^-12],moduleExitTag]]]];
RunSpheroidStudy[]:=checkedComputation[$CurrentStage="Figure 1";Block[{$Precision=35,$IterationLimit=100000000},Module[
 {rho=hp[4/5],k=hp[6/5],aspect=hp[13/10],w=hp[1/10],mus=hp/@{7/10,1,8/5,23/10},rows={},sols,all,qmix,X,abc,qe,errs,jc,jr,ss,ma,mb,err,i,lmax,fig,moduleExitTag},Catch[
 X=({#,1,1/#}&/@mus);
 Do[sols=Table[sSolve[rho,k,mu,aspect,pol,lmax],{pol,{"TM","TE"}},{mu,mus}];
 If[!FreeQ[sols,$Failed],Throw[$Failed,moduleExitTag]];all=Flatten[sols,1];
 qmix=(val[#,"q"]&/@sols[[1]])+w(val[#,"q"]&/@sols[[2]]);abc=LinearSolve[Take[X,3],Take[qmix,3]];
 qe=Abs[X[[4]].abc-qmix[[4]]]/Max[1,Abs[qmix[[4]]]];errs={};ma=mus[[1]];mb=mus[[3]];
 Do[jc=(mb val[ss[[3]],"J"]-ma val[ss[[1]],"J"])/(mb-ma);jr=ma mb(val[ss[[1]],"J"]-val[ss[[3]],"J"])/(mb-ma);
 Do[err=jc+jr/mus[[i]]-val[ss[[i]],"J"];AppendTo[errs,Sqrt[Re[sInner[err,err,val[ss[[i]],"Weights"]]]/(2val[ss[[i]],"Resource"]) ]],{i,4}],{ss,sols}];
 AppendTo[rows,{lmax,rho,k,aspect,w,Max[val[#,"ETrace"]&/@all],Max[val[#,"HTrace"]&/@all],Max[val[#,"PowerError"]&/@all],Max[errs],qe,
 abc[[1]],abc[[2]],abc[[3]],abc[[2]]^2-4abc[[1]]abc[[3]],Max[val[#,"Condition"]&/@all]}];Print["Spheroid order ",lmax,"; trace errors ",Last[rows][[{6,7}]]],{lmax,Range[1,21,2]}];
 write["fig1_spheroid_convergence",{"lmax","rho","k","aspect","w","trace_E","trace_H","real_power_error","affine_current_error","quadratic_power_error","A","B","C","Delta","max_condition"},rows];
 fig=ListLogPlot[Table[rows[[All,{1,j}]],{j,{6,7,9,10}}],Joined->True,Frame->True,FrameLabel->{"Largest interior angular order","Relative error or discrepancy"},
 PlotLegends->{"Electric trace","Magnetic trace","Affine current","Held-out power"},ImageSize->700];save["fig1_spheroid_convergence",fig],moduleExitTag]]]];

(* ================================================================== *)
(* FIGURE 8 AND CSV: NON-EXHAUSTIVE DISCRIMINANT-ROOT CRITICAL MAP         *)
(* ================================================================== *)
weightRoots[k_,rho_]:=Module[{a,b,c,d,e,ff,d0,d1,d2,det,moduleExitTag},Catch[
 {a,b,c}=abc[1,"TM",k,rho];{d,e,ff}=abc[1,"TE",k,rho];
 d0=b^2-4a c;d1=2b e-4(a ff+c d);d2=e^2-4d ff;
 If[d2==0,Throw[If[d1!=0,{-d0/d1,Indeterminate},{Indeterminate,Indeterminate}],moduleExitTag]];
 det=d1^2-4d2 d0;If[det<0,Throw[{Indeterminate,Indeterminate},moduleExitTag]];{(-d1-Sqrt[det])/(2d2),(-d1+Sqrt[det])/(2d2)},moduleExitTag]];
RunCriticalMap[]:=checkedComputation[$CurrentStage="Figure 8";Block[{$IterationLimit=100000000},Module[
 {rhos=hp/@space[1/2,3,21],ks=hp/@space[1/20,14,4001],rows={},curves={},all={},rr,ws,ds,ss,k,w,c,typ,local,j,branch,pol,fig,left,right},
 Do[local={};ws=weightRoots[#,rr]&/@ks;
 Do[ds=Table[w=ws[[j,branch]];If[finite[w]&&w>0,delta[ks[[j]],rr,w,1],Indeterminate],{j,Length[ks]}];
 Do[If[finite[ds[[j]]]&&finite[ds[[j+1]]]&&ds[[j]]ds[[j+1]]<0,
 c=Quiet[Check[RefineCriticalPoint[rr,(ks[[j]]+ks[[j+1]])/2,(ws[[j,branch]]+ws[[j+1,branch]])/2],$Failed]];
 If[ListQ[c],{k,w}=(val[c,#]&/@{"k","w"});
 If[finite[k]&&finite[w]&&1/20<k<14&&w>0&&1/10<=Abs[val[c,"epsilon"]]<=100&&1/10<=Abs[val[c,"mu"]]<=100&&
 Abs[delta[k,rr,w]]<10^-25&&Abs[delta[k,rr,w,1]]<10^-24&&Abs[val[c,"DeltaKK"]]>10^-25&&Abs[val[c,"DeltaW"]]>10^-25&&
 !AnyTrueLocal[local,Abs[#[[2]]-k]<10^-12&&Abs[#[[3]]-w]<10^-10 Max[1,w]&],
 typ=If[val[c,"DeltaKK"]>0,"saddle","definite"];
 AppendTo[local,{rr,k,w,val[c,"epsilon"],val[c,"mu"],val[c,"DeltaKK"],val[c,"DeltaW"],typ,"numerical candidate; no interval certificate"}]]]],{j,Length[ks]-1}],{branch,2}];
 rows=Join[rows,local];Print["rho = ",N[rr,5],"; candidates found = ",Length[local]];
 If[MemberQ[hp/@{1,3/2,2,3},rr],Do[Do[w=ws[[j,branch]];If[finite[w]&&0<w<200,AppendTo[curves,{rr,ks[[j]],w,branch-1}]],{j,Length[ks]}],{branch,2}]],{rr,rhos}];
 write["fig8_critical_map",{"rho","k","w","epsilon","mu","Delta_kk","Delta_w","type","status"},rows];
 write["fig8_discriminant_root_curves",{"rho","k","w","branch"},curves];
 left=ListLogPlot[Table[(#[[{2,3}]]&/@Select[curves,First[#]==rr&]),{rr,hp/@{1,3/2,2,3}}],
 Frame->True,FrameLabel->{"Internal k","Positive discriminant root w"},PlotLegends->{"rho=1","rho=1.5","rho=2","rho=3"},ImageSize->550];
 right=ListLogPlot[Table[(#[[{1,3}]]&/@Select[rows,#[[8]]==typ&]),{typ,{"saddle","definite"}}],
 Frame->True,FrameLabel->{"Exterior size rho","Critical composition w"},PlotLegends->{"Saddle candidate","Definite candidate"},ImageSize->550];
 fig=GraphicsRow[{left,right},ImageSize->1100];save["fig8_critical_map",fig]]]];

(* ================================================================== *)
(* FIGURE 9 AND CSV: FIXED-POWER RESOURCE AND RELATIVE SENSITIVITY         *)
(* ================================================================== *)
loopRows[c_,dw_,n_:161]:=Module[{rho=val[c,"rho"],kc=val[c,"k"],wc=val[c,"w"],w,scale,lo,hi,mid,half,rows={},phase,k,q,d,m,e,ang,sgn,moduleExitTag},Catch[
 w=wc+dw;scale=Sqrt[-2val[c,"DeltaW"]dw/val[c,"DeltaKK"]];
 lo=root[Function[xx,delta[xx,rho,w]],kc-2scale,kc];hi=root[Function[xx,delta[xx,rho,w]],kc,kc+2scale];
 mid=(lo+hi)/2;half=(hi-lo)/2;
 Do[k=mid+half Cos[phase];q=mix[k,rho,w];d=delta[k,rho,w];
 If[d<0&&Abs[d]<10^-28,d=0];If[d<0,Throw[fail["Negative loop discriminant outside the endpoint rounding tolerance."],moduleExitTag]];
 Do[m=(-q[[2]]+sgn Sqrt[d])/(2q[[1]]);e=k^2/(rho^2 m);ang=If[sgn==1,phase,2Pi-phase];AppendTo[rows,{ang,e,m,w}],{sgn,{1,-1}}],
 {phase,space[0,Pi,n]}];SortBy[rows,First],moduleExitTag]];
RunDesignStudy[]:=checkedComputation[$CurrentStage="Figure 9";Block[{$IterationLimit=10000000},Module[
 {c=cases[][[2]],rho,ec,mc,wc,base,rows={},vs,e,m,w,ph,res,kap,dw,matcurves={},resourcecurve={},kapcurve={},fig},
 {rho,ec,mc,wc}=(val[c,#]&/@{"rho","epsilon","mu","w"});base=resourceRaw[ec,mc,rho,wc]/(1+wc);
 Do[vs=loopRows[c,hp[dw]];AppendTo[matcurves,({100(#[[2]]/ec-1),100(#[[3]]/mc-1)}&/@vs)];
 Do[{ph,e,m,w}=v;res=resourceRaw[e,m,rho,w]/((1+w)base);kap=RelativeSensitivity[e,m,rho,w];
 AppendTo[rows,{dw,ph,e,m,100(e/ec-1),100(m/mc-1),res,kap}];
 If[dw==1/100000,AppendTo[resourcecurve,{ph/(2Pi),res}];AppendTo[kapcurve,{ph/(2Pi),kap}]],{v,vs}],{dw,{1/10000000,1/100000}}];
 write["fig9_design_interpretation",{"delta_w","phase_radians","epsilon","mu","relative_epsilon_percent","relative_mu_percent","fixed_power_resource_ratio","relative_metric_inverse_sensitivity"},rows];
 fig=GraphicsRow[{plot[matcurves,{"delta w = 10^-7","delta w = 10^-5"},{"Relative epsilon change (%)","Relative mu change (%)"}],
 plot[{resourcecurve},{"Fixed-power resource"},{"Position around loop","Resource / critical resource"}],
 ListLogPlot[kapcurve,Joined->True,Frame->True,FrameLabel->{"Position around loop","Relative inverse sensitivity"},ImageSize->650]},ImageSize->1400];
 save["fig9_design_interpretation",fig]]]];

(* ================================================================== *)
(* LAUNCHERS: GENERAL-SHAPE CHECK, SPHERICAL STUDY AND EXTENSIONS       *)
(* ================================================================== *)
RunExtendedStudy[]:=checkedComputation[Module[{f1,f8,f9,f10,moduleExitTag},Catch[
 If[!TrueQ[ensureSelfTest[]],Throw[$Failed,moduleExitTag]];
 $CurrentStage="Figure 1";Print["Computing Figure 1: spheroidal consistency check."];f1=RunSpheroidStudy[];
 If[f1===$Failed,Throw[$Failed,moduleExitTag]];
 $CurrentStage="Figure 8";Print["Computing Figure 8: critical-point map."];f8=RunCriticalMap[];
 If[f8===$Failed,Throw[$Failed,moduleExitTag]];
 $CurrentStage="Figure 9";Print["Computing Figure 9: fixed-power resource and sensitivity."];f9=RunDesignStudy[];
 If[f9===$Failed,Throw[$Failed,moduleExitTag]];
 $CurrentStage="Figure 10";Print["Computing Figure 10: weak-loss continuations."];f10=RunWeakLossStudy[];
 If[f10===$Failed,Throw[$Failed,moduleExitTag]];Column[{f1,f8,f9,f10}],moduleExitTag]]];
RunStudy[]:=checkedComputation[Module[{spherical,extended,moduleExitTag},Catch[
 Print["Computing spherical Figures 2--7."];spherical=RunSphericalStudy[];
 If[spherical===$Failed,Throw[$Failed,moduleExitTag]];extended=RunExtendedStudy[];
 If[extended===$Failed,Throw[$Failed,moduleExitTag]];
 Column[{extended[[1,1]],spherical,extended[[1,2]],extended[[1,3]],extended[[1,4]]}],moduleExitTag]]];
RevisionSelfTest[]:=checkedComputation[Block[{$IterationLimit=1000000},Module[{ok,checks={},e,m,rho,w,hh,der,ref},
 ok=CoreSelfTest[];Do[{e,m,rho,w}=hp/@pt;hh=10^-8;
 der=(resourceRaw[e+hh,m,rho,w]-resourceRaw[e-hh,m,rho,w])/(2hh);ref=2rho f[e,m,rho,w];
 AppendTo[checks,Abs[der-ref]/Max[1,Abs[ref]]],{pt,{{2,3,3/10,1/10},{-3,-17/10,2/5,169/100},{-4,2,3/10,1/2}}}];
 Print["Independent resource-variation checks: ",checks];ok=ok&&Max[checks]<10^-12;
 $LastResourceChecks=checks;Print[If[TrueQ[ok],"Extended SelfTest PASS","Extended SelfTest FAIL"]];TrueQ[ok]]]];

(* ================================================================== *)
(* WEAK-LOSS MAXWELL POWERS: HERMITIAN NORM AND CONJUGATE-PROFILE SOURCE   *)
(* ================================================================== *)
lossAbs2[z_]:=Re[z]^2+Im[z]^2;
lossHermitianNorm2[v_]:=Module[{components=Flatten[{v}]},
 If[Length[components]==0||!VectorQ[components,finite],$Failed,
  Total[(Abs[#]^2)&/@components]]];
lossProfile[t_,rr_,pol_]:=Module[{k=Sqrt[t],v,ss,ii,jj},v=waves[1,k,rr];
 {ii,jj,ss}=If[pol=="TE",{1,2,1},{2,1,0}];
 {v[[ii]]/k^ss,(v[[ii+2]]-ss v[[ii]]/k)/(2k^(ss+1)),
 k^(1-ss)v[[jj]],((1-ss)v[[jj]]+k v[[jj+2]])/(2k^(ss+1))}];
lossMaterialQ[e_,m_,rho_,w_]:=NumberQ[e]&&NumberQ[m]&&e m!=0&&realQ[rho]&&rho>0&&realQ[w]&&w>=0;
lossChannel[e_,m_,rho_,pol_,n_:36]:=Module[{t=rho^2 e m,u,ut,cu,cut,U,Ut,W,Wt,aa,bb,ext,nu=0,moment=0,nodes,weights,rv,it,Z,Ze,nue,nui,rule,term,moduleExitTag},Catch[
 $LastWeakLossNormContext={"epsilon"->e,"mu"->m,"rho"->rho,
  "Polarization"->pol,"QuadratureOrder"->n,"InternalIndexSquared"->t,
  "WorkingDigits"->$Precision};
 {u,ut,cu,cut}=lossProfile[t,hp[1],pol];{U,Ut,W,Wt}=Last/@{u,ut,cu,cut};
 ext=waves[1,rho,hp[1],True];{aa,bb}=If[pol=="TE",{Last[ext[[1]]],Last[ext[[2]]]},{Last[ext[[2]]],Last[ext[[1]]]}];
 rule=gauss[n];If[rule===$Failed,Throw[$Failed,moduleExitTag]];{nodes,weights}=rule;Do[rv=lossProfile[t,nodes[[it]],pol];u=rv[[1]];ut=rv[[2]];
 term=lossHermitianNorm2[u];
 If[term===$Failed,$LastWeakLossNormContext=Join[$LastWeakLossNormContext,
   {"NodeIndex"->it,"Profile"->u}];
  Throw[fail["Non-numerical current profile in Hermitian norm."],moduleExitTag]];
 nu+=weights[[it]]term;moment+=weights[[it]]Total[ut Conjugate[u]],{it,n}];
 $LastWeakLossNormContext=Join[$LastWeakLossNormContext,
  {"HermitianNorm"->nu,"NormPrecision"->Precision[nu],"NormAccuracy"->Accuracy[nu]}];
 If[!finite[nu]||!realQ[nu]||!TrueQ[nu>0],
  Throw[fail["Invalid Hermitian current norm; details recorded in diagnostic."],moduleExitTag]];
 Z=U bb-W aa/(rho m);Ze=rho^2 m Ut bb-rho Wt aa;nue=2Re[rho^2 m moment];nui=-2Im[rho^2 m moment];
 {"Reactive"->(Re[Ze Conjugate[Z]]/nu-lossAbs2[Z]nue/(2nu^2))/(2rho),
 "Real"->(-Im[Ze Conjugate[Z]]/nu-lossAbs2[Z]nui/(2nu^2))/(2rho),"Resource"->lossAbs2[Z]/(2nu)},moduleExitTag]];
LossyTuningPower[e_?NumericQ,m_?NumericQ,rho_?NumericQ,w_?NumericQ]:=
 If[!lossMaterialQ[e,m,rho,w],fail["Require nonzero complex material values, real rho>0 and w>=0."],
 val[lossChannel[hp[e],hp[m],hp[rho],"TM"],"Reactive"]+hp[w]val[lossChannel[hp[e],hp[m],hp[rho],"TE"],"Reactive"]];
LossySourceResource[e_?NumericQ,m_?NumericQ,rho_?NumericQ,w_?NumericQ]:=
 If[!lossMaterialQ[e,m,rho,w],fail["Invalid lossy material values."],
 val[lossChannel[hp[e],hp[m],hp[rho],"TM"],"Resource"]+hp[w]val[lossChannel[hp[e],hp[m],hp[rho],"TE"],"Resource"]];
lossDirect[e_,m_,rho_,pol_,n_:28]:=Module[{t=rho^2 e m,tb,dt,u,ut,cu,cut,ub,utb,cub,cutb,U,W,Ud,Wd,ext,aa,bb,c,hh,nn,ww,work=0,en=0,hn=0,jv,ef,hf,dd,cd,it,rv,rule,moduleExitTag},Catch[
 tb=Conjugate[t];dt=tb-t;{u,ut,cu,cut}=lossProfile[t,hp[1],pol];
 {ub,utb,cub,cutb}=lossProfile[tb,hp[1],pol];U=Last[u];W=Last[cu];
 {Ud,Wd}=If[dt==0,{Last[ut],Last[cut]},{Last[(ub-u)/dt],Last[(cub-cu)/dt]}];
 ext=waves[1,rho,hp[1],True];{aa,bb}=If[pol=="TE",{Last[ext[[1]]],Last[ext[[2]]]},{Last[ext[[2]]],Last[ext[[1]]]}];
 {c,hh}=LinearSolve[{{I rho m Ud,U},{Wd,W/(I rho m)}},{aa,-I bb}];rule=gauss[n];If[rule===$Failed,Throw[$Failed,moduleExitTag]];{nn,ww}=rule;
 Do[{u,ut,cu,cut}=lossProfile[t,nn[[it]],pol];{ub,utb,cub,cutb}=lossProfile[tb,nn[[it]],pol];
 {dd,cd}=If[dt==0,{ut,cut},{(ub-u)/dt,(cub-cu)/dt}];
 jv=c ub;ef=I rho m c dd+hh u;hf=c cd+hh cu/(I rho m);
 work-=ww[[it]]Total[ef Conjugate[jv]]/2;en+=ww[[it]]Total[lossAbs2[ef]];hn+=ww[[it]]Total[lossAbs2[hf]],{it,n}];
 {"Reactive"->Im[work],"Real"->Re[work],"PowerBalanceError"->Abs[Re[work]-1/(2rho^2)-rho(Im[e]en+Im[m]hn)/2]},moduleExitTag]];
WeakLossSelfTest[]:=checkedComputation[Block[{$Precision=Max[70,$Precision],$IterationLimit=10000000},Module[{rows={},e,m,rho,p,a,b,ok,moduleExitTag},Catch[
 If[!TrueQ[lossPrecisionSelfTest[]],Throw[False,moduleExitTag]];
 Do[{e,m,rho}=hp/@point;e+=I/100;m+=I/100;
 Do[a=lossChannel[e,m,rho,p,36];b=lossDirect[e,m,rho,p,36];
 If[a===$Failed||b===$Failed,Throw[False,moduleExitTag]];
 AppendTo[rows,{ToString[{point,p},InputForm],Abs[val[a,"Reactive"]-val[b,"Reactive"]],Abs[val[a,"Real"]-val[b,"Real"]],val[b,"PowerBalanceError"]}],{p,{"TM","TE"}}],
 {point,{{-17/10,-23/10,3},{2/5,21,3}}}];
 $LastWeakLossChecks=rows;Print[Grid[Prepend[rows,{"Lossy channel","Reactive discrepancy","Real discrepancy","Power-balance error"}],Frame->All]];
 ok=Max[Flatten[rows[[All,{2,3,4}]]]]<10^-11;Print[If[TrueQ[ok],"WeakLossSelfTest PASS","WeakLossSelfTest FAIL"]];TrueQ[ok],moduleExitTag]]]];

(* ================================================================== *)
(* FIGURE 10 AND CSV: CRITICAL TRANSITIONS IN FIXED-LOSS SECTIONS         *)
(* ================================================================== *)
(* Explicit numerical stencils avoid nested named pure functions. The
   same fourth-order stencil and exact step 1/100000 are retained. *)
lossRealStateQ[v_]:=ListQ[v]&&Length[v]==3&&
 VectorQ[v,Function[entry,finite[entry]&&realQ[entry]]]&&
 TrueQ[v[[2]]!=0&&v[[3]]>=0];
lossStatePower[v_,etae_,etam_,rho_]:=Module[{ee,mm},
 If[!lossRealStateQ[v],fail["Invalid real controls in weak-loss continuation."],
  ee=(v[[1]]/rho^2+etae etam)/v[[2]]+I etae;
  mm=v[[2]]+I etam;LossyTuningPower[ee,mm,rho,v[[3]]]]];
lossSection[v_,etae_,etam_,rho_]:=Module[{h=1/100000,dt={1,0,0},dm={0,1,0}},
 {lossStatePower[v,etae,etam,rho],
  (lossStatePower[v-2h dt,etae,etam,rho]-8lossStatePower[v-h dt,etae,etam,rho]+
   8lossStatePower[v+h dt,etae,etam,rho]-lossStatePower[v+2h dt,etae,etam,rho])/(12h),
  (lossStatePower[v-2h dm,etae,etam,rho]-8lossStatePower[v-h dm,etae,etam,rho]+
   8lossStatePower[v+h dm,etae,etam,rho]-lossStatePower[v+2h dm,etae,etam,rho])/(12h)}];
lossJacobian[v_,etae_,etam_,rho_]:=Module[{h=1/100000,j,dir},
 Transpose[Table[dir=UnitVector[3,j];
  (lossSection[v-2h dir,etae,etam,rho]-8lossSection[v-h dir,etae,etam,rho]+
   8lossSection[v+h dir,etae,etam,rho]-lossSection[v+2h dir,etae,etam,rho])/(12h),{j,3}]]];
(* A Newton iterate is a numerical guess, not a measurement carrying the
   previous correction's error budget. Keep guesses on an exact decimal
   lattice and evaluate the Maxwell functions afresh at each trial. Fixed
   physical inputs are not padded or assigned additional accuracy. *)
lossMinimumAccuracy[v_]:=Min[Accuracy/@Flatten[{v}]];
lossAccuracyUsableQ[v_,places_Integer]:=Module[{a=lossMinimumAccuracy[v]},
 TrueQ[a===Infinity]||(NumberQ[a]&&TrueQ[a>=places+4])];
lossDecimalScalar[x_,places_Integer]:=Module[{a,exponent,count,digits,sgn,moduleExitTag},Catch[
 If[x===0||TrueQ[Precision[x]==0],Throw[0,moduleExitTag]];
 a=Abs[x];exponent=RealDigits[a,10,1][[2]];count=places+exponent+1;
 If[!IntegerQ[exponent]||count<1,Throw[0,moduleExitTag]];
 digits=RealDigits[a,10,count][[1]];
 If[!VectorQ[digits,IntegerQ],Throw[fail["Insufficient digits for a weak-loss numerical guess."],moduleExitTag]];
 sgn=If[TrueQ[x<0],-1,1];
 sgn Quotient[FromDigits[digits]+5,10]/10^places,moduleExitTag]];
lossDecimalVector[v_,places_Integer]:=(lossDecimalScalar[#,places]&/@v);
lossMagnitudeBound[v_,places_Integer]:=Module[{a=lossMinimumAccuracy[v],budget},
 budget=If[a===Infinity,0,10^(-Floor[a])];
 Max[Abs[lossDecimalVector[v,places]]]+10^-places+budget];
lossPrecisionSelfTest[]:=Module[{ok,point},
 point=lossDecimalVector[{1/3,-2/7,3+1/10^35},30];
 ok=!lossAccuracyUsableQ[{0``1.09,0``1.09,0``1.09},30]&&
  lossAccuracyUsableQ[{0``45,0``45,0``45},30]&&
  lossDecimalVector[{0``45,0``45,0``45},30]=={0,0,0}&&
  And@@(TrueQ[Precision[#]===Infinity]&/@point)&&
  Max[Abs[point-{1/3,-2/7,3+1/10^35}]]<=10^-30;
 Print["Weak-loss precision regression ",If[TrueQ[ok],"PASS","FAIL"]];TrueQ[ok]];
lossCritical[c_,etae_,etam_]:=Block[{$Precision=Max[70,$Precision]},Module[
 {rho=hp[val[c,"rho"]],places=30,z,g,jac,rawStep,step,trial,alpha,alpha0,scale,gg,
  accepted,it,j,stepBound,residualBound,candidate,confirmed,moduleExitTag},Catch[
 z=lossDecimalVector[{val[c,"k"]^2,val[c,"mu"],val[c,"w"]},places];
 Do[g=lossSection[z,etae,etam,rho];jac=lossJacobian[z,etae,etam,rho];
  If[!VectorQ[g,realQ]||!VectorQ[g,finite]||!MatrixQ[jac,realQ]||!MatrixQ[jac,finite],
   Throw[fail["Non-numerical weak-loss continuation system."],moduleExitTag]];
  rawStep=LinearSolve[jac,g];
  $LastWeakLossSolverState={"Iteration"->it,"rho"->rho,"eta_e"->etae,"eta_m"->etam,
   "Controls"->z,"ControlDecimalPlaces"->places,"Residual"->g,"Jacobian"->jac,
   "Step"->rawStep,"StepAbsoluteAccuracy"->(Accuracy/@rawStep),"WorkingDigits"->$Precision};
  If[!VectorQ[rawStep,realQ]||!VectorQ[rawStep,finite],
   Throw[fail["Invalid weak-loss Newton correction."],moduleExitTag]];
  (* A zero with low absolute accuracy is not a converged correction. *)
  If[!lossAccuracyUsableQ[rawStep,places],
   Throw[fail["Weak-loss Newton correction has insufficient absolute accuracy."],moduleExitTag]];
  step=lossDecimalVector[rawStep,places];stepBound=lossMagnitudeBound[rawStep,places];
  $LastWeakLossSolverState=Join[$LastWeakLossSolverState,
   {"AppliedDecimalCorrection"->step,"CorrectionMagnitudeEstimate"->stepBound}];
  If[TrueQ[stepBound<10^-18],
   candidate=lossDecimalVector[z-step,places];confirmed=lossSection[candidate,etae,etam,rho];
   If[VectorQ[confirmed,finite]&&VectorQ[confirmed,realQ]&&lossAccuracyUsableQ[confirmed,places],
    residualBound=lossMagnitudeBound[confirmed,places];
    $LastWeakLossSolverState=Join[$LastWeakLossSolverState,
     {"ConvergenceResidual"->confirmed,"ResidualMagnitudeEstimate"->residualBound}];
    If[TrueQ[residualBound<10^-20],Throw[candidate,moduleExitTag]]]];
  scale=Max[Abs[step]/(Max[1,Abs[#]]&/@z)];
  alpha0=If[TrueQ[scale>1/4],1/(4scale),1];accepted=False;
  Do[alpha=alpha0/2^j;trial=lossDecimalVector[z-alpha step,places];
   If[lossRealStateQ[trial],gg=lossSection[trial,etae,etam,rho];
    If[VectorQ[gg,finite]&&VectorQ[gg,realQ]&&TrueQ[Norm[gg]<Norm[g]],
     z=trial;accepted=True;Break[]]],{j,0,16}];
  If[!accepted,Throw[fail["Weak-loss continuation line search failed."],moduleExitTag]],{it,1,60}];
 fail["Weak-loss continuation did not converge."],moduleExitTag]]];
RunWeakLossStudy[]:=checkedComputation[$CurrentStage="Figure 10";Block[{$Precision=Max[70,$Precision],$IterationLimit=100000000},Module[{crit=cases[],rows={},c,typ,kind,eta,ee,mm,z,e,m,w,w0,rho,rate,cs={},ratios={},ss,fig,moduleExitTag},Catch[
 Do[c=crit[[idx]];{rho,w0}=(val[c,#]&/@{"rho","w"});typ=If[idx==1,"saddle","definite"];
 Do[Do[{ee,mm}=If[kind=="electric",{hp[eta],0},{0,hp[eta]}];
 $LastWeakLossCase={"Type"->typ,"Loss"->kind,"eta_e"->ee,"eta_m"->mm};
 Print["Figure 10: ",typ," / ",kind," / loss ",N[eta,4]];
 z=lossCritical[c,ee,mm];If[z===$Failed,Throw[$Failed,moduleExitTag]];
 {m,w}=z[[{2,3}]];e=(z[[1]]/rho^2+ee mm)/m;
 AppendTo[rows,{typ,kind,ee,mm,e,m,w,e-val[c,"epsilon"],m-val[c,"mu"],w-w0,Norm[lossSection[z,ee,mm,rho]]}],
 {eta,{1/100000,1/10000,1/1000}}],{kind,{"electric","magnetic"}}],{idx,2}];
 write["fig10_weak_loss",{"type","loss","eta_e","eta_m","epsilon_real","mu_real","w","epsilon_shift","mu_shift","w_shift","residual"},rows];
 Do[ss=Select[rows,#[[1]]==typ&&#[[2]]==kind&];AppendTo[cs,({#[[3]]+#[[4]],Abs[#[[10]]]}&/@ss)],{typ,{"saddle","definite"}},{kind,{"electric","magnetic"}}];
 Do[c=crit[[idx]];typ=If[idx==1,"saddle","definite"];rho=val[c,"rho"];
 rate=-(1+val[c,"w"]+2rho^5(val[c,"w"]abc[1,"TE",val[c,"k"],rho][[2]]-abc[1,"TM",val[c,"k"],rho][[2]]))/(2rho^2 val[c,"mu"]val[c,"Fw"]);
 ss=Select[rows,#[[1]]==typ&&#[[2]]=="magnetic"&];AppendTo[ratios,({#[[4]],#[[10]]/(rate #[[4]])}&/@ss)],{idx,2}];
 fig=GraphicsRow[{ListLogLogPlot[N[cs,17],Joined->True,Frame->True,PlotLegends->{"Saddle electric","Saddle magnetic","Definite electric","Definite magnetic"},FrameLabel->{"Loss component","Critical composition shift"},ImageSize->650],
 ListLogLinearPlot[N[ratios,17],Joined->True,Frame->True,PlotLegends->{"Saddle","Definite"},FrameLabel->{"Magnetic loss","Shift / first-order prediction"},ImageSize->650]},ImageSize->1300];
 save["fig10_weak_loss",fig],moduleExitTag]]]];

ensureSelfTest[]:=If[TrueQ[$SelfTestPassed],True,SelfTest[]];
WriteDiagnosticReport[]:=Module[{path,report},
 ensureOutputDirectories[];path=FileNameJoin[{$Directory,"data","Maxwell_diagnostic.txt"}];
 report={"PackageVersion"->$PackageVersion,"KernelVersion"->$Version,"System"->$SystemID,
  "WorkingDigits"->$Precision,"Stage"->$CurrentStage,"QuadratureBackend"->$QuadratureBackend,
  "GuardedQuadratureFallbackUsed"->TrueQ[$QuadratureFallbackUsed],
  "QuadratureChecks"->$LastQuadratureChecks,"RootChecks"->$LastRootChecks,"SphericalVolumeChecks"->$LastCoreChecks,
  "CriticalResiduals"->$LastCriticalChecks,"ResourceVariationChecks"->$LastResourceChecks,
  "WeakLossChecks"->$LastWeakLossChecks,"SphereBenchmark"->$LastSpheroidChecks,
  "WeakLossCase"->$LastWeakLossCase,"WeakLossNormContext"->$LastWeakLossNormContext,
  "WeakLossSolverState"->$LastWeakLossSolverState,
  "SelfTestPassed"->TrueQ[$SelfTestPassed],"LastFailure"->$LastFailure,
  "Figure7Generated"->TrueQ[$Figure7Generated],
  "Figure7Outputs"->{FileNameJoin[{$Directory,"figures","fig7_critical_scaling.pdf"}],
   FileNameJoin[{$Directory,"data","fig7_critical_scaling.csv"}]},
  "Figure7Context"->$Figure7Context,"Figure7Checks"->$LastFigure7Checks,
  "LastRootState"->$LastRootState,"LastNewtonState"->$LastNewtonState,
  "LastRefinement"->$LastRefinement};
 Export[path,ToString[report,InputForm],"Text"];Print["Diagnostic report: ",path];path];
fullSelfTest[]:=Block[{$IterationLimit=100000000},Module[{q,r,a,b,c,g7,moduleExitTag},Catch[
 $SelfTestPassed=False;
 Print["Maxwell manifold package ",$PackageVersion,"; kernel ",$Version];
 $CurrentStage="Quadrature checks";q=QuadratureSelfTest[];
 If[!TrueQ[q],Print["SelfTest FAIL: quadrature validation."];Throw[False,moduleExitTag]];
 r=RootSolverSelfTest[];If[!TrueQ[r],Throw[False,moduleExitTag]];
 $CurrentStage="Spherical Maxwell and resource checks";a=RevisionSelfTest[];
 $CurrentStage="Weak-loss checks";b=WeakLossSelfTest[];
 $CurrentStage="Sphere benchmark";c=SpheroidSelfTest[];
 g7=Figure7SelfTest[];
 $SelfTestPassed=TrueQ[a]&&TrueQ[b]&&TrueQ[c]&&TrueQ[g7];
 Print[If[TrueQ[$SelfTestPassed],"SelfTest PASS v"<>$PackageVersion,"SelfTest FAIL v"<>$PackageVersion]];
 TrueQ[$SelfTestPassed],moduleExitTag]]];
SelfTest[]:=Block[{$IterationLimit=100000000},Module[{result},
 $SelfTestPassed=False;$LastFailure={};result=checkedComputation[fullSelfTest[]];
 If[!TrueQ[result],$SelfTestPassed=False;Print["SelfTest did not pass; study is disabled."]];
 WriteDiagnosticReport[];TrueQ[$SelfTestPassed]]];
End[];
EndPackage[];
