(* ::Package:: *)

(* ::Section:: *)
(*Intro*)


BeginPackage["SchrodingerSolver`"];


YukawaPhaseShift::usage = 
"YukawaPhaseShift[a,b,L] returns phase shift and error {\!\(\*SubscriptBox[\(\[Delta]\), \(L\)]\),\!\(\*SubscriptBox[\(\[CapitalDelta]\[Delta]\), \(L\)]\)}.\n
YukawaPhaseShift[a,b,L,options] sets options:
\t AccuracyTest->\[Epsilon] requires \!\(\*SubscriptBox[\(\[Delta]\), \(L\)]\) converges to a fractional error < \[Epsilon] (default \[Epsilon]=0.001), 
\t MaxNumIterations->num gives up after N iteractions (default N=10).\n
Sign of potential = sign of b, i.e. V(x) = sign(b)*E^(-abs(b)*x)/x. \n
Note: NDSolve Stiffness and MaxSteps errors are handled despite warnings."


Unprotect[AccuracyTest,VerboseOutput,MaxNumIterations];
Clear[AccuracyTest,VerboseOutput,MaxNumIterations];
Protect[AccuracyTest,VerboseOutput,MaxNumIterations];

Begin["`Private`"];


(* ::Section:: *)
(*Compute phase shift*)


Options[YukawaPhaseShift]={AccuracyTest->0.001,VerboseOutput->False,MaxNumIterations->10};

WaveFunctionSolver[a_,V_,L_,xrange_,initconds_]:=Block[{output,xmin,xmax,\[Chi],\[Chi]sol,xmin2,x,ErrorFlag=False,sol},

	(* begin solve for \[Chi][x] while checking for NDSolve errors *)
	xmin=xrange[[1]];xmax=xrange[[2]];
	Check[
		sol=NDSolve[{\[Chi]''[x]==-(a^2-L (L+1)/x^2-V[x])\[Chi][x],\[Chi][xmin]==initconds[[1]],\[Chi]'[xmin]==initconds[[2]]},\[Chi],{x,xmin,xmax},MaxSteps->10^6][[1]],
		ErrorFlag=True, {NDSolve::ndsz,NDSolve::mxst}
	];

	While[ErrorFlag,
		\[Chi]sol[x_]=\[Chi][x]/.sol; 
		xmin2=sol[[1,2,1,1,1]]+.9(sol[[1,2,1,1,2]]-sol[[1,2,1,1,1]]);

		Check[
			sol=NDSolve[{\[Chi]''[x]==-(a^2-L (L+1)/x^2-V[x])\[Chi][x],\[Chi][xmin2]==\[Chi]sol[xmin2],\[Chi]'[xmin2]==\[Chi]sol'[xmin2]},\[Chi],{x,xmin2,xmax},MaxSteps->10^6][[1]]; ErrorFlag=False;,
			ErrorFlag=True,{NDSolve::ndsz,NDSolve::mxst}
		];
	]; (* End While *)

	(* end solve for \[Chi][x] while checking for NDSolve errors *)

	\[Chi]sol[x_]=\[Chi][x]/.sol;
	output = {\[Chi]sol[xmax], \[Chi]sol'[xmax]};
	
	output
];

PhaseShift[a_,L_,wavefunction_,xfinal_]:=Block[{j,djdx,n,dndx,tan\[Delta]L,ax,\[Chi],d\[Chi]dx},

	\[Chi]=wavefunction[[1]]; d\[Chi]dx=wavefunction[[2]];

	ax=a*xfinal;
	j= Sqrt[\[Pi]/(2ax)]*BesselJ[L+1/2,ax];
	djdx=-(1/2) Sqrt[\[Pi]/2] (1/ax)^(3/2) BesselJ[1/2+L,ax]+1/2 Sqrt[\[Pi]/2] Sqrt[1/ax] (BesselJ[-(1/2)+L,ax]-BesselJ[3/2+L,ax]);
	n=(-1)^(L+1) Sqrt[\[Pi]/(2ax)]BesselJ[-L-1/2,ax];
	dndx=(-1)^L Sqrt[\[Pi]/2] (1/ax)^(3/2) ((1+L) BesselJ[-(1/2)-L,ax]+ax BesselJ[1/2-L,ax]);

	tan\[Delta]L=(\[Chi]*ax*djdx -(xfinal*d\[Chi]dx-\[Chi])*j)/(\[Chi]*ax*dndx-(xfinal*d\[Chi]dx-\[Chi])n);
	ArcTan[tan\[Delta]L]

];

YukawaPhaseShift[a_,b_,L_,OptionsPattern[]]:=Block[{xmin,xmid,xmax,V,initconds,wavefunction,\[Delta]Lold,\[Delta]Lnew,Test,crit,iter},

	(* Define potential *)
	(* Sign[b]= +/- for repulsive/attractive *)
	(* Abs[b] = screening length *)
	V=Function[x, Sign[b]*Exp[-x/Abs[b]]/x ];
	
	(* Initialize range *)
	xmin=(L+1)/a;
		
	(* ---------------------------- *)
	(* First test convergence with xmin *)
	(* ---------------------------- *)
	
	crit=OptionValue[AccuracyTest];
	Test[\[Delta]L1_,\[Delta]L2_]:=Abs[ \[Delta]L1/\[Delta]L2-1] > crit; (* If true, has not converged *)
	
	(* Define intermediate point *)
	xmid=10*xmin;
	initconds={1,(L+1)/xmin};
	wavefunction=WaveFunctionSolver[a,V,L,{xmin,xmid},initconds];
	\[Delta]Lold=0;
	\[Delta]Lnew=PhaseShift[a,L,wavefunction,xmid];
	
	For[ iter=0, Test[\[Delta]Lold,\[Delta]Lnew] && xmin>0.001*Min[(L+1)/a,Abs[b]] && iter < OptionValue[MaxNumIterations], iter++ ,
	
		xmin=0.5*xmin;
		initconds={1,(L+1)/xmin};
		wavefunction=WaveFunctionSolver[a,V,L,{xmin,xmid},initconds];
		\[Delta]Lold=\[Delta]Lnew;
		\[Delta]Lnew=PhaseShift[a,L,wavefunction,xmid];
		If[OptionValue[VerboseOutput],Print[{xmin,\[Delta]Lold,\[Delta]Lnew}]];
	
	]; (* End For *)
	
	(* --------------------------- *)
	(* Next test convergence with xmax *)
	(* --------------------------- *)
	
	xmax=Max[10*Abs[b],xmid];
	initconds=wavefunction; (* from previous step at x=xmid *)
	wavefunction=WaveFunctionSolver[a,V,L,{xmid,xmax},initconds];
	\[Delta]Lold=0;
	\[Delta]Lnew=PhaseShift[a,L,wavefunction,xmax];

	For[ iter=0, Test[\[Delta]Lold,\[Delta]Lnew] && xmax<100*Max[(L+1)/a,Abs[b]] && iter < OptionValue[MaxNumIterations], iter++ ,
	
		xmax=2*xmax;
		wavefunction=WaveFunctionSolver[a,V,L,{xmid,xmax},initconds];
		\[Delta]Lold=\[Delta]Lnew;
		\[Delta]Lnew=PhaseShift[a,L,wavefunction,xmax];
		If[OptionValue[VerboseOutput],Print[{xmax,\[Delta]Lold,\[Delta]Lnew}]];

	]; (* End For *)
	
	{\[Delta]Lnew,Abs[\[Delta]Lnew-\[Delta]Lold]}

];



(* ::Section:: *)
(*End*)


End[];
EndPackage[];
