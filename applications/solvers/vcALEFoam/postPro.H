// Initial output quantities 
u = x - X;
u.write();

uw = xw - X;
uw.write();

usolvedW = xsolvedW - X;
usolvedW.write();

p = model.pressure();
p.write();


if (doPostPro == "activated") {

// Additionnal output quantities
lm.write();
matJ.write();
v_.write();




// test
pointScalarField comparW(
	IOobject("comparW", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
	pMesh,
	dimensionedScalar("comparW", dimensionSet(0,1,-1,0,0,0,0), 0.0)
);
comparW = mag(solvedW - w_);
comparW.write();
 


pointScalarField comparUSolvedW(
	IOobject("comparUSolvedW", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
	pMesh,
	dimensionedScalar("comparUSolvedW", dimensionSet(0,1,0,0,0,0,0), 0.0)
);
comparUSolvedW = mag(usolvedW - uw);
comparUSolvedW.write(); 
//


// ALE - CoG
vector sum = vector::zero;
forAll(u, node)
{
	sum += u[node];
}

scalar n = sizeof( u ) / sizeof(u[0]);
Info << "CoG: " << (sum/n)[0] << " " << (sum/n)[1] << " " << (sum/n)[2] << nl;











pointScalarField normRHSLM1  
(
	IOobject("normRHSLM1", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
	pMesh,
	dimensionedScalar("normRHSLM1", rhsLm.dimensions(), 0.0)
);
normRHSLM1 = mag(rhsLm1);
normRHSLM1.write();


pointScalarField normRHSLM2   
(
	IOobject("normRHSLM2", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
	pMesh,
	dimensionedScalar("normRHSLM2", rhsLm.dimensions(), 0.0)
);
normRHSLM2 = mag(rhsLm2);
normRHSLM2.write();


pointScalarField normRHSLM3   
(
	IOobject("normRHSLM3", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
	pMesh,
	dimensionedScalar("normRHSLM3", rhsLm.dimensions(), 0.0)
);
normRHSLM3 = mag(rhsLm3);
normRHSLM3.write();






pR.write();
pTilde.write();



scalar index = 129;
scalar t = runTime.value();

Info << "time           : " << t << nl;
Info << "sample solvedW : " << solvedW[index] << nl;
Info << "sample w_      : " << w_[index] << nl;
Info << "sample comparW : " << comparW[index] << nl;
Info << "sample lm      : " << lm[index]  << nl;
Info << "sample pR      : " << pR[index]  << nl;
Info << "sample tC2     : " << tC2[index] << nl;
Info << "sample wC      : " << wC[index] << nl;
Info << "sample lambdaH : " << lambdaH[index] << nl;
Info << "sample rhsLm2  : " << rhsLm2[index] << nl;
Info << "sample rhsLm3  : " << rhsLm3[index] << nl;
Info << "sample rhsW    : " << rhsW[index] << nl;

}
