// Initial output quantities 
u = x - X;
u.write();

uw = xw - X;
uw.write();

p = model.pressure();
//forAll(p, n) { p[n] /= aleJ[n]; }
p.write();


// Additionnal output quantities
lm.write();
matJ.write();
v_.write();




// ALE - comparJ: the following field is used to 
// - matJ, that is RK integrated,
// - aleJ, that is given analytically
pointScalarField comparJ
(
	IOobject ("comparJ", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), //IOobject("comparJ", mesh),
	pMesh,
	dimensionedScalar("comparJ", dimensionSet(0,0,0,0,0,0,0), 0.0)
);
comparJ = mag(matJ - aleJ);
comparJ.write();

// ALE - comparF: the following field is used
// - matF, that is RK integrated,
// - aleF, that s given analytically
pointScalarField comparF
(
	IOobject ("comparF", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), //IOobject("comparJ", mesh),
	pMesh,
	dimensionedScalar("comparF", dimensionSet(0,0,0,0,0,0,0), 0.0)
);
comparF = mag(matF - aleF);
comparF.write();



vector sum = vector::zero;

forAll(u, node)
{
	sum += u[node];
}

scalar n = sizeof( u ) / sizeof(u[0]);
Info << "CoG: " << (sum/n)[0] << " " << (sum/n)[1] << " " << (sum/n)[2] << nl;

pointScalarField normF
(
 IOobject("normF", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
 pMesh,
 dimensionedScalar("normF", dimensionSet(0,0,0,0,0,0,0), 0.0)
 );
normF = mag(F);
normF.write();

pointScalarField normMatFandF = mag(matF - F);
normMatFandF.write();


pointScalarField normRHSF   
(
	IOobject("normRHSF", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
	pMesh,
	dimensionedScalar("normRHSF", rhsF.dimensions(), 0.0)
);
normRHSF = mag(rhsF - rhsMatF);
normRHSF.write();






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

Info << "sample lm      : " << lm[500]  << nl << nl;
Info << "sample pR      : " << pR[500]  << nl << nl;
Info << "sample tC2     : " << tC2[500] << nl << nl;
Info << "sample wC      : " << wC[500] << nl << nl;
Info << "sample lambdaH : " << lambdaH[500] << nl;
Info << "sample rhsLm2  : " << rhsLm2[500] << nl;
Info << "sample rhsLm3  : " << rhsLm3[500] << nl;
