// Initial output quantities 
u = x - X;
u.write();

uw = xw - X;
uw.write();

p = model.pressure();
p.write();


// Additionnal output quantities
lm.write();
matJ.write();
v_.write();






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

Info << "sample lm      : " << lm[500]  << nl << nl;
Info << "sample pR      : " << pR[500]  << nl << nl;
Info << "sample tC2     : " << tC2[500] << nl << nl;
Info << "sample wC      : " << wC[500] << nl << nl;
Info << "sample lambdaH : " << lambdaH[500] << nl;
Info << "sample rhsLm2  : " << rhsLm2[500] << nl;
Info << "sample rhsLm3  : " << rhsLm3[500] << nl;
