vector sum = vector::zero;

forAll(u, node)
{
	sum += u[node];
}

scalar n = sizeof( u ) / sizeof(u[0]);
Info << "CoG: " << sum/n << nl;


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
