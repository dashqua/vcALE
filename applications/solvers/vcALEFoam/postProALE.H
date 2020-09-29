//- Jacobian 
aleJ.write();

//- W
w_.write();

//- v
pointVectorField v_
(
	IOobject
	(
		"v",
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	pMesh,
	dimensionedVector("v_", dimensionSet(0,1,-1,0,0,0,0), Foam::vector::zero)
);
forAll(v_, pt)
{
	v_[pt] = lm[pt] * aleInvJ[pt] * 1./(rho.value());
}
v_.write();
//- vHat

