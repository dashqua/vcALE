Info << "Post Processing ..." << nl;
// Initial output quantities 
u = x - X;
u.write();

uw = xw - X;
uw.write();

p = model.pressure();
p.write();

E.write();
lm.write();
////solvedW.write();

// print Von Mises quantities if relevant
//model.writeVars();

DH.write();
alphaH.write();


if (doPostPro == "activated") {


Info << "sample L2 of Forces: " << op.L2normXYZ(Forces, V) << endl;


  if (tutorial == "swingingCube") {
    // Compute linear displacements or lm for small amplitudes
    const scalar& A = readScalar(runParameters.lookup("A"));
    const scalar& B = readScalar(runParameters.lookup("B"));
    const scalar& C = readScalar(runParameters.lookup("C"));
    const scalar pi = Foam::constant::mathematical::pi;

    const dimensionedScalar& E = mechanicalProperties.lookup("E");
    const dimensionedScalar& nu = mechanicalProperties.lookup("nu");
    const dimensionedScalar mu = E/(2*(1+nu));
    const dimensionedScalar lambda = (nu*E/((1+nu)*(1 - 2*nu)));
    
    const scalar cd = Foam::sqrt(((lambda.value())+(2*mu.value()))/rho.value());
    const scalar t = runTime.value();
    pointScalarField U0 (IOobject("U0", mesh), pMesh, dimensionedScalar(runParameters.lookup("U0")));
    pointVectorField lm_theo (
			     IOobject("lm_theo", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), pMesh,
			     dimensionedVector("lm_theo", u.dimensions(), vector::zero)
			     );
    forAll(lm_theo, n){
      const scalar X = mesh.points()[n][0];
      const scalar Y = mesh.points()[n][1];
      const scalar Z = mesh.points()[n][2];
    
      lm_theo[n] = - (rho.value()*cd*pi*U0[n]*Foam::sqrt(3.0)/2.) * Foam::sin((Foam::sqrt(3.0)/2.)*cd*pi*t) * vector(
				A * Foam::sin(pi*X/2.) * Foam::cos(pi*Y/2.) * Foam::cos(pi*Z/2.),
			        B * Foam::cos(pi*X/2.) * Foam::sin(pi*Y/2.) * Foam::cos(pi*Z/2.),
				C * Foam::cos(pi*X/2.) * Foam::cos(pi*Y/2.) * Foam::sin(pi*Z/2.)	  
			 );
    }
    lm_theo.write();
    //Compute the error and print it
    Info << "L2-error of LM:\n L2errorLM = " << op.L2diffXYZ(lm, lm_theo, V) << nl;
    lm_theo.write();
  }


  
// Additionnal output quantities
/*
// integration over domain of E
Info << "Energy of the system:\n intE = " << op.integrateOverDomain(E,V) << endl;

// theoretical energies
forAll(potE, n) {
  kinE[n] = 0.5 * (lm[n]&lm[n]) / rho.value();// rho.value() * (vHat[n]&vHat[n]); //(lm[n]&lm[n]) / rho.value();
  potE[n] = 0.5*mu.value()*(Foam::pow(trueJ[n],(-2.0/3.0))*(trueF[n]&&trueF[n])-3.0) + (0.5*kappa.value()*(trueJ[n]-1.0)*(trueJ[n]-1.0));
	  //(0.5*mu.value()*Foam::pow(trueJ[n],-2./3.)*(trueF[n]&&trueF[n])) - (3*mu.value()/2.) + (0.5*kappa.value()*Foam::pow(trueJ[n]-1,2));
}

// Theoretical Energy
theoE = kinE + potE;

// integration over domain of residual and theoretical energy
Info << "Theoretical kinetic energy of the system:\n kinE = " << op.integrateOverDomain(kinE,V) << endl;
Info << "Theoretical potential energy of the system:\n potE = " << op.integrateOverDomain(potE,V) << nl;
Info << "Theoretical energy of the system:\n theoE = " << op.integrateOverDomain(theoE,V) << endl;
*/
//spatJdet = det(F);
//
// ALE - CoG
/*
vector sum = vector::zero;
forAll(u, node) { sum += u[node]; } 
vector cog = sum / (sizeof( u ) / sizeof(u[0]));
*/


/*
scalar index = 129;

Info << "time               : " << runTime.value() << nl;
//Info << "CoG                : " << cog[0] << " " << cog[1] << " " << cog[2] << nl;
Info << "sample solvedW     : " << solvedW[index] << nl;
Info << "sample lm          : " << lm[index]  << nl;
Info << "sample pR          : " << pR[index]  << nl;
Info << "sample tC2         : " << tC2[index] << nl;
Info << "sample lambdaH     : " << lambdaH[index] << nl;
Info << "sample rhsLm1      : " << rhsLm1[index] << nl;
Info << "sample rhsLm2      : " << rhsLm2[index] << nl;
Info << "sample rhsLm3      : " << rhsLm3[index] << nl;
Info << "sample rhsLm4      : " << rhsLm4[index] << nl;
Info << "sample rhsW        : " << rhsW[index] << nl;
Info << "sample rhsSpatJ1   : " << rhsSpatJ1[index] << nl;
Info << "sample rhsSpatJ2   : " << rhsSpatJ2[index] << nl;
Info << "sample rhsSpatJ3   : " << rhsSpatJ3[index] << nl;
Info << "sample spatJdet    : " << spatJdet[index] << nl;
Info << "sample spatJ       : " << spatJ[index] << nl;
//(rhsSpatJ1 + rhsSpatJ2 + rhsSpatJ3
Info << "deltaT         : " << deltaT.value() << nl;
*/
}
