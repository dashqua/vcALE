// Initial output quantities 
u = x - X;
u.write();

uw = xw - X;
uw.write();

p = model.pressure();
p.write();

E.write();

if (doPostPro == "activated") {

// Additionnal output quantities
lm.write();
//matJ.write();
v_.write();
solvedW.write();

// integration over domain of E
Info << "Energy of the system:\n intE = " << op.integrateOverDomain(E,V) << endl;

// theoretical energies
//kinE = 0.5 * (lm&lm) / rho.value();
forAll(potE, n) {
  kinE[n] = 0.5 * (lm[n]&lm[n]) / rho.value();// rho.value() * (vHat[n]&vHat[n]); //(lm[n]&lm[n]) / rho.value();
  //kinE[n] += 0.5 * rho.value() * (solvedW[n]&solvedW[n]);
  potE[n] = (0.5*mu.value()*Foam::pow(trueJ[n],-2./3.)*(trueF[n]&&trueF[n])) - (3*mu.value()/2.) + (0.5*kappa.value()*Foam::pow(trueJ[n]-1,2));
}

// Theoretical Energy
theoE = kinE + potE;

// Residual
residualE = mag(E - (kinE + potE));

// integration over domain of residual and theoretical energy
Info << "Theoretical of the system:\n theoE = " << op.integrateOverDomain(theoE,V) << endl;
Info << "Residual of the system:\n residualE = " << op.integrateOverDomain(residualE,V) << endl;


// ALE - CoG
vector sum = vector::zero;
forAll(u, node) { sum += u[node]; } 
vector cog = sum / (sizeof( u ) / sizeof(u[0]));

scalar index = 129;

Info << "time           : " << runTime.value() << nl;
Info << "CoG            : " << cog[0] << " " << cog[1] << " " << cog[2] << nl;
Info << "sample solvedW : " << solvedW[index] << nl;
Info << "sample lm      : " << lm[index]  << nl;
Info << "sample pR      : " << pR[index]  << nl;
Info << "sample tC2     : " << tC2[index] << nl;
Info << "sample wC      : " << wC[index] << nl;
Info << "sample lambdaH : " << lambdaH[index] << nl;
Info << "sample rhsLm2  : " << rhsLm2[index] << nl;
Info << "sample rhsLm3  : " << rhsLm3[index] << nl;
Info << "sample rhsLm   : " << rhsLm[index] << nl;
Info << "sample rhsW    : " << rhsW[index] << nl;
Info << "deltaT         : " << deltaT.value() << nl;

}
