// Initial output quantities 
u = x - X;
u.write();

uw = xw - X;
uw.write();

p = model.pressure();
p.write();

matE.write();

if (doPostPro == "activated") {

// Additionnal output quantities
lm.write();
matJ.write();
v_.write();

//integration over domain of matE
Info << "Energy of the material system:\n intMatE = " << op.integrateOverDomain(matE,V) << endl;

// ALE - CoG
vector sum = vector::zero;
forAll(u, node) { sum += u[node]; } 
vector cog = sum / (sizeof( u ) / sizeof(u[0]));

scalar index = 129;

Info << "time           : " << runTime.value() << nl;
Info << "CoG            : " << cog[0] << " " << cog[1] << " " << cog[2] << nl;
Info << "sample solvedW : " << solvedW[index] << nl;
Info << "sample wDot_   : " << wDot_[index] << nl;
Info << "sample rhsW    : " << rhsW[index] << nl;
Info << "sample lm      : " << lm[index]  << nl;
Info << "sample pR      : " << pR[index]  << nl;
Info << "sample tC2     : " << tC2[index] << nl;
Info << "sample wC      : " << wC[index] << nl;
Info << "sample lambdaH : " << lambdaH[index] << nl;
Info << "sample rhsLm2  : " << rhsLm2[index] << nl;
Info << "sample rhsLm3  : " << rhsLm3[index] << nl;
Info << "sample rhsW    : " << rhsW[index] << nl;

}
