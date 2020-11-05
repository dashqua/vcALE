// F = F, Cp = C, epsilon_p = eps

// track F, Cp and epsilon

// Given F_nplus1, Cm1_pn, Eps_pn

scalar deltaGamma = 0, nuplus1_Alphan = 0;     // step 2
pointScalarField J_nplus1 = det(F_nplus1);      // step 3
pointScalarField p = kappa_.value() * Foam::log(J_nplus1) / J_nplus1;  //step 4
// this is final pressure for this timestep






pointVectorField bTrial_enplus1 = (F_nplus1 & Cm1_pn) & F_nplus1.T(); // step 5
//spectral decomposition of bTrial => scalarList lambdaTrial ??  AND pVectorFIeld nTrial//step 6
//
// based on 
//         op.eigenStructure(C_[nodeID]);
//         vector eigVal_ = op.eigenValue();
//         tensor eigVec_ = op.eigenVector();
//
op.eigenStructure( bTrial_enplus1 );
vector eigVal_ = op.eigenValue();
tensor eigVec_ = op.eigenVector();  // step 6
//
/*
pointVectorField nAlpha = nTrial; //step 7 
scalar tauPrimeTrial = 2 * mu_.value() * Foam::log(lambdaTrial) - ((2./3.) * mu_.value()*Foam::log(J)); // step 8
if ( f(tauPrimeTrial, eps) > 0  )  {  // define f
	nuAlpha =  ;// step 9
	deltaGamma =  ; // step 10
}
lambdaNE = Foam::exp(Foam::log(lambdaTrial) - deltaGamma*nuAlpha);// step 11
tauPrime = (1 - ())*tauPrimeTrial;// step 12
tauAA = tauPrime + (Joutside * pOutside);// step 13
tau = //step 13 bis
P  = tau * Foutside.invT(); // step 14
*/
