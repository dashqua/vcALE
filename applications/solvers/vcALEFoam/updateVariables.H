// Strain measures
stretch = mech.stretch();

// Wave speeds
forAll(mesh.points(), node)
{
    Up_time[node] = model.Up().value()/stretch[node];
}

// Constitutive model
aleModel_.correct();
trueF = F & op.inverse(matF);
/*
forAll(trueF, n)
{
	trueF[n] = Foam::tensor::I;
}
*/
trueJ = det(trueF);
trueH = trueJ * op.invT(trueF);
model.correct(trueF, trueH, trueJ);
P = model.piola();

pR = op.inverseScalar(matJ) * lm;
/*
forAll(pR, n) {
  pR[n][0] = 1000;
  pR[n][1] = 1000;
  pR[n][2] = 0;
}
*/
v_ = pR/rho; 
vHat = v_ + (trueF & w_);
PHat = P & aleH;

// Linear reconstruction
vGrad = grad.gradient(v_);
grad.reconstruct(v_, vGrad, vM, vP);

lambdaH = sqrt((inter.pointToEdge(aleH) & N) & (inter.pointToEdge(aleH) & N)); 
scalarField aleUp___ = inter.pointToEdge( op.inverseScalar(matJ) ) *
       (
         (lambdaH * inter.pointToEdge(Up))
       - (inter.pointToEdge(w_) & (inter.pointToEdge(aleH) & N))
       );

scalarField aleUs___ = inter.pointToEdge( op.inverseScalar(matJ) ) *
       (
         (lambdaH * inter.pointToEdge(Us))
       - (inter.pointToEdge(w_) & (inter.pointToEdge(aleH) & N))
       );

// Average terms
tC = (inter.pointToEdge(PHat) & N)
  + inter.pointToEdge(matJ*rho)*(0.5*mech.SmatrixMod(aleUp___, aleUs___) & (vP - vM));
if (GCLtechnique == "withInvolution") {
    tC2 = inter.pointToEdge( pR )   * N;
} else if (GCLtechnique == "withoutInvolution") {
  HTwC = inter.pointToEdge( aleH.T()&w_ ) & N;
} else {
  FatalErrorIn("updateVariables.H") << "GCL technique is not properly defined." << abort(FatalError);
}

vC  = inter.pointToEdge( vHat ) * N;
wC  = inter.pointToEdge( w_ )   * N;

