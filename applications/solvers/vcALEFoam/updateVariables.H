// Strain measures
stretch = mech.stretch();

// Wave speeds
forAll(mesh.points(), node)
{
    Up_time[node] = model.Up().value()/stretch[node];
}

// Update material quantities
materialModel_.correct();
aleH    = matJ * op.invT(matF); 
matP    = materialModel_.piola(matF, aleH, matJ);
//matPres = materialModel_.getMaterialPressure(matF, aleH, matJ);
lambdaH = sqrt((inter.pointToEdge(aleH) & N) & (inter.pointToEdge(aleH) & N)); 
scalarField aleUp___ = inter.pointToEdge( op.inverseScalar(matJ) ) *
       (
         (lambdaH * inter.pointToEdge(Up))
       - (inter.pointToEdge(solvedW) & (inter.pointToEdge(aleH) & N))
       );

scalarField aleUs___ = inter.pointToEdge( op.inverseScalar(matJ) ) *
       (
         (lambdaH * inter.pointToEdge(Us))
       - (inter.pointToEdge(solvedW) & (inter.pointToEdge(aleH) & N))
       );

// Update spatial quantities
spatialModel_.correct();
spatJ   = det(F);
spatH   = spatJ * op.invT(F);
spatP   = spatialModel_.piola(F,  spatH, spatJ);

// Update 'true'/solid quantities
starP   = ((1.0-alpha)*matP) + (alpha*(trueF.T()&spatP));//(alpha*matP) + ((1.0-alpha)*(trueF.T()&spatP)); 
trueF = F & op.inverse(matF);
trueJ = det(trueF);
trueH = trueJ * op.invT(trueF);
model.correct(trueF, trueH, trueJ);
P = model.piola();

// Update velocities
pR = op.inverseScalar(matJ) * lm;
v_ = pR/rho; 
vHat = v_ + (trueF & solvedW);
PHat = P & aleH;

// Linear reconstruction
vGrad = grad.gradient(v_);
grad.reconstruct(v_, vGrad, vM, vP);


// Average terms
stabC = inter.pointToEdge(matJ*rho)*(0.5*mech.SmatrixMod(aleUp___, aleUs___) & (vP - vM));
if (PHatTerm == "withInvolution") {
  forAll(Px, n) { 
    Px[n] = vector(P[n].xx(),P[n].yx(),P[n].zx()); 
    Py[n] = vector(P[n].xy(),P[n].yy(),P[n].zy());
    Pz[n] = vector(P[n].xz(),P[n].yz(),P[n].zz());
  }
  P1C = inter.pointToEdge(Px) * N;
  P2C = inter.pointToEdge(Py) * N;
  P3C = inter.pointToEdge(Pz) * N;
} else if (PHatTerm == "withoutInvolution") {
  tC = (inter.pointToEdge(PHat) & N);
}
if (GCLtechnique == "withInvolution") {
    tC2   = inter.pointToEdge( pR )   * N;
} else if (GCLtechnique == "withoutInvolution") {
    HTwC  = inter.pointToEdge( aleH.T() & solvedW ) & N;
    pRwHC = inter.pointToEdge( (pR*solvedW) & aleH ) & N;					      
} else {
    FatalErrorIn("updateVariables.H") << "GCL technique is not properly defined." << abort(FatalError);
}

vC  = inter.pointToEdge( vHat )      * N;
wC  = inter.pointToEdge( solvedW )   * N;

PmatHC    = inter.pointToEdge( matP  & aleH ) & N;
PstarHC   = inter.pointToEdge( starP & aleH ) & N;

Pmat_pmat = inter.pointToEdge( matP.T()&solvedW ) & N;
