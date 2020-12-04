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
matP    = materialModel_.piola(matF, aleH, matJ);           //P(material)
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
spatJdet = det(F); //spatJ   = det(F);
spatH    = spatJdet * op.invT(F);//spatJ * op.invT(F);
spatP    = spatialModel_.piola(F,  spatH, spatJdet);             //P(spatial)

// Update 'true'/solid quantities
trueF = F & op.inverse(matF);
trueJ = det(trueF);
trueH = trueJ * op.invT(trueF);
if (solidModel__ == "vonMises") {
    D = DH / matJ; alpha = alphaH / matJ; 
    model.correct(trueF, trueH, trueJ, spatJ, matJ, D, alpha);
    forAll(D, n) { DH[n] = D[n] * matJ[n];  alphaH[n] = alpha[n] * matJ[n]; }
} else {
    model.correct(trueF, trueH, trueJ, spatJ, matJ);
}
starP   = ((1.0-alpha_)*matP) + (alpha_*(trueF.T()&spatP));    // P*
P = model.piola();

// Update velocities
pR = op.inverseScalar(matJ) * lm;
v_ = pR/rho; 
vHat = v_ + (trueF & solvedW);
PHat = P & aleH;

// Linear reconstruction
vGrad = grad.gradient(v_);
grad.reconstruct(v_, vGrad, vM, vP);
//
p = model.pressure();
pGrad = grad.gradient(p);
grad.reconstruct(p, pGrad, pM, pP);
//
wGrad = grad.gradient(solvedW);
grad.reconstruct(solvedW, wGrad, wM, wP);

// Average terms
stabC    = inter.pointToEdge(matJ*rho)*(0.5*mech.SmatrixMod(aleUp___, aleUs___) & (vP - vM));
stabJacC = ((0.5*mech.SmatrixJacobianMod(aleUp___, aleUs___)) * (pP - pM))/rho.value();
stabWC   = rho.value() * 0.5 * (mech.SmatrixMod(aleUp___, aleUs___) & (wP - wP));

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
} 
if (solidModel__ == "vonMises") {	
    //D = DH / matJ;    alpha = alpha / matJ;
    forAll(D1, n) {
        D1[n] = vector(D[n].xx(), D[n].xy(), D[n].xz());
        D2[n] = vector(D[n].yx(), D[n].yy(), D[n].yz());
        D3[n] = vector(D[n].zx(), D[n].zy(), D[n].zz());
    }
}
