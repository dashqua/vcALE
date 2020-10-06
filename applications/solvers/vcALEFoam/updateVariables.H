// Strain measures
//H = det(F)*op.invT(F);
//J = det(F);
stretch = mech.stretch();

// Wave speeds
forAll(mesh.points(), node)
{
    Up_time[node] = model.Up().value()/stretch[node];
}

// Constitutive model
Info << "aleF sample: " << aleF[20] << nl;
trueF = F & op.inverse(aleF);
trueJ = det(trueF);
trueH = trueJ * op.invT(trueF);
model.correct(trueF, trueH, trueJ);
aleModel_.correct();
P = model.piola(); // 'trueP' obtained from true F, H and J
PHat = ( P + (  op.inverseScalar(aleJ) *(lm * w_)) ) & aleH;

// Linear reconstruction
lmGrad = grad.gradient(lm);
grad.reconstruct(lm, lmGrad, lmM, lmP);

v_ = op.inverseScalar(aleJ) * lm/rho;
vHat = v_ + ((F & op.inverse(aleF)) & w_);
vC = inter.pointToEdge(vHat); 

//alefinvt and lambdah are updated but we need to update aleF too.....
lambdaH = sqrt((inter.pointToEdge(  op.invT(aleF) ) & N) & (inter.pointToEdge(   op.invT(aleF)  ) & N));
scalarField aleUp___ = inter.pointToEdge( op.inverseScalar(aleJ) ) *
       (
         (lambdaH * inter.pointToEdge(Up))
       - (inter.pointToEdge(w_) & (inter.pointToEdge(aleH) & N))
       );

scalarField aleUs___ = inter.pointToEdge( op.inverseScalar(aleJ) ) *
       (
         (lambdaH * inter.pointToEdge(Us))
       - (inter.pointToEdge(w_) & (inter.pointToEdge(aleH) & N))
       );

tC = (inter.pointToEdge(PHat) & N)  + (0.5*mech.SmatrixMod(aleUp___, aleUs___) & (lmP - lmM));
