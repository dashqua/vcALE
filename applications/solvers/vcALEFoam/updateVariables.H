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

//forAll(trueF, n)
//{
//	trueF[n] = Foam::tensor::I;
//}
//
//trueF = F & op.inverse(aleF);
trueJ = det(trueF);
trueH = trueJ * op.invT(trueF);
model.correct(trueF, trueH, trueJ);
P = model.piola();

v_ = op.inverseScalar(aleJ) * lm/rho;
//forAll(v_, n)
//{
//	v_[n][0] = 1000./1100.;
//	v_[n][1] = 1000./1100.;
//	v_[n][2] = 0;
//}

vHat = v_ + (trueF & w_);
//vHat = v_ + ((F & op.inverse(matF)) & w_);
//vHat = v_ + w_;
//vHat = v_ + ((F & op.inverse(aleF)) & w_);
PHat = (P + ((v_*rho) * w_)) & aleH;

// Linear reconstruction
lmGrad = grad.gradient(v_);
grad.reconstruct(v_, lmGrad, lmM, lmP);

vC = inter.pointToEdge(vHat);
wC = inter.pointToEdge(w_);

lambdaH = sqrt((inter.pointToEdge(aleH) & N) & (inter.pointToEdge(aleH) & N));
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

tC = (inter.pointToEdge(PHat) & N) + inter.pointToEdge(aleJ*rho)*(0.5*mech.SmatrixMod(aleUp___, aleUs___) & (lmP - lmM));
