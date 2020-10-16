// Strain measures
stretch = mech.stretch();

// Wave speeds
forAll(mesh.points(), node)
{
    Up_time[node] = model.Up().value()/stretch[node];
}

// Constitutive model
Info << "aleF sample: " << aleF[20] << nl;

aleModel_.correct();
trueF = F & op.inverse(aleF);
trueJ = det(trueF);
trueH = trueJ * op.invT(trueF);
model.correct(trueF, trueH, trueJ);
P = model.piola(); // 'trueP' obtained from true F, H and J
v_ = op.inverseScalar(aleJ) * lm/rho;
vHat = v_ + ((F & op.inverse(aleF)) & w_);
PHat = (P + ((v_*rho) * w_)) & aleH;

// Linear reconstruction
lmGrad = grad.gradient(lm);
grad.reconstruct(lm, lmGrad, lmM, lmP);
// /!\ this computation of left and right terms lmM and lmP 
//     has to be replaced by the spatial equivalent 
//     (see formula).
//     -> Since this lm is now the spatial one, there is no 
//     problem :ok:
vC = inter.pointToEdge(vHat); 

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

tC = (inter.pointToEdge(PHat) & N)  + (0.5*mech.SmatrixMod(aleUp___, aleUs___) & (lmP - lmM));
