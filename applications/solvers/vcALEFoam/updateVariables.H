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
trueF = F & aleFinv;
trueJ = det(trueF);
trueH = trueJ * op.invT(trueF);
model.correct(trueF, trueH, trueJ);
P = model.piola(); // 'trueP' obtained from true F, H and J
PHat = ( P + (aleInvJ*(lm * w_)) ) & aleH;

// Linear reconstruction
lmGrad = grad.gradient(lm);
grad.reconstruct(lm, lmGrad, lmM, lmP);

v_ = aleInvJ * lm/rho;
vHat = v_ + ((F & aleFinv) & w_);
vC = inter.pointToEdge(vHat); 

//tC = (inter.pointToEdge((P + aleInvJ * (p_ * w_)) & aleH ) & N) + (0.5*mech.Smatrix(aleUp, aleUs) & (lmP - lmM));/////////
scalarField aleUp___ = inter.pointToEdge(aleInvJ) *
       (
         (lambdaH * inter.pointToEdge(Up))
       - (inter.pointToEdge(w_) & (inter.pointToEdge(aleH) & N))
       );

scalarField aleUs___ = inter.pointToEdge(aleInvJ) *
       (
         (lambdaH * inter.pointToEdge(Us))
       - (inter.pointToEdge(w_) & (inter.pointToEdge(aleH) & N))
       );

tC = (inter.pointToEdge(PHat) & N)  + (0.5*mech.SmatrixMod(aleUp___, aleUs___) & (lmP - lmM));
