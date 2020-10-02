// Strain measures
H = det(F)*op.invT(F);
J = det(F);
stretch = mech.stretch();

// Wave speeds
forAll(mesh.points(), node)
{
    Up_time[node] = model.Up().value()/stretch[node];
}

// Constitutive model
trueJ = det(trueF);
trueH = trueJ * op.invT(trueF);
model.correct(trueF, trueH, trueJ);
P = model.piola();

// Linear reconstruction
lmGrad = grad.gradient(lm);
grad.reconstruct(lm, lmGrad, lmM, lmP);


lmC = inter.pointToEdge(lm);
//tC = (inter.pointToEdge((P + aleInvJ * (p_ * w_)) & aleH ) & N) + (0.5*mech.Smatrix(aleUp, aleUs) & (lmP - lmM));/////////
tC = (inter.pointToEdge(P) & N)  + (0.5*mech.Smatrix(Up, Us) & (lmP - lmM));
