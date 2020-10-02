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
model.correct();
P = model.piola();

// update of vHat
vHat = (aleInvJ * lm / rho) + ( F & inv(aleF) & w_);

// Linear reconstruction
lmGrad = grad.gradient(lm);
grad.reconstruct(lm, lmGrad, lmM, lmP);

//ALE model
aleModel_.correct();
aleRecLM = lm * aleJ;


lmC = inter.pointToEdge(lm);
tC = (inter.pointToEdge((P + aleInvJ * (p_ * w_)) & aleH ) & N) + (0.5*mech.Smatrix(aleUp, aleUs) & (lmP - lmM));
