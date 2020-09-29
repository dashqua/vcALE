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

// Linear reconstruction
lmGrad = grad.gradient(lm);
grad.reconstruct(lm, lmGrad, lmM, lmP);

// Contact fluxes
//Info << "  Contact fluxes" << nl;
//Info << " dim P: " << P.dimensions() << nl;
//Info << " dim aleInvJ: " << aleInvJ.dimensions() << nl;
//Info << " dim w_: " << w_.dimensions() << nl;
//Info << " dim p_: " << p_.dimensions() << nl;
//Info << " dim N: " << "no numerical dimensions" << nl; //N.dimensions();

//ALE model
aleModel_.correct();
aleRecLM = lm * aleJ;

lmC = inter.pointToEdge(lm);
tC = (inter.pointToEdge((P + aleInvJ * (p_ * w_)) & aleH ) & N) + (0.5*mech.Smatrix(aleUp, aleUs) & (lmP - lmM));
