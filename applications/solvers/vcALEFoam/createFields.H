/*
    ____    ____  ________   ______   ____  ____  
   |_   \  /   _||_   __  |.' ____ \ |_   ||   _| 
     |   \/   |    | |_ \_|| (___ \_|  | |__| |   
     | |\  /| |    |  _| _  _.____`.   |  __  |   
    _| |_\/_| |_  _| |__/ || \____) | _| |  | |_  
   |_____||_____||________| \______.'|____||____| 
                                                     */
Info << "Creating mesh parameters ..." << endl;

// Operations class
operations op(mesh);

// Dual mesh class
dualMesh dual(mesh);

// Interpolation schemes class
interpolationSchemes inter(mesh);

// Point mesh
pointMesh pMesh(mesh);

// Dual control volume
const pointScalarField& V = dual.volume();

// Primal mesh material face area normal vector
const surfaceVectorField& Sf_p = mesh.Sf();

// Primal mesh material face area
const surfaceScalarField& magSf_p = mesh.magSf();

// Primal mesh material normal
const surfaceVectorField N_p = Sf_p/magSf_p;

// Dual area vector
const vectorList Sf = dual.areaVector();

// Dual face area
const scalarList magSf = mag(Sf);

// Dual material normals
const vectorList N = Sf/mag(Sf);

// Minimum edge length
const dimensionedScalar& h = op.minimumEdgeLength();

// Continuum mechanics class
mechanics mech(mesh, Sf);

// Boundary patches
//const polyBoundaryMesh& bm = mesh.boundaryMesh();

/*
       ______    ___      ___   _______     ______   _____  ____  _____       _     _________  ________   ______   
     .' ___  | .'   `.  .'   `.|_   __ \   |_   _ `.|_   _||_   \|_   _|     / \   |  _   _  ||_   __  |.' ____ \
    / .'   \_|/  .-.  \/  .-.  \ | |__) |    | | `. \ | |    |   \ | |      / _ \  |_/ | | \_|  | |_ \_|| (___ \_| 
    | |       | |   | || |   | | |  __ /     | |  | | | |    | |\ \| |     / ___ \     | |      |  _| _  _.____`.  
    \ `.___.'\\  `-'  /\  `-'  /_| |  \ \_  _| |_.' /_| |_  _| |_\   |_  _/ /   \ \_  _| |_    _| |__/ || \____) | 
     `.____ .' `.___.'  `.___.'|____| |___||______.'|_____||_____|\____||____| |____||_____|  |________| \______.' 
                                                                                                                      */
Info << "Creating mesh coordinate fields ..." << endl;

// Spatial nodal coordinate
pointVectorField x (
    IOobject("x", mesh),
    pMesh,
    dimensionedVector("x", dimensionSet(0,1,0,0,0,0,0), vector::zero)
);
x.primitiveFieldRef() = mesh.points();

// Material nodal coordinate
pointVectorField xw (
    IOobject("xw", mesh),
    pMesh,
    dimensionedVector("xw", dimensionSet(0,1,0,0,0,0,0), vector::zero)
);
xw.primitiveFieldRef() = mesh.points();

// Referential (previously Material) nodal coordinate
pointVectorField X( IOobject("X", mesh), x);

// Material edge centre coordinate
const vectorList& Xe = dual.edgeCentre();

// Spatial edge centre coordiantes
vectorList xe = Xe;

/*
    ______     ___   _____  ____   ____  ________  ______     ____   ____  _       _______     ______   
  .' ____ \  .'   `.|_   _||_  _| |_  _||_   __  ||_   _ `.  |_  _| |_  _|/ \     |_   __ \  .' ____ \
  | (___ \_|/  .-.  \ | |    \ \   / /    | |_ \_|  | | `. \   \ \   / / / _ \      | |__) | | (___ \_| 
   _.____`. | |   | | | |   _ \ \ / /     |  _| _   | |  | |    \ \ / / / ___ \     |  __ /   _.____`.  
  | \____) |\  `-'  /_| |__/ | \ ' /     _| |__/ | _| |_.' /     \ ' /_/ /   \ \_  _| |  \ \_| \____) | 
   \______.' `.___.'|________|  \_/     |________||______.'       \_/|____| |____||____| |___|\______.' 
                                                                                                           */
Info << "Creating strain measure fields ..." << endl;

// ALE - Spatial F (to be solved)
pointTensorField F ( IOobject("F", mesh), pMesh, Foam::tensor::I );

// ALE - Material J (to be solved)
pointScalarField matJ (
	IOobject ("matJ", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
	pMesh,
	dimensionedScalar("matJ", dimensionSet(0,0,0,0,0,0,0), 1)
);

// ALE - Material F (to be solved)
pointTensorField matF (
        IOobject("matF", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
	pMesh,
	dimensionedTensor("matF", F.dimensions(), Foam::tensor::I)
);

Info << "Creating linear momentum field ..." << endl;
// ALE - Referential LM (to be solved)
pointVectorField lm ( IOobject("lm", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), pMesh );

Info << "Creating material velocity field ..." << endl;
// ALE - Material velocity (to be solved)
pointVectorField solvedW(
		IOobject("solvedW", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), pMesh
		/*
        IOobject("solvedW",runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
        pMesh,
        dimensionedVector("solvedW", dimensionSet(0,1,-1,0,0,0,0), vector::zero)
	*/
);

Info << "Creating material energy field..." << endl;
pointScalarField E (
	IOobject("E", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), pMesh
		/*
        IOobject("E", runTime.timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),
	pMesh,
	mechanicalProperties.lookupOrDefault<dimensionedScalar>("E",0)
        //dimensionedScalar("E", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
	*/
);
/*
     ______     ___   _____     _____  ______     ____    ____   ___   ______   ________  _____     
   .' ____ \  .'   `.|_   _|   |_   _||_   _ `.  |_   \  /   _|.'   `.|_   _ `.|_   __  ||_   _|    
   | (___ \_|/  .-.  \ | |       | |    | | `. \   |   \/   | /  .-.  \ | | `. \ | |_ \_|  | |      
    _.____`. | |   | | | |   _   | |    | |  | |   | |\  /| | | |   | | | |  | | |  _| _   | |   _  
   | \____) |\  `-'  /_| |__/ | _| |_  _| |_.' /  _| |_\/_| |_\  `-'  /_| |_.' /_| |__/ | _| |__/ | 
    \______.' `.___.'|________||_____||______.'  |_____||_____|`.___.'|______.'|________||________| 
                                                                                                      */
Info << "Creating constitutive model ..." << endl;
// Deformation gradient tensor for P material model
pointTensorField trueF ( IOobject("trueF", mesh), pMesh, Foam::tensor::I );
// cofactor tensor for P material model
pointTensorField trueH ( IOobject("trueH", mesh), pMesh, Foam::tensor::I );
// Jacobian for P material model
pointScalarField trueJ ( IOobject("trueJ", mesh), pMesh, 1.0 );
// ALE - Initialization of trueF and Piola Hat // NOTE : trueJ and trueH are not necessary here
trueF = F & op.inverse(matF);
trueJ = det(trueF);
trueH = trueJ * op.invT(trueF);

// ALE - spatial J and H
pointScalarField spatJ = det(F);
pointTensorField spatH = spatJ * op.invT(F);

// Solid Model
solidModel model(trueF, mechanicalProperties, mesh);
model.correct(trueF, trueH, trueJ);

// Model quantities
const dimensionedScalar& rho   = model.density();
const dimensionedScalar& mu    = model.mu() ;
const dimensionedScalar& kappa = model.kappa() ;

// Pressure
pointScalarField p = model.pressure();
// First Piola Kirchhoff stress tensor
pointTensorField P = model.piola();

// ALE - CoFactor computed from solved quantities 
pointTensorField aleH = matJ * op.invT(matF);

// TO MOVE EVENTUALLY
forAll(E,n){
  E[n] = 0.5 * (lm[n]&lm[n])/rho.value()  +  0.5*mu.value()*(Foam::pow(trueJ[n],(-2.0/3.0))*(trueF[n]&&trueF[n])-3.0) + (0.5*kappa.value()*(trueJ[n]-1.0)*(trueJ[n]-1.0));
}


//vectorField n = mech.spatialNormal();

/*
  ____      ____  _  ____   ____  ________      ______   _______  ________  ________  ______     ______   
 |_  _|    |_  _|/ \|_  _| |_  _||_   __  |   .' ____ \ |_   __ \|_   __  ||_   __  ||_   _ `. .' ____ \
   \ \  /\  / / / _ \ \ \   / /    | |_ \_|   | (___ \_|  | |__) | | |_ \_|  | |_ \_|  | | `. \| (___ \_| 
    \ \/  \/ / / ___ \ \ \ / /     |  _| _     _.____`.   |  ___/  |  _| _   |  _| _   | |  | | _.____`.  
     \  /\  /_/ /   \ \_\ ' /     _| |__/ |   | \____) | _| |_    _| |__/ | _| |__/ | _| |_.' /| \____) | 
      \/  \/|____| |____|\_/     |________|    \______.'|_____|  |________||________||______.'  \______.' 
	                                                                                                      */
Info << "Creating fields for wave speeds ..." << endl;
// Stretch
pointScalarField stretch = mech.stretch();

// Wave speed for time increment
pointScalarField Up_time ( IOobject("Up_time", mesh), pMesh, model.Up() );
Up_time = Up_time/stretch;

// Longitudinal wave speed
pointScalarField Up
(
    IOobject("Up", mesh),
    pMesh,
    model.Up()/beta
);

// Shear wave speed
pointScalarField Us
(
    IOobject("Us", mesh),
    pMesh,
    model.Us()*beta
);

/*
        _       _____     ________    ____    ____   ___   ______   ________  _____     
       / \     |_   _|   |_   __  |  |_   \  /   _|.'   `.|_   _ `.|_   __  ||_   _|    
      / _ \      | |       | |_ \_|    |   \/   | /  .-.  \ | | `. \ | |_ \_|  | |      
     / ___ \     | |   _   |  _| _     | |\  /| | | |   | | | |  | | |  _| _   | |   _  
   _/ /   \ \_  _| |__/ | _| |__/ |   _| |_\/_| |_\  `-'  /_| |_.' /_| |__/ | _| |__/ | 
  |____| |____||________||________|  |_____||_____|`.___.'|______.'|________||________| 		                                                                                         
                                                                                           */

//- ALE: Declaration and initialization of ALE Objects
//pointVectorField& w_    = materialModel_.w();
//pointVectorField& wDot_ = materialModel_.wDot();

// Depending on tim mapping, we have to provide with good Init Cond
//solvedW.primitiveFieldRef() = w_.primitiveFieldRef();

// ALE - Material Model
Info << "Creating ALE Material Model ..." << nl;
aleModel materialModel_(mechanicalProperties, mesh, pMesh, "materialModel");
materialModel_.correct();
pointTensorField matP    = materialModel_.piola(matF, aleH, matJ);
//pointScalarField matPres = materialModel_.getMaterialPressure(matF, aleH, matJ);
//
scalarField lambdaH = sqrt((inter.pointToEdge(aleH) & N) & (inter.pointToEdge(aleH) & N));
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

// ALE - Spatial Model
Info << "Creating ALE Spatial Model ..." << nl;
aleModel spatialModel_(mechanicalProperties, mesh, pMesh, "spatialModel");
spatialModel_.correct();
pointTensorField spatP   = spatialModel_.piola(F,  spatH, spatJ);

// ALE - get Piola Star
/*
pointTensorField starP(
		       IOobject("starP", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), pMesh,
		       dimensionedTensor("starP", matP.dimensions(), tensor::zero)
		       );
starP   = ((1.0-alpha)*matP) + (alpha*(trueF.T()&spatP));
*/
pointTensorField starP = ((1.0-alpha)*matP) + (alpha*(trueF.T()&spatP));

/*
    ____   ____  ________  _____       ___      ______  _____  _________  _____  ________   ______   
   |_  _| |_  _||_   __  ||_   _|    .'   `.  .' ___  ||_   _||  _   _  ||_   _||_   __  |.' ____ \
     \ \   / /    | |_ \_|  | |     /  .-.  \/ .'   \_|  | |  |_/ | | \_|  | |    | |_ \_|| (___ \_| 
      \ \ / /     |  _| _   | |   _ | |   | || |         | |      | |      | |    |  _| _  _.____`.  
       \ ' /     _| |__/ | _| |__/ |\  `-'  /\ `.___.'\ _| |_    _| |_    _| |_  _| |__/ || \____) |
        \_/     |________||________| `.___.'  `.____ .'|_____|  |_____|  |_____||________| \______.' 	                                                                                                                                                                                                               */

// ALE: pR and pTilde, step quantities for rhsLM
pointVectorField pR
(
    IOobject ("pR", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    pMesh,
    dimensionedVector("pR", dimensionSet(1,-2,-1,0,0,0,0), vector::zero)
);
pR = op.inverseScalar(matJ) * lm;

pointVectorField pTilde
(
    IOobject ("pTilde", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    pMesh,
    dimensionedVector("pTilde", dimensionSet(1,-2,-1,0,0,0,0), vector::zero)
);
pTilde = matJ * pR;

//- ALE: v
pointVectorField v_ (
	IOobject("v_", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
	pMesh,
	dimensionedVector("v_", dimensionSet(0,1,-1,0,0,0,0), vector::zero)
);
v_ = pR/rho;

//- ALE:  vHat
pointVectorField vHat (
	IOobject("vHat", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
	pMesh,
	dimensionedVector("vHat", dimensionSet(0,1,-1,0,0,0,0), vector::zero)
);
vHat = v_ + (trueF & solvedW);

/*
    ____  ____       _       _____     ________  ____      ____  _     ____  ____  
   |_   ||   _|     / \     |_   _|   |_   __  ||_  _|    |_  _|/ \   |_  _||_  _| 
     | |__| |      / _ \      | |       | |_ \_|  \ \  /\  / / / _ \    \ \  / /   
     |  __  |     / ___ \     | |   _   |  _|      \ \/  \/ / / ___ \    \ \/ /    
    _| |  | |_  _/ /   \ \_  _| |__/ | _| |_        \  /\  /_/ /   \ \_  _|  |_    
   |____||____||____| |____||________||_____|        \/  \/|____| |____||______|   
                                                                                       */

//- ALE: Piola Hat
pointTensorField PHat (
    IOobject("PHat", runTime.timeName(), P.db(), IOobject::NO_READ, IOobject::NO_WRITE),
    P.mesh(),
    dimensionedTensor("PHat", dimensionSet(1,-1,-2,0,0,0,0), tensor::zero)
);
PHat = P & aleH ;

//- ALE: kinetic Energy
pointScalarField kinE (
    IOobject("kinE", runTime.timeName(), E.db(), IOobject::NO_READ, IOobject::NO_WRITE),
    E.mesh(),
    dimensionedScalar("kinE", E.dimensions(), 0.0)
);
//kinE =  0.5 * (lm&lm) / rho.value();

//- ALE: potential Energy
pointScalarField potE (
    IOobject("potE", runTime.timeName(), E.db(), IOobject::NO_READ, IOobject::NO_WRITE),
    E.mesh(),
    dimensionedScalar("potE", E.dimensions(), 0.0)
);
forAll(potE, n) {
    kinE[n] = 0.5 * (lm[n]&lm[n]) / rho.value();//rho.value() * (vHat[n]&vHat[n]) ;//(lm[n]&lm[n]) / rho.value();
    kinE[n] += 0.5 * rho.value() * (solvedW[n]&solvedW[n]);
    potE[n] = (0.5*mu.value()*Foam::pow(trueJ[n],-2./3.)*(trueF[n]&&trueF[n])) - (3*mu.value()/2.) + (0.5*kappa.value()*Foam::pow(trueJ[n]-1,2));
}

// - ALE: theoretical Energy
pointScalarField theoE (
    IOobject("theoE", runTime.timeName(), E.db(), IOobject::NO_READ, IOobject::NO_WRITE),
    E.mesh(),
    dimensionedScalar("theoE", E.dimensions(), 0.0)
);
theoE = kinE + potE;


/*
       ______  _______          _       ______   _____  ________  ____  _____  _________  
     .' ___  ||_   __ \        / \     |_   _ `.|_   _||_   __  ||_   \|_   _||  _   _  | 
    / .'   \_|  | |__) |      / _ \      | | `. \ | |    | |_ \_|  |   \ | |  |_/ | | \_| 
    | |   ____  |  __ /      / ___ \     | |  | | | |    |  _| _   | |\ \| |      | |     
    \ `.___]  |_| |  \ \_  _/ /   \ \_  _| |_.' /_| |_  _| |__/ | _| |_\   |_    _| |_    
     `._____.'|____| |___||____| |____||______.'|_____||________||_____|\____|  |_____|                                                                                         
                                                                                            */
Info << "Creating fields for gradient ..." << endl;

// Gradient class
gradientSchemes grad(mesh, Sf, V);

// Gradient of nodal linear momentum
tensorField vGrad = grad.gradient(v_);

// Reconstruction of nodal linear momentum
vectorField vM(mesh.edges().size());
vectorField vP(mesh.edges().size());
grad.reconstruct(v_, vGrad, vM, vP);

/*
      ________  _____  _____  _____  ____  ____  ________   ______   
     |_   __  ||_   _||_   _||_   _||_  _||_  _||_   __  |.' ____ \
       | |_ \_|  | |    | |    | |    \ \  / /    | |_ \_|| (___ \_| 
       |  _|     | |   _| '    ' |     > `' <     |  _| _  _.____`.  
      _| |_     _| |__/ |\ \__/ /    _/ /'`\ \_  _| |__/ || \____) | 
     |_____|   |________| `.__.'    |____||____||________| \______.' 
                                                                      
                                                                        */
Info << "Creating fields for fluxes ..." << endl;

// Edge Stabilisation
vectorField stabC = inter.pointToEdge(matJ*rho)*(0.5*mech.SmatrixMod(aleUp___, aleUs___) & (vP - vM));
// Edge traction
vectorField tC    = inter.pointToEdge(PHat) & N;  // case without Involution
//
pointVectorField Px ( IOobject("Px", mesh), pMesh, vector::zero );
pointVectorField Py ( IOobject("Py", mesh), pMesh, vector::zero );
pointVectorField Pz ( IOobject("Pz", mesh), pMesh, vector::zero );
forAll(Px, n) { 
	Px[n] = vector(P[n].xx(),P[n].yx(),P[n].zx()); 
	Py[n] = vector(P[n].xy(),P[n].yy(),P[n].zy()); 
	Pz[n] = vector(P[n].xz(),P[n].yz(),P[n].zz());
}
pointVectorField vector_100 ( IOobject("vector_100", mesh), pMesh, vector(1,0,0) );
pointVectorField vector_010 ( IOobject("vector_010", mesh), pMesh, vector(0,1,0) );
pointVectorField vector_001 ( IOobject("vector_001", mesh), pMesh, vector(0,0,1) );
//
tensorField P1C = inter.pointToEdge(Px) * N;   
tensorField P2C = inter.pointToEdge(Py) * N;
tensorField P3C = inter.pointToEdge(Pz) * N;     // case with Involution
// Edge linear momentum
tensorField vC = inter.pointToEdge(vHat) * N;
// Edge traction two
tensorField tC2 = inter.pointToEdge( pR ) * N;
// Edge material velocity
tensorField wC = inter.pointToEdge( solvedW )  * N;
// Edge alternative jacobian term
scalarField HTwC = inter.pointToEdge( aleH.T() & solvedW ) & N;
// ALE - The GCL technique distinction (commented code below) could
// be done, but for the sake of simplicity, initial values for both
// cases are declared anyway.
// It is also about not cluttering the code with IFs verifications.
/*
if (GCLtechnique == "withInvolution") {
    tC2 = inter.pointToEdge( pR )   * N;
} else if (GCLtechnique == "withoutInvolution") {
    HTwC = inter.pointToEdge( aleH.T()&solvedW ) & N;
} else {
    FatalErrorIn("updateVariables.H") << "GCL technique is not properly defined." << abort(FatalError);
}
*/
vectorField PmatHC    = inter.pointToEdge( matP & aleH ) & N;
vectorField pRwHC     = inter.pointToEdge( (pR*solvedW) & aleH ) & N;

scalarField Pmat_pmat = inter.pointToEdge( matP.T()&solvedW ) & N; // for energy, deprecated ?
vectorField PstarHC   = inter.pointToEdge( starP & aleH ) & N; // alternative Pstar term for PmatHc

/*
    _______     ____  ____   ______  
   |_   __ \   |_   ||   _|.' ____ \
     | |__) |    | |__| |  | (___ \_| 
     |  __ /     |  __  |   _.____`.  
    _| |  \ \_  _| |  | |_ | \____) | 
   |____| |___||____||____| \______.' 
                                         */
// ALE - RHS of jacobian equation
pointScalarField rhsJ (
	IOobject("rhsJ", mesh),
	pMesh,
	dimensionedScalar("rhsJ", dimensionSet(0,0,-1,0,0,0,0), 0)
);

// RHS of linear momentum equation
pointVectorField rhsLm (
    IOobject ("rhsLm", mesh),
    pMesh,
    dimensionedVector("rhsLm", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
);
pointVectorField rhsLm1 (
    IOobject ("rhsLm1", mesh),
    pMesh,
    dimensionedVector("rhsLm1", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
);
pointVectorField rhsLm2 (
    IOobject ("rhsLm2", mesh),
    pMesh,
    dimensionedVector("rhsLm2", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
);
pointVectorField rhsLm3 (
    IOobject ("rhsLm3", mesh),
    pMesh,
    dimensionedVector("rhsLm3", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
);
pointVectorField rhsLm4 (
    IOobject ("rhsLm4", mesh),
    pMesh,
    dimensionedVector("rhsLm4", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
);

rhsLm4.primitiveFieldRef() = op.surfaceSum(stabC*magSf) ; 
if (PHatTerm == "withoutInvolution") { 
    rhsLm1.primitiveFieldRef() = op.surfaceSum(tC*magSf);
} else if (PHatTerm == "withInvolution") {
    rhsLm1.primitiveFieldRef() = 
		    ((aleH && op.surfaceSum(P1C*magSf))*vector_100)
		  + ((aleH && op.surfaceSum(P2C*magSf))*vector_010)
                  + ((aleH && op.surfaceSum(P3C*magSf))*vector_001);  
}
if (GCLtechnique == "withInvolution") {
    rhsJ.primitiveFieldRef() = aleH && op.surfaceSum(wC*magSf);
    rhsLm2.primitiveFieldRef() = op.surfaceSum(tC2*magSf) & (aleH.T() & solvedW);
    rhsLm3.primitiveFieldRef() = pR * (aleH && (op.surfaceSum(wC*magSf)));
    //rhsLm.primitiveFieldRef() = rhsLm1 + rhsLm2  + rhsLm3;
} else if (GCLtechnique == "withoutInvolution") {
    rhsJ.primitiveFieldRef() = op.surfaceSum(HTwC*magSf);
    rhsLm2.primitiveFieldRef() = op.surfaceSum(pRwHC*magSf);
    //rhsLm.primitiveFieldRef() = rhsLm1 + rhsLm2;
} else {
    FatalErrorIn("createFields.H") << "GCL technique is not properly defined." << abort(FatalError);
}

//pointVectorField trac ( IOobject("trac",  runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), pMesh );

// RHS of F equation
pointTensorField rhsF (
    IOobject ("rhsF", mesh),
    pMesh,
    dimensionedTensor("rhsF",F.dimensions()/dimTime,tensor::zero)
);

// ALE - RHS of material deformation gradient equation
pointTensorField rhsMatF (
    IOobject ("rhsMatF", mesh),
    pMesh,
    dimensionedTensor("rhsMatF", matF.dimensions()/dimTime, tensor::zero)
);

// ALE - RHS of solvedW velocity
pointVectorField rhsW (
    IOobject ("rhsW",runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    pMesh,
    dimensionedVector("rhsW", dimensionSet(0,1,-2,0,0,0,0), vector::zero)
);

// ALE - RHS of E energy
pointScalarField rhsE1 ( IOobject ("rhsE1", mesh), pMesh,
    dimensionedScalar("rhsE1",E.dimensions()/dimTime, 0.0)
);
pointScalarField rhsE2 ( IOobject ("rhsE2", mesh), pMesh,
    dimensionedScalar("rhsE2",E.dimensions()/dimTime, 0.0)
);
pointScalarField rhsE ( IOobject ("rhsE", mesh), pMesh,
    dimensionedScalar("rhsE",E.dimensions()/dimTime, 0.0)
);
rhsE1.primitiveFieldRef() = op.surfaceSum( inter.pointToEdge(PHat.T()&(lm/rho.value()))             & Sf );
rhsE2.primitiveFieldRef() = op.surfaceSum( inter.pointToEdge((aleH.T()&solvedW) * E)  & Sf );
rhsE = rhsE1 + rhsE2;

/*
   ______   _____   ______   _______  _____          _        ______  ________  ____    ____  ________  ____  _____  _________   ______   
   |_   _ `.|_   _|.' ____ \ |_   __ \|_   _|        / \     .' ___  ||_   __  ||_   \  /   _||_   __  ||_   \|_   _||  _   _  |.' ____ \
     | | `. \ | |  | (___ \_|  | |__) | | |         / _ \   / .'   \_|  | |_ \_|  |   \/   |    | |_ \_|  |   \ | |  |_/ | | \_|| (___ \_| 
     | |  | | | |   _.____`.   |  ___/  | |   _    / ___ \  | |         |  _| _   | |\  /| |    |  _| _   | |\ \| |      | |     _.____`.  
    _| |_.' /_| |_ | \____) | _| |_    _| |__/ | _/ /   \ \_\ `.___.'\ _| |__/ | _| |_\/_| |_  _| |__/ | _| |_\   |_    _| |_   | \____) | 
   |______.'|_____| \______.'|_____|  |________||____| |____|`.____ .'|________||_____||_____||________||_____|\____|  |_____|   \______.' 
                                                                                                                                              */

Info << "Creating fields for post-processing ..." << endl;

// Nodal displacements
pointVectorField u (
    IOobject ("u", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    pMesh,
    dimensionedVector("u", dimLength, vector::zero)
);
// w Nodal displacements
pointVectorField uw (
    IOobject ("uw", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    pMesh,
    dimensionedVector("uw", dimLength, vector::zero)
);
X -= u;

/*
  _________  _____  ____    ____  ________   ____   ____  _       _______     ______   
 |  _   _  ||_   _||_   \  /   _||_   __  | |_  _| |_  _|/ \     |_   __ \  .' ____ \
 |_/ | | \_|  | |    |   \/   |    | |_ \_|   \ \   / / / _ \      | |__) | | (___ \_| 
     | |      | |    | |\  /| |    |  _| _     \ \ / / / ___ \     |  __ /   _.____`.  
    _| |_    _| |_  _| |_\/_| |_  _| |__/ |     \ ' /_/ /   \ \_  _| |  \ \_| \____) | 
   |_____|  |_____||_____||_____||________|      \_/|____| |____||____| |___|\______.' 
	                                                                                  */
Info << "Creating variables for time ..." << endl;

scalar tstep = 0;
dimensionedScalar deltaT("deltaT", dimTime, runTime.deltaTValue());
dimensionedScalar t("t", dimTime, 0.0);
if (timeStepping == "variable") {
    deltaT = (cfl*h)/model.Up();
    runTime.setDeltaT(deltaT);
}


scalarList RKstage(2);
RKstage[0] = 0;
RKstage[1] = 1;


/*
      _______  _______     _____  ____  _____  _________  
     |_   __ \|_   __ \   |_   _||_   \|_   _||  _   _  | 
       | |__) | | |__) |    | |    |   \ | |  |_/ | | \_| 
       |  ___/  |  __ /     | |    | |\ \| |      | |     
      _| |_    _| |  \ \_  _| |_  _| |_\   |_    _| |_    
     |_____|  |____| |___||_____||_____|\____|  |_____|   
                                                            */

Info << "Printing data ..." << endl;
model.printMaterialProperties();
materialModel_.printProperties();
spatialModel_.printProperties();

//Info << "Printing global momentum ..." << endl;
Info << "Printing mesh centroids" << nl;
dual.printPrimalMeshCentroid();
dual.printDualMeshCentroid();

Info << "\nObjects all created !" << nl << nl << nl;
