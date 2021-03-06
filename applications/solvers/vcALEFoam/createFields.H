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

// Initialize paper Functions
paperFunctions pf(mechanicalProperties, paperFunctionsProperties, mesh, pMesh);

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
    IOobject ("u", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    pMesh,
    dimensionedVector("u", dimLength, vector::zero)
);
// w Nodal displacements
pointVectorField uw (
    IOobject ("uw", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    pMesh,
    dimensionedVector("uw", dimLength, vector::zero)
);
X -= u;

const word& tutorial(runParameters.lookupOrDefault<word>("tutorial", "None"));
if (tutorial == "swingingCube") {
    const scalar& A = readScalar(runParameters.lookup("A"));
    const scalar& B = readScalar(runParameters.lookup("B"));
    const scalar& C = readScalar(runParameters.lookup("C"));
    const scalar pi = Foam::constant::mathematical::pi;
    pointScalarField U0 (IOobject("U0", mesh), pMesh, dimensionedScalar(runParameters.lookup("U0")));

    forAll(mesh.points(), n) {
        const scalar X = mesh.points()[n][0];
        const scalar Y = mesh.points()[n][1];
        const scalar Z = mesh.points()[n][2];
        u[n] = U0[n] * vector(
             A * Foam::sin(pi*X/2.) * Foam::cos(pi*Y/2.) * Foam::cos(pi*Z/2.),
             B * Foam::cos(pi*X/2.) * Foam::sin(pi*Y/2.) * Foam::cos(pi*Z/2.),
             C * Foam::cos(pi*X/2.) * Foam::cos(pi*Y/2.) * Foam::sin(pi*Z/2.)
                                     );
   }
   x = u + X; 
}
/*
    ______     ___   _____  ____   ____  ________  ______     ____   ____  _       _______     ______   
  .' ____ \  .'   `.|_   _||_  _| |_  _||_   __  ||_   _ `.  |_  _| |_  _|/ \     |_   __ \  .' ____ \
  | (___ \_|/  .-.  \ | |    \ \   / /    | |_ \_|  | | `. \   \ \   / / / _ \      | |__) | | (___ \_| 
   _.____`. | |   | | | |   _ \ \ / /     |  _| _   | |  | |    \ \ / / / ___ \     |  __ /   _.____`.  
  | \____) |\  `-'  /_| |__/ | \ ' /     _| |__/ | _| |_.' /     \ ' /_/ /   \ \_  _| |  \ \_| \____) | 
   \______.' `.___.'|________|  \_/     |________||______.'       \_/|____| |____||____| |___|\______.' 
                                                                                                           */
Info << "Creating strain measure fields ..." << endl;

// ALE - SPatial J (to be solved)
pointScalarField spatJ (
        IOobject ("spatJ", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), pMesh,
        dimensionedScalar("spatJ", dimensionSet(0,0,0,0,0,0,0), 1)	
);

// ALE - Spatial F (to be solved)
pointTensorField F ( IOobject("F", mesh), pMesh, Foam::tensor::I );

// ALE - Material J (to be solved)
pointScalarField matJ (
	IOobject ("matJ", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), pMesh,
	dimensionedScalar("matJ", dimensionSet(0,0,0,0,0,0,0), 1)
);

// ALE - Material F (to be solved)
pointTensorField matF (
        IOobject("matF", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), pMesh,
	dimensionedTensor("matF", F.dimensions(), Foam::tensor::I)
);

Info << "Creating linear momentum field ..." << endl;
// ALE - Referential LM (to be solved)
pointVectorField lm ( IOobject("lm", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), pMesh );

Info << "Creating material velocity field ..." << endl;
// ALE - Material velocity (to be solved)
pointVectorField solvedW(
	IOobject("solvedW", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE), pMesh,
	dimensionedVector("solvedW", dimLength/dimTime, vector::zero)
);

Info << "Creating energy field ..." << endl;
pointScalarField E (
	IOobject("E", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE), pMesh,
	dimensionedScalar("E", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
);

Info << "Creating Plasticity solved variables ..." << endl;
pointTensorField DH (
	IOobject("DH", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), pMesh,
	dimensionedTensor("DH", dimless, tensor::I)
);
// plastic strain epsilon
pointScalarField alphaH ( 
	IOobject("alphaH", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), 
	pMesh,
	dimensionedScalar("alphaH", dimless, 0.0)
);


/*
     ______     ___   _____     _____  ______     ____    ____   ___   ______   ________  _____     
   .' ____ \  .'   `.|_   _|   |_   _||_   _ `.  |_   \  /   _|.'   `.|_   _ `.|_   __  ||_   _|    
   | (___ \_|/  .-.  \ | |       | |    | | `. \   |   \/   | /  .-.  \ | | `. \ | |_ \_|  | |      
    _.____`. | |   | | | |   _   | |    | |  | |   | |\  /| | | |   | | | |  | | |  _| _   | |   _  
   | \____) |\  `-'  /_| |__/ | _| |_  _| |_.' /  _| |_\/_| |_\  `-'  /_| |_.' /_| |__/ | _| |__/ | 
    \______.' `.___.'|________||_____||______.'  |_____||_____|`.___.'|______.'|________||________| 
                                                                                                      */
Info << "Creating Plasticity scaled variables ..." << nl;
pointTensorField D = DH / matJ;
pointScalarField alpha ( 
	IOobject("alpha", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), alphaH / matJ
		);

Info << "Creating constitutive model ..." << endl;
// Deformation gradient tensor for P material model
pointTensorField trueF ( IOobject("trueF", mesh), pMesh, Foam::tensor::I );
// ALE - Initialization of trueF and Piola Hat // NOTE : trueJ and trueH are not necessary here
trueF = F & op.inverse(matF);



if (tutorial == "swingingCube") { 
    const scalar& A = readScalar(runParameters.lookup("A"));
    const scalar& B = readScalar(runParameters.lookup("B"));
    const scalar& C = readScalar(runParameters.lookup("C"));
    const scalar pi = Foam::constant::mathematical::pi;
    pointScalarField U0 (IOobject("U0", mesh), pMesh, dimensionedScalar(runParameters.lookup("U0")));
    forAll(mesh.points(), n) {// += because F is already I
      const scalar X = mesh.points()[n][0];
      const scalar Y = mesh.points()[n][1];
      const scalar Z = mesh.points()[n][2];
      trueF[n] += pi/2. * U0[n] * tensor(
	            A * Foam::cos(pi*X/2.) * Foam::cos(pi*Y/2.) * Foam::cos(pi*Z/2.),
		    -A * Foam::sin(pi*X/2.) * Foam::sin(pi*Y/2.) * Foam::cos(pi*Z/2.),
		    -A * Foam::sin(pi*X/2.) * Foam::cos(pi*Y/2.) * Foam::sin(pi*Z/2.),
		    //
		    -B * Foam::sin(pi*X/2.) * Foam::sin(pi*Y/2.) * Foam::cos(pi*Z/2.),
		    B * Foam::cos(pi*X/2.) * Foam::cos(pi*Y/2.) * Foam::cos(pi*Z/2.),
		    -B * Foam::cos(pi*X/2.) * Foam::sin(pi*Y/2.) * Foam::sin(pi*Z/2.),
		    //
		    -C * Foam::sin(pi*X/2.) * Foam::cos(pi*Y/2.) * Foam::sin(pi*Z/2.),
		    -C * Foam::cos(pi*X/2.) * Foam::sin(pi*Y/2.) * Foam::sin(pi*Z/2.),
		    C * Foam::cos(pi*X/2.) * Foam::cos(pi*Y/2.) * Foam::cos(pi*Z/2.)
		      );
      F = trueF;
    }
} 

// cofactor tensor for P material model
pointTensorField trueH ( IOobject("trueH", mesh), pMesh, Foam::tensor::I );
// Jacobian for P material model
pointScalarField trueJ ( IOobject("trueJ", mesh), pMesh, 1.0 );
trueJ = det(F);//trueJ = spatJ; trueJ *= op.inverseScalar(matJ);//trueJ = /* det(trueF) */ spatJ * op.inverseScalar(matJ);
trueH = trueJ * op.invT(trueF);

// ALE - spatial J and H
//pointScalarField spatJ = det(F);
//
pointScalarField spatJdet(
        IOobject ("spatJdet", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        pMesh,
        dimensionedScalar("spatJdet", dimensionSet(0,0,0,0,0,0,0), 1)
);
spatJdet = det(F);
pointTensorField spatH = spatJdet * op.invT(F);//spatJ * op.invT(F);

// Solid Model
solidModel model(trueF, mechanicalProperties, mesh);
const word& solidModel__ = mechanicalProperties.subDict("solidModel").lookup("solidModel");
if (solidModel__ == "vonMises") {
    D = DH / matJ; alpha = alphaH / matJ;
    model.correct(trueF, trueH, trueJ, spatJ, matJ, D, alpha);
    forAll(D, n) { DH[n] = D[n] * matJ[n];  alphaH[n] = alpha[n] * matJ[n]; }
} else {
    model.correct(trueF, trueH, trueJ, spatJ, matJ);
}

// Model quantities
const dimensionedScalar& rho   = model.density();
const dimensionedScalar& mu    = model.mu() ;
const dimensionedScalar& kappa = model.kappa() ;

// Pressure
pointScalarField p = model.pressure(); p.write();
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
pointScalarField Up (
    IOobject("Up", mesh), pMesh, model.Up()/beta
);

// Shear wave speed
pointScalarField Us (
    IOobject("Us", mesh), pMesh, model.Us()*beta
);

/*
        _       _____     ________    ____    ____   ___   ______   ________  _____     
       / \     |_   _|   |_   __  |  |_   \  /   _|.'   `.|_   _ `.|_   __  ||_   _|    
      / _ \      | |       | |_ \_|    |   \/   | /  .-.  \ | | `. \ | |_ \_|  | |      
     / ___ \     | |   _   |  _| _     | |\  /| | | |   | | | |  | | |  _| _   | |   _  
   _/ /   \ \_  _| |__/ | _| |__/ |   _| |_\/_| |_\  `-'  /_| |_.' /_| |__/ | _| |__/ | 
  |____| |____||________||________|  |_____||_____|`.___.'|______.'|________||________| 		                                                                                         
                                                                                           */
// ALE - Material Model
Info << "Creating ALE Material Model ..." << nl;
aleModel materialModel_(mechanicalProperties, mesh, pMesh, "materialModel");
materialModel_.correct();
pointTensorField matP    = materialModel_.piola(matF, aleH, matJ);
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
// spatJdet and spatH are declared in 'solid model' up here
pointTensorField spatP   = spatialModel_.piola(F,  spatH, spatJdet);

// ALE - get Piola Star
pointTensorField starP = (0.5*muMaterial_*matP) + (0.5*muSpatial_*(pf.Fbar(trueF).T()&spatP));
//((1.0-alpha_)*matP) + (alpha_*(trueF.T()&spatP));

/*
    ____   ____  ________  _____       ___      ______  _____  _________  _____  ________   ______   
   |_  _| |_  _||_   __  ||_   _|    .'   `.  .' ___  ||_   _||  _   _  ||_   _||_   __  |.' ____ \
     \ \   / /    | |_ \_|  | |     /  .-.  \/ .'   \_|  | |  |_/ | | \_|  | |    | |_ \_|| (___ \_| 
      \ \ / /     |  _| _   | |   _ | |   | || |         | |      | |      | |    |  _| _  _.____`.  
       \ ' /     _| |__/ | _| |__/ |\  `-'  /\ `.___.'\ _| |_    _| |_    _| |_  _| |__/ || \____) |
        \_/     |________||________| `.___.'  `.____ .'|_____|  |_____|  |_____||________| \______.' 	                                                                                                                                                                                                                                                                     */
//- ALE: pR
pointVectorField pR ( IOobject ("pR", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), 
    op.inverseScalar(matJ) * lm 
);

//- ALE: v
pointVectorField v_ ( IOobject("v_", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    pR/rho 
);

//- ALE:  vHat
pointVectorField vHat ( IOobject("vHat", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), 
    v_ + (trueF & solvedW) 
);

/*
    ____  ____       _       _____     ________  ____      ____  _     ____  ____  
   |_   ||   _|     / \     |_   _|   |_   __  ||_  _|    |_  _|/ \   |_  _||_  _| 
     | |__| |      / _ \      | |       | |_ \_|  \ \  /\  / / / _ \    \ \  / /   
     |  __  |     / ___ \     | |   _   |  _|      \ \/  \/ / / ___ \    \ \/ /    
    _| |  | |_  _/ /   \ \_  _| |__/ | _| |_        \  /\  /_/ /   \ \_  _|  |_    
   |____||____||____| |____||________||_____|        \/  \/|____| |____||______|   
                                                                                       */

//- ALE: Piola Hat
pointTensorField PHat ( IOobject("PHat", runTime.timeName(), P.db(), IOobject::NO_READ, IOobject::NO_WRITE),
    P & aleH
);

//- ALE: kinetic Energy
pointScalarField kinE ( IOobject("kinE", runTime.timeName(), E.db(), IOobject::NO_READ, IOobject::NO_WRITE), //E.mesh(),
    0.5 * (lm&lm) /rho //dimensionedScalar("kinE", E.dimensions(), 0.0)
);

//- ALE: potential Energy
pointScalarField potE ( IOobject("potE", runTime.timeName(), E.db(), IOobject::NO_READ, IOobject::NO_WRITE), E.mesh(),
    //(0.5*mu*Foam::pow(trueJ.primitiveField(), -2./3.) * (trueF&&trueF)) //- (3.*mu/2.) + (0.5*kappa*Foam::pow(trueJ-1,2.)) 
    dimensionedScalar("potE", E.dimensions(), 0.0)
);
//potE.internalField() = (0.5*mu*Foam::pow(trueJ.internalField(), -2./3.) );//* (trueF&&trueF)) ;//- (3.*mu/2.) + (0.5*kappa*Foam::pow(trueJ-1,2.))

forAll(potE, n) {
    //kinE[n] = 0.5 * (lm[n]&lm[n]) / rho.value();//rho.value() * (vHat[n]&vHat[n]) ;//(lm[n]&lm[n]) / rho.value();
    //kinE[n] += 0.5 * rho.value() * (solvedW[n]&solvedW[n]);
    potE[n] = (0.5*mu.value()*Foam::pow(trueJ[n],-2./3.)*(trueF[n]&&trueF[n])) - (3*mu.value()/2.) + (0.5*kappa.value()*Foam::pow(trueJ[n]-1,2));
}


// - ALE: theoretical Energy
pointScalarField theoE (
    IOobject("theoE", runTime.timeName(), E.db(), IOobject::NO_READ, IOobject::NO_WRITE), //E.mesh(),
    kinE + potE //dimensionedScalar("theoE", E.dimensions(), 0.0)
);
//theoE = kinE + potE;


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

// Gradient and reconstruction of spatial velocity
tensorField vGrad = grad.gradient(v_);
vectorField vM(mesh.edges().size());
vectorField vP(mesh.edges().size());
grad.reconstruct(v_, vGrad, vM, vP);

// Gradient and reconstruction of pressure
vectorField pGrad = grad.gradient(p);
scalarField pM(mesh.edges().size());
scalarField pP(mesh.edges().size());
grad.reconstruct(p, pGrad, pM, pP);

// Gradient and reconstruction of nodal material velocity
tensorField wGrad = grad.gradient(solvedW);
vectorField wM(mesh.edges().size());
vectorField wP(mesh.edges().size());
grad.reconstruct(solvedW, wGrad, wM, wP);

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
vectorField stabC    = inter.pointToEdge(matJ*rho)*(0.5*mech.SmatrixMod(aleUp___, aleUs___) & (vP - vM));
tensorField stabJacC =   ((0.5*mech.SmatrixJacobianMod(aleUp___, aleUs___)) * (pP - pM))/rho.value(); 
vectorField stabWC   = rho.value() * 0.5 * (mech.SmatrixMod(aleUp___, aleUs___) & (wP - wP));

// Edge traction
vectorField tC    = inter.pointToEdge(PHat) & N;

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
vector v100(1,0,0), v010(0,1,0), v001(0,0,1); // used in gEqns.H
//
tensorField P1C = inter.pointToEdge(Px) * N;   
tensorField P2C = inter.pointToEdge(Py) * N;
tensorField P3C = inter.pointToEdge(Pz) * N;

// Edge linear momentum
//tensorField vC = inter.pointToEdge(vHat) * N;   // not used anymore
// Edge traction two
tensorField tC2 = inter.pointToEdge( pR ) * N;
// Edge material velocity
//tensorField wC = inter.pointToEdge( solvedW )  * N;   // not used anymore
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
//vectorField PmatHC    = inter.pointToEdge( matP & aleH ) & N;          // not used anymore
vectorField pRwHC     = inter.pointToEdge( (pR*solvedW) & aleH ) & N;

vectorField PstarHC   = inter.pointToEdge( starP & aleH ) & N;              // not used anymore

pointVectorField D1 ( IOobject("D1", mesh), pMesh, vector::zero );
pointVectorField D2 ( IOobject("D2", mesh), pMesh, vector::zero );
pointVectorField D3 ( IOobject("D3", mesh), pMesh, vector::zero );
forAll(D1, n) { 
	D1[n] = vector(D[n].xx(), D[n].xy(), D[n].xz());
	D2[n] = vector(D[n].yx(), D[n].yy(), D[n].yz());
	D3[n] = vector(D[n].zx(), D[n].zy(), D[n].zz());
}

/*
    _______     ____  ____   ______  
   |_   __ \   |_   ||   _|.' ____ \
     | |__) |    | |__| |  | (___ \_| 
     |  __ /     |  __  |   _.____`.  
    _| |  \ \_  _| |  | |_ | \____) | 
   |____| |___||____||____| \______.' 
                                         */
// ALE - RHS of jacobian equation
pointScalarField rhsJ ( IOobject("rhsJ", mesh), pMesh,
	dimensionedScalar("rhsJ", dimensionSet(0,0,-1,0,0,0,0), 0)
);

// RHS of linear momentum equation
pointVectorField rhsLm ( IOobject ("rhsLm", mesh), pMesh,
    dimensionedVector("rhsLm", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
); //useless

pointVectorField rhsLm1 ( IOobject ("rhsLm1", mesh), pMesh,
    dimensionedVector("rhsLm1", lm.dimensions()/dimTime, vector::zero)
);
pointVectorField rhsLm2 ( IOobject ("rhsLm2", mesh), pMesh,
    dimensionedVector("rhsLm2", lm.dimensions()/dimTime, vector::zero)
);
pointVectorField rhsLm3 ( IOobject ("rhsLm3", mesh), pMesh,
    dimensionedVector("rhsLm3", lm.dimensions()/dimTime, vector::zero)
);
pointVectorField rhsLm4 ( IOobject ("rhsLm4", mesh), pMesh,
    dimensionedVector("rhsLm4", lm.dimensions()/dimTime, vector::zero)
);
rhsLm4.primitiveFieldRef() = op.surfaceSum(stabC*magSf) ; 

if (PHatTerm == "withInvolution") {
    rhsLm1.primitiveFieldRef() = 
	((aleH && op.surfaceSum(P1C*magSf))*vector_100) + ((aleH && op.surfaceSum(P2C*magSf))*vector_010) + ((aleH && op.surfaceSum(P3C*magSf))*vector_001);  
} else if (PHatTerm == "withoutInvolution") { 
    rhsLm1.primitiveFieldRef() = op.surfaceSum(tC*magSf);
}

if (GCLtechnique == "withInvolution") {
    rhsJ.primitiveFieldRef() = aleH && op.surfaceSum( inter.pointToEdge(solvedW) * Sf);
    rhsLm2.primitiveFieldRef() = op.surfaceSum(tC2*magSf) & (aleH.T() & solvedW);
    rhsLm3.primitiveFieldRef() = pR * (aleH && (op.surfaceSum( inter.pointToEdge(solvedW) * Sf)));
    //rhsLm.primitiveFieldRef() = rhsLm1 + rhsLm2  + rhsLm3;
} else if (GCLtechnique == "withoutInvolution") {
    rhsJ.primitiveFieldRef() = op.surfaceSum(HTwC*magSf);
    rhsLm2.primitiveFieldRef() = op.surfaceSum(pRwHC*magSf);
    //rhsLm.primitiveFieldRef() = rhsLm1 + rhsLm2;
} else {
    FatalErrorIn("createFields.H") << "GCL technique is not properly defined." << abort(FatalError);
}

// RHS of spatial J
pointScalarField rhsSpatJ1 ( IOobject("rhsSpatJ1", mesh), pMesh,
    dimensionedScalar("rhsSpatJ1", spatJ.dimensions()/dimTime, 0)
);
rhsSpatJ1.primitiveFieldRef() = aleH && op.surfaceSum( inter.pointToEdge( v_ ) * Sf );	

pointScalarField rhsSpatJ2 ( IOobject("rhsSpatJ2", mesh), pMesh,
    dimensionedScalar("rhsSpatJ2", spatJ.dimensions()/dimTime, 0)
);
rhsSpatJ2.primitiveFieldRef() = op.surfaceSum( inter.pointToEdge( (spatH.T()&trueF)&solvedW ) & Sf );

pointScalarField rhsSpatJ3 ( IOobject("rhsSpatJ3", mesh), pMesh,
    dimensionedScalar("rhsSpatJ3", spatJ.dimensions()/dimTime, 0)
);
rhsSpatJ3.primitiveFieldRef() = op.surfaceSum( (stabJacC & (inter.pointToEdge(spatH)&N) ) & (inter.pointToEdge(spatH) & Sf) ) ;

pointScalarField rhsSpatJalt ( IOobject("rhsSpatJalt", mesh), pMesh,
    dimensionedScalar("rhsSpatJalt", spatJ.dimensions()/dimTime, 0)
);
rhsSpatJalt.primitiveFieldRef() = op.surfaceSum( inter.pointToEdge(spatH&vHat) & Sf );

// RHS of spatial F equation
pointTensorField rhsF ( IOobject ("rhsF", mesh), pMesh,
    dimensionedTensor("rhsF", F.dimensions()/dimTime,tensor::zero)
);
rhsF.primitiveFieldRef()    = op.surfaceSum( inter.pointToEdge(vHat)    * Sf );

// ALE - RHS of material deformation gradient equation
pointTensorField rhsMatF ( IOobject ("rhsMatF", mesh), pMesh,
    dimensionedTensor("rhsMatF", matF.dimensions()/dimTime, tensor::zero)
);
rhsMatF.primitiveFieldRef() = op.surfaceSum( inter.pointToEdge(solvedW) * Sf );

// ALE - RHS of solvedW velocity
pointVectorField rhsW ( IOobject ("rhsW", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), pMesh,
    dimensionedVector("rhsW", solvedW.dimensions()/dimTime, vector::zero)
);
rhsW.primitiveFieldRef() = op.surfaceSum( inter.pointToEdge(starP &aleH ) & Sf ) / rho.value();

pointVectorField rhsW2 ( IOobject ("rhsW", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), pMesh,
    dimensionedVector("rhsW2", solvedW.dimensions()/dimTime, vector::zero)
);
rhsW2.primitiveFieldRef() = op.surfaceSum( stabWC*magSf );

// ALE - RHS of E energy
pointScalarField rhsE1 ( IOobject ("rhsE1", mesh), pMesh,
    dimensionedScalar("rhsE1", E.dimensions()/dimTime, 0.0)
);
pointScalarField rhsE2 ( IOobject ("rhsE2", mesh), pMesh,
    dimensionedScalar("rhsE2", E.dimensions()/dimTime, 0.0)
);
pointScalarField rhsE ( IOobject ("rhsE", mesh), pMesh,
    dimensionedScalar("rhsE", E.dimensions()/dimTime, 0.0)
);
rhsE1.primitiveFieldRef() = op.surfaceSum( inter.pointToEdge(PHat.T()&(lm/rho.value()))             & Sf );
rhsE2.primitiveFieldRef() = op.surfaceSum( inter.pointToEdge((aleH.T()&solvedW) * E/matJ)  & Sf );
rhsE = rhsE1 + rhsE2;

// ALE - RHS of plasticity variables
pointTensorField rhsD ( IOobject("rhsD", mesh), pMesh,
    dimensionedTensor("rhsD", DH.dimensions()/dimTime, tensor::zero)
);
rhsD.primitiveFieldRef()     = //op.surfaceSum( inter.pointToEdge((D&aleH.T()) & solvedW) * Sf );
             ( op.surfaceSum( inter.pointToEdge(D1*(aleH.T()&solvedW)) & Sf ) * vector_100)
	   + ( op.surfaceSum( inter.pointToEdge(D2*(aleH.T()&solvedW)) & Sf ) * vector_010)
	   + ( op.surfaceSum( inter.pointToEdge(D3*(aleH.T()&solvedW)) & Sf ) * vector_001);

pointScalarField rhsAlpha ( IOobject("rhsAlpha", mesh), pMesh,
    dimensionedScalar("rhsAlpha", alphaH.dimensions()/dimTime, 0.0)
);
rhsAlpha.primitiveFieldRef() = op.surfaceSum( inter.pointToEdge(alpha*(aleH.T()&solvedW)) & Sf );


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

// Alternative implementation of RK coefficient table
// the first coef in the table is for U^n 
// the second coef is for the solved term
//int nbRKstages = 3; // 2 or 3;
double RKcoefTable[nbRKstages][2];
//int RKcoefTable[nbRKstages][2];
if (nbRKstages == 2) {
    //double RKcoefTable[nbRKstages][2] = {{0,1}, {.5,.5}};
    RKcoefTable[0][0] = 0;
    RKcoefTable[0][1] = 1;
    RKcoefTable[1][0] = .5;
    RKcoefTable[1][1] = .5;
} else if (nbRKstages == 3) {
    //double RKcoefTable[nbRKstages][2] = {{0,1}, {.75,.25}, {1./3.,2./3.}};
    RKcoefTable[0][0] = 0;
    RKcoefTable[0][1] = 1;
    RKcoefTable[1][0] = .75;
    RKcoefTable[1][1] = .25;
    RKcoefTable[2][0] = 1./3.;
    RKcoefTable[2][1] = 2./3.;    
} else {
   FatalErrorIn("createFields.H") << "Wrong number of RK steps" << abort(FatalError);
}
Info << "Integration method: " << nbRKstages << " stages Runge Kutta." << nl;

/*
      _______  _______     _____  ____  _____  _________  
     |_   __ \|_   __ \   |_   _||_   \|_   _||  _   _  | 
       | |__) | | |__) |    | |    |   \ | |  |_/ | | \_| 
       |  ___/  |  __ /     | |    | |\ \| |      | |     
      _| |_    _| |  \ \_  _| |_  _| |_\   |_    _| |_    
     |_____|  |____| |___||_____||_____|\____|  |_____|   
                                                            */

Info << "Printing data ..." << endl;
model.printProperties();
materialModel_.printProperties();
spatialModel_.printProperties();

//Info << "Printing global momentum ..." << endl;
Info << "Printing mesh centroids" << nl;
dual.printPrimalMeshCentroid();
dual.printDualMeshCentroid();

Info << "\nObjects all created !" << nl << nl << nl;


/*
pointVectorField Forces (IOobject("Forces", mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE), pMesh, dimensionedVector("Forces", dimensionSet(1,1,-2,0,0,0,0), vector::zero));
*/
