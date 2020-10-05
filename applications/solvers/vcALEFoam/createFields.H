Info << "\nCreating mesh parameters ..." << endl;

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

// Boundary patches
const polyBoundaryMesh& bm = mesh.boundaryMesh();
const label& freePatchID = bm.findPatchID("free");
const label& tractionPatchID = bm.findPatchID("traction");


Info << "Creating mesh coordinate fields ..." << endl;

// Spatial nodal coordinate
pointVectorField x
(
    IOobject("x", mesh),
    pMesh,
    dimensionedVector("x", dimensionSet(0,1,0,0,0,0,0), vector::zero)
);
x.primitiveFieldRef() = mesh.points();

// Material nodal coordinate
pointVectorField xw
(
    IOobject("xw", mesh),
    pMesh,
    dimensionedVector("xw", dimensionSet(0,1,0,0,0,0,0), vector::zero)
);
xw.primitiveFieldRef() = mesh.points();

// Referential (previously Material) nodal coordinate
pointVectorField X
(
    IOobject ("X", mesh),
    x
);

// Material edge centre coordinate
const vectorList& Xe = dual.edgeCentre();

// Spatial edge centre coordiantes
vectorList xe = Xe;


Info << "Creating linear momentum fields ..." << endl;

// Nodal linear momentum
pointVectorField lm
(
    IOobject ("lm", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    pMesh
);


Info << "Creating strain measure fields ..." << endl;

// Deformation gradient tensor
pointTensorField F
(
    IOobject("F", mesh),
    pMesh,
    Foam::tensor::I
);

// RHS of F equation
pointTensorField rhsF
(
    IOobject ("rhsF", mesh),
    pMesh,
    dimensionedTensor("rhsF",F.dimensions()/dimTime,tensor::zero)
);

// Cofactor of deformation
pointTensorField H
(
    IOobject ("H", mesh),
    det(F)*op.invT(F)
);

// Jacobian of deformation
pointScalarField J
(
    IOobject("J", mesh),
    det(F)
);

//Deformation gradient tensor for P material model
pointTensorField trueF
(
 	IOobject("trueF", mesh),
	pMesh,
	Foam::tensor::I
);

//cofactor tensor for P material model
pointTensorField trueH
(
 	IOobject("trueH", mesh),
	pMesh,
	Foam::tensor::I
);

//Jacobian for P material model
pointScalarField trueJ
(
 	IOobject("trueJ", mesh),
	pMesh,
	1.0
);



Info << "Creating constitutive model ..." << endl;

// Solid model class
solidModel model(trueF, mechanicalProperties, mesh);
model.correct(trueF, trueH, trueJ);

// Density
const dimensionedScalar& rho = model.density();

// Pressure
pointScalarField p = model.pressure();
p.write();


// First Piola Kirchhoff stress tensor
pointTensorField P = model.piola();

//- ALE: Piola Hat
pointTensorField PHat
(
    IOobject
    (
        "PHat",
	P.time().timeName(),
	P.db(),
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    P.mesh(),
    dimensionedTensor("PHat", dimensionSet(1,-1,-2,0,0,0,0), tensor::zero)
);

Info << "Creating fields for wave speeds ..." << endl;

// Continuum mechanics class
mechanics mech(mesh, Sf);

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

////////////////////////////////////////////////////////////
//- ALE: Declaration and initialization of ALE Objects
Info << "Creating ALE objects ..." << nl;
aleModel aleModel_(mechanicalProperties, mesh, pMesh, N);
aleModel_.correct(); 
pointVectorField w_    = aleModel_.w();
pointTensorField aleF = aleModel_.defGrad(); 
pointScalarField aleJ = aleModel_.jacobian();
pointTensorField aleH = aleModel_.coFactor();

pointTensorField aleFinv = aleModel_.invDefGrad();
pointScalarField aleInvJ = aleModel_.invJacobian();
pointTensorField aleFinvT = aleModel_.defGradinvT();

//- ALE: Initialization of trueF and Piola Hat
trueF = F & aleFinv;
PHat = (P + ((aleInvJ * lm) * w_) ) & aleH;

//- ALE: Initialization of new wave speeds 
// Note : that part is not useful anymore ?
//- ALE: Longitudinal wave speed
pointScalarField aleUp
(
    IOobject("aleUp", mesh),
    pMesh,
    model.Up()/beta
);

//- ALE: Shear wave speed
pointScalarField aleUs
(
    IOobject("Us", mesh),
    pMesh,
    model.Us()*beta
);

//- ALE: initialization of wave speeds
scalarField lambdaH = sqrt((inter.pointToEdge(aleFinvT) & N) & (inter.pointToEdge(aleFinvT) & N));
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


/*
pointScalarField aleJ = aleModel_.jacobian();  Info << "aleJ created\n;";
pointScalarField aleInvJ = aleModel_.invJacobian(); Info << "aleInvJ created\n;";
pointTensorField aleFinvT = aleModel_.defGradinvT();//inv(aleModel_.defGrad().T());
Info << "aleFinvT created\n;";
pointTensorField aleF = aleModel_.defGrad(); Info << "aleF created\n;";
pointTensorField aleH = aleModel_.coFactor(); Info << "aleH created\n;";
scalarField lambdaH = sqrt((inter.pointToEdge(aleFinvT) & N) & (inter.pointToEdge(aleFinvT) & N)); Info << "lambdaH created\n;";
pointScalarField lambdaH_(IOobject("lambdaH_", mesh), pMesh, 0.0); Info << "lambdaH_ created\n;";
lambdaH_.primitiveFieldRef() = lambdaH; 
pointVectorField w_    = aleModel_.w(); Info << "w_ created\n;";
pointVectorField p_    = aleModel_.aleP(); Info << "p_ created\n;";

pointScalarField aleUp = aleInvJ; Info << "aleUp created\n;";
pointScalarField aleUs = aleInvJ; Info << "aleUs created\n;";
forAll(aleUp, pt)
{
  aleUp[pt] *= (lambdaH_[pt] * Up[pt]) - (w_[pt] & (aleF[pt] & N[pt] ));
  aleUs[pt] *= (lambdaH_[pt] * Us[pt]) - (w_[pt] & (aleF[pt] & N[pt] ));
}
Info << "aleUs and aleUp initialized\n;";
*/


/*
pointVectorField aleRecLM
(
 IOobject("aleRecLM", mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
 pMesh,
 dimensionedVector("aleRecLM", dimensionSet(1,-2,-1,0,0,0,0), vector::zero)
 
 //     This odd dimension is required because "lm" has the same one
  
);
Info << aleRecLM.dimensions() << nl;
Info << lm.dimensions() << nl;
Info << aleJ.dimensions() <<nl;
aleRecLM = lm * aleJ;
Info << "reconstructed Linear Momentum initialized\n";
*/












// Stretch
pointScalarField stretch = mech.stretch();

// Wave speed for time increment
pointScalarField Up_time
(
    IOobject("Up_time", mesh),
    pMesh,
    model.Up()
);
Up_time = Up_time/stretch;


Info << "Creating fields for gradient ..." << endl;

// Gradient class
gradientSchemes grad(mesh, Sf, V);

// Gradient of nodal linear momentum
tensorField lmGrad = grad.gradient(lm);

// Reconstruction of nodal linear momentum
vectorField lmM(mesh.edges().size());
vectorField lmP(mesh.edges().size());
grad.reconstruct(lm, lmGrad, lmM, lmP);


Info << "Creating fields for fluxes ..." << endl;

//- ALE: v
pointVectorField v_
(
	IOobject("v_", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
	pMesh,
	dimensionedVector("v_", dimensionSet(0,1,-1,0,0,0,0), vector::zero)
);
v_ = lm/rho;

//- ALE:  vHat
pointVectorField vHat
(
	IOobject("vHat", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
	pMesh,
	dimensionedVector("vHat", dimensionSet(0,1,-1,0,0,0,0), vector::zero)
);
vHat = (aleInvJ * v_) + ((F & aleFinv) & w_);



// Edge linear momentum
vectorList vC = inter.pointToEdge(v_);

// Edge traction
vectorField tC = (inter.pointToEdge(PHat) & N) + (0.5*mech.SmatrixMod(aleUp___, aleUs___) & (lmP-lmM));

Info << "Creating fields for angular momentum ..." << endl;

// Angular momentum class
angularMomentum am(mesh, mechanicalProperties);

// RHS of linear momentum equation
pointVectorField rhsLm
(
    IOobject ("rhsLm", mesh),
    pMesh,
    dimensionedVector("rhsLm", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
);

// RHS of linear momentum equation (stage 1)
pointVectorField rhsLm1 = rhsLm;

// RHS of angular momentum equation
pointVectorField rhsAm
(
    IOobject ("rhsAm", mesh),
    pMesh,
    dimensionedVector("rhsAm", dimensionSet(1,-1,-2,0,0,0,0), vector::zero)
);


Info << "Creating fields for post-processing ..." << endl;

// Nodal displacements
pointVectorField u
(
    IOobject ("u", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    pMesh,
    dimensionedVector("u", dimLength, vector::zero)
);
// w Nodal displacements
pointVectorField uw
(
    IOobject ("uw", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    pMesh,
    dimensionedVector("uw", dimLength, vector::zero)
);

u.write();
X -= u;


Info << "Creating variables for time ..." << endl;

scalar tstep = 0;
dimensionedScalar deltaT("deltaT", dimTime, runTime.deltaTValue());
dimensionedScalar t("t", dimTime, 0.0);

if (timeStepping == "variable")
{
    deltaT = (cfl*h)/model.Up();
    runTime.setDeltaT(deltaT);
}

scalarList RKstage(2);
RKstage[0] = 0;
RKstage[1] = 1;


Info << "Printing data ..." << endl;

// Print material properties
model.printMaterialProperties();
aleModel_.printMaterialProperties();

Info << "\nPrinting global momentum ..." << endl;

// Print global linear and angular momentum
am.printGlobalMomentum(lm,x,V);

// Print centroid of mesh
dual.printPrimalMeshCentroid();
dual.printDualMeshCentroid();







