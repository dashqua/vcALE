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
//const polyBoundaryMesh& bm = mesh.boundaryMesh();
//const label& freePatchID = bm.findPatchID("free");
//const label& tractionPatchID = bm.findPatchID("traction");
//const label& fixedPatchID = bm.findPatchID("fixed");

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
pointVectorField X( IOobject("X", mesh), x);

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


// ALE - Material J (to be solved)
pointScalarField matJ
(
	IOobject ("matJ", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
	pMesh,
	dimensionedScalar("matJ", dimensionSet(0,0,0,0,0,0,0), 1)
);
matJ.write();

// ALE - Material F (to be solved)
pointTensorField matF
(
        IOobject("matF", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
	pMesh,
	dimensionedTensor("mat", F.dimensions(), Foam::tensor::I)
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
aleModel aleModel_(mechanicalProperties, mesh, pMesh);
aleModel_.correct();
pointVectorField& w_    = aleModel_.w();
//pointTensorField& aleF = aleModel_.defGrad(); 
//pointScalarField& aleJ = aleModel_.jacobian();
pointTensorField& aleH = aleModel_.coFactor();

//- ALE: Initialization of trueF and Piola Hat
trueF = F & op.inverse(matF);               // NOTE : this aleF has to be ideally replaced by matF


//- ALE: initialization of wave speeds
scalarField lambdaH = sqrt((inter.pointToEdge(aleH) & N) & (inter.pointToEdge(aleH) & N));
scalarField aleUp___ = inter.pointToEdge( op.inverseScalar(matJ) ) * 
       (
         (lambdaH * inter.pointToEdge(Up)) 
       - (inter.pointToEdge(w_) & (inter.pointToEdge(aleH) & N))
       );

scalarField aleUs___ = inter.pointToEdge( op.inverseScalar(matJ) ) * 
       (
         (lambdaH * inter.pointToEdge(Us)) 
       - (inter.pointToEdge(w_) & (inter.pointToEdge(aleH) & N))
       );













// Stretch
pointScalarField stretch = mech.stretch();

// Wave speed for time increment
pointScalarField Up_time
(
    IOobject("Up_time", mesh), pMesh, model.Up()
);
Up_time = Up_time/stretch;


Info << "Creating fields for fluxes ..." << endl;


// ALE: pR and pTilde, step quantities for rhsLM
pointVectorField pR
(
    IOobject ("pR", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    pMesh,
    dimensionedVector("pR", dimensionSet(1,-2,-1,0,0,0,0), vector::zero)
);
pR = op.inverseScalar(matJ) * lm;
pR.write();

pointVectorField pTilde
(
    IOobject ("pTilde", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    pMesh,
    dimensionedVector("pTilde", dimensionSet(1,-2,-1,0,0,0,0), vector::zero)
);
pTilde = matJ * pR;
pTilde.write(); 


//- ALE: v
pointVectorField v_
(
	IOobject("v_", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
	pMesh,
	dimensionedVector("v_", dimensionSet(0,1,-1,0,0,0,0), vector::zero)
);
v_ = pR/rho;

//- ALE:  vHat
pointVectorField vHat
(
	IOobject("vHat", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
	pMesh,
	dimensionedVector("vHat", dimensionSet(0,1,-1,0,0,0,0), vector::zero)
);
vHat = v_ + (trueF & w_);

// ALE: PHat
PHat = P & aleH ;



Info << "Creating fields for gradient ..." << endl;

// Gradient class
gradientSchemes grad(mesh, Sf, V);

// Gradient of nodal linear momentum
tensorField vGrad = grad.gradient(v_);

// Reconstruction of nodal linear momentum
vectorField vM(mesh.edges().size());
vectorField vP(mesh.edges().size());
grad.reconstruct(v_, vGrad, vM, vP);



// Edge traction
vectorField tC = (inter.pointToEdge(PHat) & N )
	+ inter.pointToEdge(matJ*rho)*(0.5*mech.SmatrixMod(aleUp___, aleUs___) & (vP - vM));
// Edge linear momentum
tensorField vC = inter.pointToEdge(vHat) * N;
// Edge traction two
tensorField tC2 = inter.pointToEdge( pR ) * N;
// Edge material velocity
tensorField wC = inter.pointToEdge( w_ )  * N;
// Edge alternative jacobian term
scalarField HTwC = inter.pointToEdge(aleH.T()&w_) & N;
// ALE - The GCL technique distinction (commented code below) could
// be done, but for the sake of simplicity, initial values for both
// cases are declared anyway.
// It is also about not cluttering the code with IFs verifications.
/*
if (GCLtechnique == "withInvolution") {
    tC2 = inter.pointToEdge( pR )   * N;
} else if (GCLtechnique == "withoutInvolution") {
    HTwC = inter.pointToEdge( aleH.T()&w_ ) & N;
} else {
    FatalErrorIn("updateVariables.H") << "GCL technique is not properly defined." << abort(FatalError);
}
*/




// ALE - RHS of jacobian equation
pointScalarField rhsJ
(
	IOobject("rhsJ", mesh),
	pMesh,
	dimensionedScalar("rhsJ", dimensionSet(0,0,-1,0,0,0,0), 0)
);


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
pointVectorField rhsLm1
(
    IOobject ("rhsLm1", mesh),
    pMesh,
    dimensionedVector("rhsLm1", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
);
pointVectorField rhsLm2
(
    IOobject ("rhsLm2", mesh),
    pMesh,
    dimensionedVector("rhsLm2", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
);
pointVectorField rhsLm3
(
    IOobject ("rhsLm3", mesh),
    pMesh,
    dimensionedVector("rhsLm3", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
);

// ALE - Here: "IF statement" like in gEqns.H to attribute value to rhsLm*
rhsLm1.primitiveFieldRef() = op.surfaceSum(tC*magSf);
if (GCLtechnique == "withInvolution") {
    rhsJ.primitiveFieldRef() = aleH && op.surfaceSum(wC*magSf);
    rhsLm2.primitiveFieldRef() = op.surfaceSum(tC2*magSf) & (aleH.T() & w_);
    rhsLm3.primitiveFieldRef() = pR * (aleH && (op.surfaceSum(wC*magSf)));
    rhsLm.primitiveFieldRef() = rhsLm1 + rhsLm2  + rhsLm3;
} else if (GCLtechnique == "withoutInvolution") {
    rhsJ.primitiveFieldRef() = op.surfaceSum(HTwC*magSf);
    rhsLm2.primitiveFieldRef() = pR.primitiveFieldRef() * rhsJ.primitiveFieldRef();
    rhsLm.primitiveFieldRef() = rhsLm1 + rhsLm2;
} else {
    FatalErrorIn("createFields.H") << "GCL technique is not properly defined." << abort(FatalError);
}

// ALE - RHS of material deformation gradient equation
pointTensorField rhsMatF
(
    IOobject("rhsMatF", mesh),
    pMesh,
    dimensionedTensor("rhsMatF", F.dimensions()/dimTime, tensor::zero)
);

// RHS of linear momentum equation (stage 1)
//pointVectorField rhsLmBis = rhsLm;

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
uw.write();
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







