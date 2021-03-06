Info << "Reading dictionaries ..." << endl;

// Mechanical properties
IOdictionary mechanicalProperties
(
     IOobject
    (
        "mechanicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Run parameters
IOdictionary runParameters
(
     IOobject
    (
        "runParameters",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Control dictionary
IOdictionary controlDict
(
    IOobject
    (
        "controlDict",
        runTime.system(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Finite volume solution dictionary
IOdictionary fvSolution
(
    IOobject
    (
        "fvSolution",
        runTime.system(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// papr Functions dictionary
IOdictionary paperFunctionsProperties (
  IOobject ("paperFunctionsProperties", runTime.constant(), mesh, 
	  IOobject::MUST_READ_IF_MODIFIED, IOobject::NO_WRITE)
);


Info << "Reading data from dictionaries ..." << endl;

// Boundary traction
const vector& traction =
    runParameters.lookupOrDefault<vector>("traction", vector::zero);
Info << "value of 'traction': " << traction << nl;

// ALE - GCL technique:
// 1. "withInvolution" = Full GCL coverage with involution
// 2. "withoutInvolution" = GCL coverage without involution
const word& GCLtechnique = mechanicalProperties.lookupOrDefault<word>("GCLtechnique", "withoutInvolution");
if (GCLtechnique != "withoutInvolution" && GCLtechnique != "withInvolution") { 
    FatalErrorIn("readControls.H") << "Invalid keyword for 'GCLtechnique' in constant/mechanicalProperties."
		                       << abort(FatalError);
}
Info << "GCL technique: " << GCLtechnique << nl;

// ALE - PHatTerm: with or without involution?
const word& PHatTerm  = mechanicalProperties.lookupOrDefault<word>("PHatTerm", "withoutInvolution");
if (PHatTerm != "withoutInvolution" && PHatTerm != "withInvolution") {
  FatalErrorIn("readControls.H") << "Invalid keyword for 'PHatTerm' in constant/mechanicalProperties."
	                         << abort(FatalError); 
}
Info << "PHat term: " << PHatTerm << nl;

// Addition of an option to do postProcessing:
const word& doPostPro = mechanicalProperties.lookupOrDefault<word>("doPostPro", "activated");
Info << "Post Processing: " << doPostPro << nl;

// Type of time stepping
const word& timeStepping = controlDict.lookup("timeStepping");
if (timeStepping != "constant" && timeStepping != "variable") {
    FatalErrorIn("readControls.H") << "Valid type entries are 'constant' or 'variable' for timeStepping"
                                   << abort(FatalError);
}

// Stabilisation parameter for near incompressibility
const scalar& beta =
    readScalar(fvSolution.lookup("incompressiblilityCoefficient"));

// ALE - For P*: use of manufactured motion or use of P*-based RHS for solvedW
bool usePstar = mechanicalProperties.lookupOrDefault<bool>("usePstar", false);
Info << "Computation of solvedW's RHS: use of Pstar (1) or wDot (0): " << usePstar << nl;
const scalar& muMaterial_ = readScalar(mechanicalProperties.lookup("muMaterial"));
const scalar& muSpatial_  = readScalar(mechanicalProperties.lookup("muSpatial")) ;
Info << "Coefficient for Pstar -> muMaterial: " << muMaterial_ << " ___ muSpatial: " << muSpatial_ << nl;


// Courant-Friedrichs-Lewy Number
const scalar& cfl = readScalar(controlDict.lookup("cfl"));
if (cfl <= 0.0 || cfl > 1.0) {
    FatalErrorIn("readControls.H") << "Valid type entries are '<= 1' or '> 0' for cfl"
                                   << abort(FatalError);
}


// ALE - For update of solvedW
// update material motion every x iteration
//const scalar updateEvery = mechanicalProperties.lookupOrDefault<scalar>("updateEvery", 1);
//const scalar& updateEvery = readScalar(mechanicalProperties.lookup("updateEvery"));
//scalar curTimeStep = 0;



// Number of Runge-Kutta Stages
//  Note: Current possible options are 2 or 3
// int nbRKstages = 3;// 2 or 3;
const scalar& nbRKstages_ = readScalar(controlDict.lookup("RungeKuttaStages"));
const int nbRKstages(nbRKstages_);
/*
Info << "nbRKstages  : " << nbRKstages  << nl;
if (nbRKstages != 2 || nbRKstages != 3) {
    FatalErrorIn("readControls.H") << "'RungeKuttaStages' has to be '2' or '3'." << abort(FatalError);
}
*/

// Force Roller BCs for solvedW ?
bool rollerBCforW = mechanicalProperties.lookupOrDefault<bool>("rollerBCforW", false);
Info << "Do we apply roller boundary conditions for w ? " << (rollerBCforW ? "yes" : "no") << nl;
//
//
//
// To remove
//scalar tauBC  = 0; mechanicalProperties.lookup("tauBC") >> tauBC;
//scalar uMaxBC = 0; mechanicalProperties.lookup("uMaxBC") >> uMaxBC;
