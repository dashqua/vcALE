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


Info << "Reading data from dictionaries ..." << endl;

// Boundary traction
const vector& traction =
    runParameters.lookupOrDefault<vector>("traction", vector::zero);

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

const scalar& alpha = mechanicalProperties.lookupOrDefault<dimensionedScalar>("alpha", 0.0).value();
Info << "Average P* coefficient alpha: " << alpha << nl;

// Courant-Friedrichs-Lewy Number
const scalar& cfl = readScalar(controlDict.lookup("cfl"));
if (cfl <= 0.0 || cfl > 1.0) {
    FatalErrorIn("readControls.H") << "Valid type entries are '<= 1' or '> 0' for cfl"
                                   << abort(FatalError);
}

// Stabilisation parameter for near incompressibility
const scalar& beta =
    readScalar(fvSolution.lookup("incompressiblilityCoefficient"));

