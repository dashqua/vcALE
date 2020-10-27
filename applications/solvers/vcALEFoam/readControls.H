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


// ALE - GCL Gradient method
//const word& GCLgradientMethod = controlDict.lookup("GCLgradientMethod");

//if (GCLgradientMethod != "gradientScheme" && GCLgradientMethod != "gradientEvaluation")
//{
//    FatalErrorIn("readControls.H")
//        << "Valid type entries are 'gradientScheme' or 'gradientEvaluation' for GCLgradientMethod"
//        << abort(FatalError);
//}

// ALE - GCL technique:
// 1. "withInvolution" = Full GCL coverage with involution
// 2. "withoutInvolution" = GCL coverage without involution
const word& GCLtechnique = mechanicalProperties.lookupOrDefault<word>("GCLtechnique", "withInvolution");
Info << "GCL technique: " << GCLtechnique << nl;

// Type of time stepping
const word& timeStepping = controlDict.lookup("timeStepping");

if (timeStepping != "constant" && timeStepping != "variable")
{
    FatalErrorIn("readControls.H")
        << "Valid type entries are 'constant' or 'variable' for timeStepping"
        << abort(FatalError);
}

// Courant-Friedrichs-Lewy Number
const scalar& cfl = readScalar(controlDict.lookup("cfl"));

if (cfl <= 0.0 || cfl > 1.0)
{
    FatalErrorIn("readControls.H")
        << "Valid type entries are '<= 1' or '> 0' for cfl"
        << abort(FatalError);
}

// Stabilisation parameter for near incompressibility
const scalar& beta =
    readScalar(fvSolution.lookup("incompressiblilityCoefficient"));

// Angular momentum conservation
const word& angularMomentumConservation =
    fvSolution.lookup("angularMomentumConservation");

if
(
    angularMomentumConservation != "yes" && angularMomentumConservation != "no"
)
{
    FatalErrorIn("readControls.H")
        << "Valid type entries are 'yes' or 'no' "
        << "for angularMomentumConservation"
        << abort(FatalError);
}
