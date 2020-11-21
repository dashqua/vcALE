/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    initialConditions_vc

Description
    Generates non-standard initial conditions for test cases.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pointFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // Point mesh
    pointMesh pMesh(mesh);

    // Read mechanical properties dictionary
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

    // Read run parameters dictionary
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

    // Test case name
    const word& tutorial(runParameters.lookupOrDefault<word>("tutorial", "None"));
    Info << "Tutorial name: " << tutorial;

    // Read density
    const dimensionedScalar& rho(mechanicalProperties.lookup("rho"));

    // Read energy
    const dimensionedScalar& E(mechanicalProperties.lookup("E"));

    // Read nu
    const dimensionedScalar& nu(mechanicalProperties.lookup("nu"));

    // Compute mu
    const dimensionedScalar mu = E/(2*(1+nu));

    // Nodal coordinates
    const vectorField& X = mesh.points();

    // Read linear momentum field
    pointVectorField lm
    (
        IOobject
        (
            "lm",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh,
        dimensionedVector("lm", dimensionSet(1,-2,-1,0,0,0,0), vector::zero)
    );

    if (tutorial == "bendingColumnLM")
    {
    	forAll(mesh.points(), node)
	{
		vector res(1,0,0);	
		lm[node] = rho.value() * X[node].y()*res;
	}
    }

    // Non uniform angular velocity initialsed testcases
    else if (tutorial == "twistingColumn" || tutorial == "spinningCube" || tutorial == "spinningTop" )
    {
        pointVectorField omega
        (
            IOobject("omega", mesh),
            pMesh,
            dimensionedVector(runParameters.lookup("initialAngularVelocity"))
        );

        if (tutorial == "twistingColumn")
        {
            const scalar& PI = Foam::constant::mathematical::pi;
            dimensionedScalar height("height", dimensionSet(0,1,0,0,0,0,0), 6.0);

            forAll(mesh.points(), node)
            {
                //lm[node] =
                //    rho.value()*(omega[node]
                //   *Foam::sin(PI*X[node].y()/(2*height.value()))) ^ X[node];

                omega[node] *= Foam::sin(PI*X[node].y()/(2*height.value()));
            }
        }

        forAll(mesh.points(), node)
        {
            lm[node] = rho.value()*omega[node] ^ X[node];
        }
    }

    else if (tutorial == "taylorImpact") {
        // Read initial velocity
	pointVectorField v (
	  IOobject("v", mesh), 
	  pMesh,
	  dimensionedVector(runParameters.lookup("initialVelocity"))
	);
	lm = rho*v;
    } else if (tutorial == "punchTest") {   
        // Read initial velocity
	pointVectorField v ( IOobject("v", mesh), pMesh, dimensionedVector(runParameters.lookup("initialVelocity")) );
	forAll(mesh.points(), n){
	    const scalar X = mesh.points()[n][0];
	    const scalar Y = mesh.points()[n][1];
	    //Info << X << " , " << Y << nl;
	    if ((0.0<X) && (0.0<Y)) {
		lm[n] = rho.value() * v[n];
	    }
	}
    } else if (tutorial == "swingingCube") {
	// Read parameters
	const scalar& A  = readScalar(runParameters.lookup("A"));
	const scalar& B  = readScalar(runParameters.lookup("B"));
	const scalar& C  = readScalar(runParameters.lookup("C"));
	pointScalarField U0 (
	         IOobject("U0", mesh),
	         pMesh,
	         dimensionedScalar(runParameters.lookup("U0"))
	);

	//const scalar cd  = Foam::sqrt(mu/rho); // useless
	const scalar pi  = Foam::constant::mathematical::pi;

	// Initialization of velocity
	forAll(mesh.points(), n){
          const scalar X = mesh.points()[n][0];
	  const scalar Y = mesh.points()[n][1];
	  const scalar Z = mesh.points()[n][2];
          lm[n][0] = U0[n] * A * Foam::sin(pi*X/2.) * Foam::cos(pi*Y/2.) * Foam::cos(pi*Z/2.);
	  lm[n][1] = U0[n] * B * Foam::cos(pi*X/2.) * Foam::sin(pi*Y/2.) * Foam::cos(pi*Z/2.);
	  lm[n][2] = U0[n] * C * Foam::cos(pi*X/2.) * Foam::cos(pi*Y/2.) * Foam::sin(pi*Z/2.);
	}
    } else {
	Info << "Initial Conditions: nothing is applied." << nl;
    }

    lm.write();

    Info<< "\n end\n";

    return 0;
}
