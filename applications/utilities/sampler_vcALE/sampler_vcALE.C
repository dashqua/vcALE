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
    sampler.C

Description
    Samples.

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
    Info << "Tutorial name: " << tutorial << nl;



    while (runTime.loop())  {
	    if (mesh.time().outputTime()){
        Info << "t = " << mesh.time().value() ; 
        pointVectorField u(
           IOobject("u", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), pMesh 
        );
	forAll(mesh.points(), n) {
	  scalar X = mesh.points()[n][0],   Y =  mesh.points()[n][1],     Z =  mesh.points()[n][2];
	  if ( X == 0.006413 && Y == 0 && Z == 0 ) {
	       Info << " - uBOT = " << u[n] ;
	  }
	  if ( X == 0 && Y == 0.026667 && Z == 0 ) {
	       Info << " - uTOP = " << u[n] ;
	  }
	}
        Info << nl;	
	    }
    }
    return 0;
}
