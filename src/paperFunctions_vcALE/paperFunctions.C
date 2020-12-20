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

\*---------------------------------------------------------------------------*/

#include "paperFunctions.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(paperFunctions, 0);


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

paperFunctions::paperFunctions
(
    const dictionary& dict,
    const dictionary& paperdict,
    const fvMesh& vm,
    pointMesh& pMesh
)
:
    mesh_(vm),
    pMesh_(pMesh),    
    rho_(dict.lookup("rho")),
    E_(dict.lookup("E")),
    nu_(dict.lookup("nu")),
    mu_(E_/(2.0*(1.0 + nu_))),
    lambda_(nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))),
    kappa_(lambda_ + (2.0/3.0)*mu_),

    FbarFunction(paperdict.lookup("FbarFunction"))
      
{ 
    Info << "paperFunction:" << nl
	 << "-> FbarFunction: " << FbarFunction <<nl;
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //
paperFunctions::~paperFunctions() {}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

pointTensorField paperFunctions::Fbar (pointTensorField& F){
    pointTensorField Fbar = F;
    if (FbarFunction == "Function1"){
      forAll(mesh_.points(), n) {
	Fbar[n] = tensor(1,0,0,
			 0,F[n].yy(),0,
			 0,0,1);
      }
    }
    if (FbarFunction == "Function2"){
      forAll(mesh_.boundary(), patch) {
      const fvPatch& patx = mesh_.boundary()[patch];
      if (patx.name().find("roller") != string::npos) {

        forAll(mesh_.boundaryMesh()[patch], facei) {
          const label& face = mesh_.boundaryMesh()[patch].start() + facei;
          forAll(mesh_.faces()[face], nodei) {
            const label& node = mesh_.faces()[face][nodei];
            Fbar[node]      = tensor(1,0,0,
			             0,F[node].yy(),0,
				     0,0,1
			    );
          }
        }
      }
      }
    }
    return Fbar;
  }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
}

// ************************************************************************* //
