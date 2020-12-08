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

Class
    Foam::paperFunctions

Description
    .

SourceFiles
    paperFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef paperFunctions_H
#define paperFunctions_H

#include "fvMesh.H"
#include "pointFields.H"
#include "vectorList.H"
#include "operations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class paperFunctions
{
  // Mesh
  const fvMesh& mesh_;

  // pointMesh
  const pointMesh& pMesh_;

  //- Density
  dimensionedScalar rho_;

  //- Young's modulus of Elasticity
  dimensionedScalar E_;

  //- Poisson's ratio
  dimensionedScalar nu_;

  //- Lame coefficient
  dimensionedScalar mu_;

  //- Lame coefficient
  dimensionedScalar lambda_;

  //- Bulk modulus
  dimensionedScalar kappa_;
  
private:

        //- Disallow default bitwise copy construct
        paperFunctions(const paperFunctions&);

        //- Disallow default bitwise assignment
        void operator=(const paperFunctions&);

public:

    //- Runtime type information
    TypeName("paperFunctions");


  //- Constructor
  paperFunctions ( const dictionary&, const dictionary&, const fvMesh&, pointMesh& );

  //- Destructor
  virtual ~paperFunctions();


  // Member functions

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
