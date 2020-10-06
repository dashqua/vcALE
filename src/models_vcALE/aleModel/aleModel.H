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
    Foam::aleModel

Description
    Defines the ale model.

SourceFiles
    aleModel.C

\*---------------------------------------------------------------------------*/

#ifndef aleModel_H
#define aleModel_H

#include "fvMesh.H"
#include "pointFields.H"
#include "vectorList.H"
#include "interpolationSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class aleModel
{
    // Private data

        // Mesh
         const fvMesh& mesh_;

  // pointMesh
  const pointMesh& pMesh_;
  
        //- Constitutive model
        word model_;

  //- ALE Motion Mapping
  pointVectorField motMap_;

  //- Time derivative of motMap
  pointVectorField w_;
  
  //- ALE Deformation Gradient
  pointTensorField defGrad_;

  //- ALE Jacobian
  pointScalarField J_;

  //- ALE coFactor
  pointTensorField H_;
  
  //- ALE motion type
  word fictitiousMotionType_;
  
  //- ALE beta coefficient for Motion Mapping
  scalar beta_;

  //- ALE T period coefficient for Mot Map
  scalar T_;

  
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        aleModel(const aleModel&);

        //- Disallow default bitwise assignment
        void operator=(const aleModel&);

public:

    //- Runtime type information
    TypeName("aleModel");


  //- Constructor
  aleModel
  (
   const dictionary&,
   const fvMesh&,
   pointMesh&
  );

  //- Destructor
  virtual ~aleModel();


  // Member functions

  // Access

  // ALE mappings

  //- Motion Mapping
  pointVectorField& motMap()
  {
    return motMap_;
  }

  //- Time derivative of motMap
  pointVectorField& w()
  {
    return w_;
  }
  
  //- Deformation Gradient
  pointTensorField& defGrad()
  {
    return defGrad_;
  }

  //- Jacobian
  pointScalarField& jacobian()
  {
    return J_; 
  }
  
  //- CoFactor
  pointTensorField& coFactor()
  {
    pointScalarField& jac = jacobian();
    pointTensorField& dg  = defGrad();
    pointTensorField& cof = defGrad();
    forAll(cof, pt)
    {
   	cof[pt] = jac[pt] * inv( dg[pt].T() );
    }
    return cof; // jac * inv( dg.T() );
  }

  word fictitiousMotionType()
  {
    return fictitiousMotionType_;
  }
      
  // Edit

  //- Correct material model
  virtual void correct();

  // Print

  //- Print material properties
  void printMaterialProperties();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
