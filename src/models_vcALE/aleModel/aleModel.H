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
#include "operations.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class aleModel
{
  const word name_;
  //const word& timeStepping = controlDict.lookup("timeStepping");

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

  //- ALE wDot;
  pointVectorField wDot_;
  
  //- ALE motion type
  word fictitiousMotionType_;
  
  //- ALE beta coefficient for Motion Mapping
  scalar beta_;

  //- ALE XR and YR for fictitious motion
  scalar XR;
  scalar YR;

  //- ALE T period coefficient for Mot Map
  scalar T_;

  ////////////////////////////////////////////
  // Take quantities for solidModel for this fictitious model
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

  /////////////////////////////
  // Operations definition
  operations op;

  // in an attempt to save up some time
  bool usePstar_;
  
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
   pointMesh&,
   const word name
  );

  //- Destructor
  virtual ~aleModel();


  // Member functions

  // Access

  // ALE mappings

  //- Motion Mapping
  pointVectorField& motMap() { return motMap_; }

  //- Time derivative of motMap
  pointVectorField& w() { return w_; }

  // Compute wDot
  pointVectorField& wDot() { return wDot_; }
  
  //- Deformation Gradient
  pointTensorField& defGrad() { return defGrad_; }

  //- Jacobian
  pointScalarField& jacobian() { return J_; }  

  /////////////////////////////////////////////////////////////
  
  word fictitiousMotionType() { return fictitiousMotionType_; }

  //- Correct material model
  virtual void correct();

  //- Print material properties
  void printProperties();

  // get P_{Mat} from F_{Mat}, according to the chosen model
  pointTensorField piola (pointTensorField&, pointTensorField&, pointScalarField&);

  // get Pres_{Mat} from J_{Mat}, according to the chosen model
  //pointScalarField getMaterialPressure (pointTensorField&, pointTensorField&, pointScalarField&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
