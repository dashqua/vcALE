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
  const pointMesh pMesh_;
  
        //- Constitutive model
        word model_;

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

        //- Pressure wave speed
        dimensionedScalar Up_;

        //- Shear wave speed
        dimensionedScalar Us_;

  //- (for ALE) material normals
  const vectorList N_;
  
  //- ALE Motion Mapping
  pointVectorField motMap_;

  //- Time derivative of motMap
  pointVectorField w_;
  
  //- ALE Deformation Gradient
  pointTensorField defGrad_;

  //- ALE density
  pointScalarField aleRho_;

  //- ALE motion type
  //word fictitiousMotionType_;
  
  //- ALE pressure wave speed
  //dimensionedScalar aleUp_;
  pointScalarField aleUp_;
  
  //- ALE Shear wave speed
  //dimensionedScalar aleUs_;
  pointScalarField aleUs_;  
  
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
   const pointTensorField&,
   const dictionary&,
   const fvMesh&,
   const vectorList&
   );

  //- Destructor
  virtual ~aleModel();


  // Member functions

  // Access

  const dimensionedScalar& Up() const
  {
    return Up_;
  }

  const dimensionedScalar& Us() const
  {
    return Us_;
  }

  // ALE mappings

  //- Motion Mapping
  pointVectorField motMap()
  {
    return motMap_;
  }

  //- Time derivative of motMap
  pointVectorField w()
  {
    return w_;
  }
  
  //- Deformation Gradient
  pointTensorField defGrad()
  {
    return defGrad_;
  }

  //- Jacobian
  pointScalarField jacobian()
  {
    return det(defGrad_);
  }

  //- inverse of the Jacobian
  pointScalarField invJacobian()
  {
    pointScalarField jac = jacobian();
    forAll(jac, x){ jac[x] = 1./jac[x]; }
    return jac;
  }
  
  //- CoFactor
  pointTensorField coFactor()
  {
    return jacobian() * inv(defGrad_.T());
  }

  //- Density
  pointScalarField aleRho()
  {
    return aleRho_;
  }

  //- Momentum
  pointVectorField aleP()
  {
    return aleRho_ * w_;
  }
  /*
  word fictitiousMotionType()
  {
    return fictitiousMotionType_;
  }
  */  


  //- Lambda coefficient
  /*
  pointScalarField lambdaCoeff()
  {
  //return (defGrad() & mesh_.magSf()/mesh_.Sf()) & (defGrad() & mesh_.magSf()/mesh_.Sf());
    interpolationSchemes inter(mesh_);
    vectorField han = (inter.pointToEdge(coFactor()) & N_);
    scalarField ent = han & han;
    pointScalarField hanhan(ent);
    
    return aleUp_;//(inter.pointToEdge(coFactor()) & N_) & (inter.pointToEdge(coFactor()) & N_);

    // problem with H * N
    // because N cannot be reached somehow

    // solution: H would have to be interpolated
    // see updateVariables.H
  }
  */
  
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
