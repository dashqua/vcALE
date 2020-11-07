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
    Foam::solidModel

Description
    Defines the solid model.

SourceFiles
    solidModel.C

\*---------------------------------------------------------------------------*/

#ifndef solidModel_H
#define solidModel_H

#include "fvMesh.H"
#include "pointFields.H"
#include "operations.H"
// #include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class solidModel
{
    // Private data

        // Mesh
        const fvMesh& mesh_;   

        //- First Piola Kirchhoff stress tensor
        pointTensorField P_;

        //- Pressure
        pointScalarField p_;

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

        // operations
        operations op;  

        // Von Mises vars
        //- Elastic left strain tensor
        pointTensorField b_;
  
        //- Hardening modulus
        dimensionedScalar Hm_ = dimensionedScalar("Hm", dimensionSet(1,-1,-2,0,0,0,0), 0.0);

        //- Initial yield stress
        dimensionedScalar Ys0_ = dimensionedScalar("Ys", dimensionSet(0,-1,-2,0,0,0,0), 0.0);

        //- Current yield stress
        pointScalarField Ys_;

        //- Plastic strain
        pointScalarField strain_p_;

        //- Inverse of plastic Right Cauchy Green strain tensor
        pointTensorField CpInv_;

        //- First Piola Kirchhoff stress tensor
        pointTensorField tau_;

        //- Von Mises stresses
        pointScalarField vMises_;
  
private:

        //- Disallow default bitwise copy construct
        solidModel(const solidModel&);

        //- Disallow default bitwise assignment
        void operator=(const solidModel&);

public:

    //- Runtime type information
    TypeName("solidModel");


    // Constructors

        //- Construct from field and dictionary
        solidModel
        (
            const pointTensorField&,
            const dictionary&,
            const fvMesh&
        );


    //- Destructor
    virtual ~solidModel();


    // Member functions

        // Access

            const pointTensorField& piola() const { return P_; }

            const pointScalarField& pressure() const { return p_; }

            const dimensionedScalar& density() const { return rho_; }

            const dimensionedScalar& lambda() const { return lambda_; }

            const dimensionedScalar& mu() const { return mu_; }

            const dimensionedScalar& nu() const { return nu_; }

            const dimensionedScalar& kappa() const { return kappa_; }

            const dimensionedScalar& Up() const { return Up_; }

            const dimensionedScalar& Us() const { return Us_; }

        // Edit

            //- Correct material model
            virtual void correct(pointTensorField&, pointTensorField&, pointScalarField&);

        // Print

            //- Print material properties
            void printMaterialProperties();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
