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

#include "aleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(aleModel, 0);


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

aleModel::aleModel
(
    const pointTensorField& F,
    const dictionary& dict,
    const fvMesh& vm
)
:
    mesh_(vm),
    pMesh_(mesh_),

    model_(dict.lookup("aleModel")),

    rho_(dict.lookup("rho")),
    E_(dict.lookup("E")),
    nu_(dict.lookup("nu")),
    mu_(E_/(2.0*(1.0 + nu_))),
    lambda_(nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))),
    kappa_(lambda_ + (2.0/3.0)*mu_),

    Up_(sqrt((lambda_+2.0*mu_)/rho_)),
    Us_(sqrt(mu_/rho_)),

    motMap_
    (
     IOobject("motMap", mesh_),
     pMesh_,
     dimensionedVector("motMap", dimensionSet(0,1,0,0,0,0,0), vector::zero)
    ),

    defGrad_
    (
     IOobject("defGrad", mesh_),
     pMesh_,
     Foam::tensor::I
    )
{

  // The structure of the solidModel is kept so far.
  // The new wave speeds are initialized here.

  //aleUp_ = inv(this->jacobian(pMesh_.points()));
  //aleUs_ = aleUp_;
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //
aleModel::~aleModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void aleModel::correct()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void aleModel::printMaterialProperties()
{
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
