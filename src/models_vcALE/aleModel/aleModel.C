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
    const fvMesh& vm,
    const vectorList& N
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

    N_(N),
    
    motMap_
    (
     IOobject("motMap", mesh_),
     pMesh_,
     dimensionedVector("motMap", dimensionSet(0,1,0,0,0,0,0), vector::zero)
    ),

    w_
    (
     IOobject("w", mesh_),
     pMesh_,
     dimensionedVector("w", dimensionSet(0,1,-1,0,0,0,0), vector::zero)     
    ),
    
    defGrad_
    (
     IOobject("defGrad", mesh_),
     pMesh_,
     Foam::tensor::I
    ),

    aleRho_
    (
     IOobject("aleRho", mesh_),
     pMesh_,
     dimensionedScalar("aleRho", dimensionSet(1,-3,0,0,0,0,0), 1.0)
    ),

    //fictitiousMotionType_(dict.lookup("fictitiousMotionType")),
    
    // to remove
    aleUp_
    (
     jacobian()
    ),

    aleUs_
    (
     jacobian()
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
{
  
  if (fictitiousMotionType() == "sinusoid")
  {
 //XE = xe + 1./20. * 2.0 * sin(2*pi*xe/1.) * sin(2*pi*ye/6.) * sin(2*pi*t/T)
 //YE = ye + 3./5. * 1.5 * sin(2*pi*xe/1.) * sin(2*pi*ye/6.) * sin(4*pi*t/T)
    scalar t = mesh_.time().value();
    scalar T = 2.0;
    scalar pi = Foam::constant::mathematical::pi;
    forAll(motMap_, p)
    {
      //motion Mapping
      scalar x = mesh_.points()[p][0];//pMesh_[p][0];
      scalar y = mesh_.points()[p][1];//pMesh_[p][1];
      motMap_[p][0] = x + 1.0/10.0 * Foam::sin(2*pi*x)
	* Foam::sin(pi*y/3.0)
	* Foam::sin(2*pi*t/T);
      motMap_[p][1] = y + 9.0/10.0 * Foam::sin(2*pi*x)
	* Foam::sin(pi*y/3.0)
	* Foam::sin(4*pi*t/T);
      //Velocity
      w_[p][0] = pi/(5.0*T) * Foam::sin(2*pi*x)
	* Foam::sin(pi*y/3.0)
	* Foam::cos(2*pi*t/T);
      w_[p][1] = 2.0*9.0*pi/(5.0*T) * Foam::sin(2*pi*x)
	* Foam::sin(pi*y/3.0)
	* Foam::cos(4*pi*t/T) ;
      //deformation Gradient
      defGrad_[p] = tensor
	(
	 1 + pi/5.0 * Foam::cos(2*pi*x)
	 * Foam::sin(pi*y/3.0)
	 * Foam::sin(2*pi*t/T),
	 pi/30.0 * Foam::sin(2*pi*x)
	 * Foam::cos(pi*y/3.0)
	 * Foam::sin(2*pi*t/T),
	 0,
	 9.0*pi/5.0 * Foam::cos(2*pi*x)
	 * Foam::sin(pi*y/3.0)
	 * Foam::sin(4*pi*t/T),
	 1 + 3.0*pi/10.0 * Foam::sin(2*pi*x)
	 * Foam::cos(pi*y/3.0)
	 * Foam::sin(4*pi*t/T),
	 0,
	 0,
	 0,
	 0
	);
    }   
  }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void aleModel::printMaterialProperties()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
