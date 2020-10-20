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
    const dictionary& dict,
    const fvMesh& vm,
    pointMesh& pMesh
)
:
    mesh_(vm),
    pMesh_(pMesh),

    model_(dict.lookup("aleModel")),
    
    motMap_
    (
     IOobject("motMap", mesh_),
     pMesh_,
     dimensionedVector("motMap", dimensionSet(0,1,0,0,0,0,0), vector::zero)
    ),

    w_
    (
     IOobject
     (
      "w",
      mesh_.time().timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     pMesh_,
     dimensionedVector("w", dimensionSet(0,1,-1,0,0,0,0), vector::zero)     
    ),
    
    defGrad_
    (
     IOobject("defGrad", mesh_),
     pMesh_,
     Foam::tensor::I
    ),

    J_
    (
     IOobject
     (
      "J",
      mesh_.time().timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     pMesh_,
     1
    ),

    H_
    (
     IOobject("H", mesh_),
     pMesh_,
     Foam::tensor::I
    ),

    fictitiousMotionType_(dict.lookup("fictitiousMotionType")),

    beta_( 0.0 ),
    T_( 0.0 )
    //T_(dict.lookup("T"))
    //const scalar& cfl = readScalar(controlDict.lookup("cfl"));
  
{
  dict.lookup("beta") >> beta_;
  dict.lookup("T") >> T_;
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //
aleModel::~aleModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void aleModel::correct()
{
    //Info << fictitiousMotionType() << " .. t=" << mesh_.time().value() << nl;
  /*
    To update:
      - motMap_
      - w_
      - F (defGrad_)
      - J_
      - H_
   */
 if (fictitiousMotionType() == "sinusoid_order2")
  {
    // XE = xe +       beta * sin(2*pi*xe/1.) * sin(2*pi*ye/6.) * sin(pi*t/T)
    // YE = ye + 5.0 * beta * sin(2*pi*xe/1.) * sin(2*pi*ye/6.) * sin(2*pi*t/T)
    scalar t = mesh_.time().value();
    scalar pi = Foam::constant::mathematical::pi;
    forAll(motMap_, p)
    {
      //Info << "correcting ALE Motion Mapping\n";
      scalar X = mesh_.points()[p][0];
      scalar Y = mesh_.points()[p][1];
      scalar Z = mesh_.points()[p][2];
      motMap_[p][0] =
	X +   beta_ * Foam::sin(2*pi*X) * Foam::sin(pi*Y/3.0) * Foam::sin(pi*t/T_);
      motMap_[p][1] =
	Y + 5*beta_ * Foam::sin(2*pi*X) * Foam::sin(pi*Y/3.0) * Foam::sin(2*pi*t/T_);
      motMap_[p][2] = Z;
      //Info << "correcting ALE Velocity\n";
      w_[p][0] =
	beta_*pi/T_ * Foam::sin(2*pi*X) * Foam::sin(2*pi*X) * Foam::sin(pi*Y/3.0) * Foam::sin(pi*Y/3.0) * Foam::cos(pi*t/T_);
      w_[p][1] =
	10.0*beta_*pi/T_ * Foam::sin(2*pi*X) * Foam::sin(pi*Y/3.0) * Foam::sin(2*pi*X) * Foam::sin(pi*Y/3.0) * Foam::cos(2*pi*t/T_);
      w_[p][2] = 0.0;
      //Info << "correcting ALE defGrad\n";
      defGrad_[p] = tensor
	(
	 1 + 4*beta_*pi * Foam::sin(2*pi*X) * Foam::cos(2*pi*X) * Foam::sin(pi*Y/3.0) * Foam::sin(pi*Y/3.0) * Foam::sin(pi*t/T_),
	   2*beta_*pi/3.0 * Foam::sin(2*pi*X) * Foam::sin(2*pi*X) * Foam::cos(pi*Y/3.0) * Foam::sin(pi*Y/3.0) * Foam::sin(pi*t/T_),
	 0,//
	  20.0*beta_*pi * Foam::cos(2*pi*X) * Foam::sin(2*pi*X) * Foam::sin(pi*Y/3.0) * Foam::sin(pi*Y/3.0) * Foam::sin(2*pi*t/T_),
	 1 + 10.0*beta_*pi/3.0 * Foam::sin(2*pi*X) * Foam::sin(2*pi*X) * Foam::cos(pi*Y/3.0) * Foam::sin(pi*Y/3.0) * Foam::sin(2*pi*t/T_),
	 0,//
	 0,
	 0,
	 1//
	);
      J_[p] = det(defGrad_[p]);
      H_[p] = J_[p] * inv(defGrad_[p]).T();
    }   
  }
  if (fictitiousMotionType() == "sinusoid")
  {
    // XE = xe +       beta * sin(2*pi*xe/1.) * sin(2*pi*ye/6.) * sin(pi*t/T)
    // YE = ye + 5.0 * beta * sin(2*pi*xe/1.) * sin(2*pi*ye/6.) * sin(2*pi*t/T)
    scalar t = mesh_.time().value();
    scalar pi = Foam::constant::mathematical::pi;
    forAll(motMap_, p)
    {
      //Info << "correcting ALE Motion Mapping\n";
      scalar X = mesh_.points()[p][0];
      scalar Y = mesh_.points()[p][1];
      scalar Z = mesh_.points()[p][2];
      motMap_[p][0] =
	X +   beta_ * Foam::sin(2*pi*X) * Foam::sin(pi*Y/3.0) * Foam::sin(pi*t/T_);
      motMap_[p][1] =
	Y + 5*beta_ * Foam::sin(2*pi*X) * Foam::sin(pi*Y/3.0) * Foam::sin(2*pi*t/T_);
      motMap_[p][2] = Z;
      //Info << "correcting ALE Velocity\n";
      w_[p][0] =
	beta_*pi/T_ * Foam::sin(2*pi*X) * Foam::sin(pi*Y/3.0) * Foam::cos(pi*t/T_);
      w_[p][1] =
	10.0*beta_*pi/T_ * Foam::sin(2*pi*X) * Foam::sin(pi*Y/3.0) * Foam::cos(2*pi*t/T_);
      w_[p][2] = 0.0;
      //Info << "correcting ALE defGrad\n";
      defGrad_[p] = tensor
	(
	 1 + 2*beta_*pi * Foam::cos(2*pi*X) * Foam::sin(pi*Y/3.0) * Foam::sin(pi*t/T_),
	   beta_*pi/3.0 * Foam::sin(2*pi*X) * Foam::cos(pi*Y/3.0) * Foam::sin(pi*t/T_),
	 0,//
	  10.0*beta_*pi * Foam::cos(2*pi*X) * Foam::sin(pi*Y/3.0) * Foam::sin(2*pi*t/T_),
	 1 + 5.0*beta_*pi/3.0 * Foam::sin(2*pi*X) * Foam::cos(pi*Y/3.0) * Foam::sin(2*pi*t/T_),
	 0,//
	 0,
	 0,
	 1//
	);
      J_[p] = det(defGrad_[p]);
      H_[p] = J_[p] * inv(defGrad_[p]).T();
    }   
  }  
  
  if (fictitiousMotionType() == "sinusoidOLD") //DEPRECATED
  {
 //XE = xe + 1./20. * 2.0 * sin(2*pi*xe/1.) * sin(2*pi*ye/6.) * sin(2*pi*t/T)
 //YE = ye + 3./5. * 1.5 * sin(2*pi*xe/1.) * sin(2*pi*ye/6.) * sin(4*pi*t/T)
    scalar t = mesh_.time().value();
    scalar T = 2.0;
    scalar pi = Foam::constant::mathematical::pi;
    forAll(motMap_, p)
    {
      //Info << "correcting ALE Motion Mapping\n";
      scalar x = mesh_.points()[p][0];//pMesh_[p][0];
      scalar y = mesh_.points()[p][1];//pMesh_[p][1];
      motMap_[p][0] = x + 1.0/10.0 * Foam::sin(2*pi*x)
	* Foam::sin(pi*y/3.0)
	* Foam::sin(2*pi*t/T);
      motMap_[p][1] = y + 9.0/10.0 * Foam::sin(2*pi*x)
	* Foam::sin(pi*y/3.0)
	* Foam::sin(4*pi*t/T);
      //Info << "correcting ALE Velocity\n";
      w_[p][0] = pi/(5.0*T) * Foam::sin(2*pi*x)
	* Foam::sin(pi*y/3.0)
	* Foam::cos(2*pi*t/T);
      w_[p][1] = 2.0*9.0*pi/(5.0*T) * Foam::sin(2*pi*x)
	* Foam::sin(pi*y/3.0)
	* Foam::cos(4*pi*t/T) ;
      //Info << "correcting ALE defGrad\n";
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
	 1
	);
      //Info << "correcting ALE Jacobian\n";
      J_[p] = det(defGrad_[p]);
      //Info << "J[p] = " << J_[p] << nl;
      //Info << "correcting ALE invJacobian\n";
      //invJ_[p] = 1.0 / J_[p];
    }   
  }
  //Info << "ALE correction done\n";
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void aleModel::printMaterialProperties()
{
  Info << "aleModel: " << model_ << nl;
  Info << "fictitious Motion Type: " << fictitiousMotionType() << nl;
  Info << "beta: " << beta_ << nl;
  Info << "T: " << T_ << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
