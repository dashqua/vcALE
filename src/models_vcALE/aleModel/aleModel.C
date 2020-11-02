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
    
    motMap_(
     IOobject("motMap", mesh_),
     pMesh_,
     dimensionedVector("motMap", dimensionSet(0,1,0,0,0,0,0), vector::zero)
    ),

    w_(
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
    
    defGrad_(
     IOobject("defGrad", mesh_),
     pMesh_,
     Foam::tensor::I
    ),

    J_(
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

    H_(
     IOobject("H", mesh_),
     pMesh_,
     Foam::tensor::I
    ),

    wDot_(
     IOobject
     (
      "wDot",
      mesh_.time().timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
     ),
     pMesh_,
     dimensionedVector("wDot", dimensionSet(0,1,-2,0,0,0,0), vector::zero)     
    ),
    
    fictitiousMotionType_(dict.lookup("fictitiousMotionType")),

    beta_( 0.0 ),
    T_( 0.0 ),
    
    rho_(dict.lookup("rho")),
    E_(dict.lookup("E")),
    nu_(dict.lookup("nu")),
    mu_(E_/(2.0*(1.0 + nu_))),
    lambda_(nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))),
    kappa_(lambda_ + (2.0/3.0)*mu_),

    op(mesh_)
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
 if (fictitiousMotionType() == "sinusoid_order2inTime") {
    // XE = xe +       beta * sin(2*pi*xe/1.) * sin(2*pi*ye/6.) * sin(pi*t/T)
    // YE = ye + 5.0 * beta * sin(2*pi*xe/1.) * sin(2*pi*ye/6.) * sin(2*pi*t/T)
    scalar t = mesh_.time().value();
    scalar pi = Foam::constant::mathematical::pi, pi2 = pi*pi, T2=T_*T_;
    forAll(motMap_, p)
    {
      //Info << "correcting ALE Motion Mapping\n";
      scalar X = mesh_.points()[p][0];
      scalar Y = mesh_.points()[p][1];
      scalar Z = mesh_.points()[p][2];
      motMap_[p][0] =
	X +   beta_ * Foam::sin(2*pi*X) * Foam::sin(pi*Y/3.0) * Foam::sin(pi*t/T_)*Foam::sin(pi*t/T_);
      motMap_[p][1] =
	Y + 5*beta_ * Foam::sin(2*pi*X) * Foam::sin(pi*Y/3.0) * Foam::sin(2*pi*t/T_)*Foam::sin(2*pi*t/T_);
      motMap_[p][2] = Z;
      //Info << "correcting ALE Velocity\n";
      w_[p][0] =
	2*beta_*pi/T_ * Foam::sin(2*pi*X) * Foam::sin(pi*Y/3.0) * Foam::sin(pi*t/T_) * Foam::cos(pi*t/T_);
      w_[p][1] =
	20.0*beta_*pi/T_ * Foam::sin(pi*Y/3.0) * Foam::sin(2*pi*X) * Foam::sin(pi*Y/3.0) * Foam::cos(2*pi*t/T_) * Foam::sin(2*pi*t/T_);
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
      wDot_[p][0] = 2*beta_*pi2/T2 * Foam::sin(2*pi*X) * Foam::sin(pi*Y/3.0)
	* ((Foam::cos(pi*t/T_)*Foam::cos(pi*t/T_)) - (Foam::sin(pi*t/T_)*Foam::sin(pi*t/T_)));
      wDot_[p][1] = 40*beta_*pi2/T2 * Foam::sin(2*pi*X) * Foam::sin(pi*Y/3.0)
	* ((Foam::cos(2*pi*t/T_)*Foam::cos(2*pi*t/T_)) - (Foam::sin(2*pi*t/T_)*Foam::sin(2*pi*t/T_)));
      wDot_[p][2] = 0; 
      
    }   
} else if (fictitiousMotionType() == "sinusoid_order2") {
    // XE = xe +       beta * sin(2*pi*xe/1.) * sin(2*pi*ye/6.) * sin(pi*t/T)
    // YE = ye + 5.0 * beta * sin(2*pi*xe/1.) * sin(2*pi*ye/6.) * sin(2*pi*t/T)
    scalar t = mesh_.time().value();
    scalar pi = Foam::constant::mathematical::pi, pi2 = pi*pi, T2=T_*T_;
    forAll(motMap_, p)
    {
      //Info << "correcting ALE Motion Mapping\n";
      scalar X = mesh_.points()[p][0];
      scalar Y = mesh_.points()[p][1];
      scalar Z = mesh_.points()[p][2];
      motMap_[p][0] =
	X +   beta_ * Foam::sin(2*pi*X) * Foam::sin(2*pi*X) * Foam::sin(pi*Y/3.0) * Foam::sin(pi*Y/3.0) * Foam::sin(pi*t/T_);
      motMap_[p][1] =
	Y + 5*beta_ * Foam::sin(2*pi*X) * Foam::sin(2*pi*X) * Foam::sin(pi*Y/3.0) * Foam::sin(pi*Y/3.0) * Foam::sin(2*pi*t/T_);
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
      wDot_[p][0] = (X - motMap_[p][0]) * (pi2/T2);
      wDot_[p][1] = (Y - motMap_[p][1]) * (4*pi2/T2);
      wDot_[p][2] = 0; 
      
    }   
  } else if (fictitiousMotionType() == "sinusoid") {
    // XE = xe +       beta * sin(2*pi*xe/1.) * sin(2*pi*ye/6.) * sin(pi*t/T)
    // YE = ye + 5.0 * beta * sin(2*pi*xe/1.) * sin(2*pi*ye/6.) * sin(2*pi*t/T)
    scalar t = mesh_.time().value();
    scalar pi = Foam::constant::mathematical::pi, pi2 = pi*pi, T2=T_*T_;
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
      wDot_[p][0] = (X - motMap_[p][0]) * (pi2/T2);
      //-beta_*pi2/T2    * Foam::sin(2*pi*X) * Foam::sin(pi*Y/3.0) * Foam::sin(pi*t/T_) ; //(mesh_.points()[n][0] - wD[n][0]) * (pi2/T2) ;
      wDot_[p][1] = (Y - motMap_[p][1]) * (4*pi2/T2);
      //-20*beta_*pi2/T2 * Foam::sin(2*pi*X) * Foam::sin(pi*Y/3.0) * Foam::sin(2*pi*t/T_) ;//(mesh_.points()[n][1] - wD[n][0]) * (4*pi2/T2) ;
      wDot_[p][2] = 0; 
    }   
  } else {
   FatalErrorIn("aleModel.C") << "Problem in updating prescribed material motion." << abort(FatalError) ;
  }
 //
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

pointTensorField aleModel::getMaterialPiola (pointTensorField& matF, pointTensorField& matH, pointScalarField& matJ)
{
  pointTensorField matP
  (
     IOobject ("matP", mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::AUTO_WRITE),
     pMesh_,
     dimensionedTensor("matP", dimensionSet(1,-1,-2,0,0,0,0), tensor::zero)  // 1,-1,-2
  );
  
  forAll(mesh_.points(), n) {
    matP[n] = mu_.value()*pow(matJ[n], -2.0/3.0)*matF[n]
      - ((mu_.value()/3.0)*pow(matJ[n],(-5.0/3.0))*(matF[n] && matF[n])*matH[n])
      + ( kappa_.value()*(matJ[n]-1.0) )*matH[n];
  }
  
  return matP;
}

pointVectorField aleModel::bodyForces ()
{
  scalar pi = Foam::constant::mathematical::pi, pi2 = pi*pi, T2 = T_*T_;  
  pointVectorField w_dot = motMap();
  forAll(w_dot, n){
    w_dot[n][0] = (mesh_.points()[n][0] - w_dot[n][0]) * (pi2/T2);
    w_dot[n][1] = (mesh_.points()[n][1] - w_dot[n][1]) * (pi2/T2);
    w_dot[n][2] = 0;
  }

  //test -> op is recognized OK
  pointScalarField test = det(H_);
  return op.inverseScalar(test) * rho_  * w_dot ;// /!\ Not ready yet
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//pointVectorField& aleModel::wDot(){
  /*
  scalar pi = Foam::constant::mathematical::pi, pi2 = pi*pi, T2 = T_*T_;
  scalar t = mesh_.time().value(); 
  pointVectorField wD( 
     IOobject("wD", mesh_),
     pMesh_,
     dimensionedVector("wD", dimensionSet(0,1,-2,0,0,0,0), vector::zero)       // /!\ dimensions are wrong  0,1,-2
  );
 
  //wD.primitiveFieldRef() = motMap_.primitiveFieldRef();
  
  forAll(wD, n){
    scalar X = mesh_.points()[n][0];
    scalar Y = mesh_.points()[n][1];
    //scalar Z = mesh_.points()[n][2];    
    wD[n][0] = -beta_*pi2/T2    * Foam::sin(2*pi*X) * Foam::sin(pi*Y/3.0) * Foam::sin(pi*t/T_) ; //(mesh_.points()[n][0] - wD[n][0]) * (pi2/T2) ;
    wD[n][1] = -10*beta_*pi2/T2 * Foam::sin(2*pi*X) * Foam::sin(pi*Y/3.0) * Foam::sin(2*pi*t/T_) ;//(mesh_.points()[n][1] - wD[n][0]) * (4*pi2/T2) ;
    wD[n][2] = 0;
  }

  return wD;
  */
//  return wDot_;
//  
//}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
}

// ************************************************************************* //
