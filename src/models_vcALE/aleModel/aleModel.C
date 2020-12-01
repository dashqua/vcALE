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
    pointMesh& pMesh,
    const word name
)
:
    name_(name),
    mesh_(vm),
    pMesh_(pMesh),

    model_(dict.subDict(name_).lookup("aleModel")),
    
    motMap_( IOobject("motMap", mesh_), pMesh_,
     dimensionedVector("motMap", dimLength, vector::zero)
    ),

    w_( IOobject("w", mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE), pMesh_,
     dimensionedVector("w", dimLength/dimTime, vector::zero)     
    ),
    
    defGrad_( IOobject("defGrad", mesh_), pMesh_, Foam::tensor::I ),

    J_( IOobject("J", mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE), pMesh_, 1),

    H_( IOobject("H", mesh_), pMesh_, Foam::tensor::I),

    wDot_( IOobject("wDot", mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE), pMesh_,
	   dimensionedVector("wDot", w_.dimensions()/dimTime, vector::zero)     
    ),
    
    fictitiousMotionType_("Void"),

    beta_( 0.0 ),
    XR(1.0),
    YR(6.0),
    T_(0.0),
    
    rho_(dict.lookup("rho")),
    E_(dict.lookup("E")),
    nu_(dict.lookup("nu")),
    mu_(E_/(2.0*(1.0 + nu_))),
    lambda_(nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))),
    kappa_(lambda_ + (2.0/3.0)*mu_),
    
    op(mesh_),
    usePstar_(dict.lookup("usePstar"))
  
{
  dict.lookup("usePstar") >> usePstar_;
  if (model_ == "neoHookean") {
    if (!usePstar_) {
      dict.subDict(name_).subDict("neoHookeanDict").lookup("fictitiousMotionType") >> fictitiousMotionType_;
      dict.subDict(name_).subDict("neoHookeanDict").lookup("beta") >> beta_;
      
      dict.subDict(name_).subDict("neoHookeanDict").lookup("XR") >> XR;
      dict.subDict(name_).subDict("neoHookeanDict").lookup("YR")>> YR;
    
      dict.subDict(name_).subDict("neoHookeanDict").lookup("T") >> T_;
    }
  }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //
aleModel::~aleModel() {}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void aleModel::correct()
{
 if (usePstar_) { return; }
 if (fictitiousMotionType() == "sinusoid_order2inTime") {
    scalar t = mesh_.time().value();
    scalar pi = Foam::constant::mathematical::pi, pi2 = pi*pi, T2=T_*T_;
    forAll(wDot_, p)
    {	/*    
      scalar X = mesh_.points()[p][0];
      scalar Y = mesh_.points()[p][1];
      wDot_[p][0] = 2*beta_*pi2/T2 * Foam::sin(2*pi*X/1.) * Foam::sin(2*pi*Y/6.)
	* ((Foam::cos(pi*t/T_)*Foam::cos(pi*t/T_)) - (Foam::sin(pi*t/T_)*Foam::sin(pi*t/T_)));
      wDot_[p][1] = 40*beta_*pi2/T2 * Foam::sin(2*pi*X/1.) * Foam::sin(2*pi*Y/6.)
	* ((Foam::cos(2*pi*t/T_)*Foam::cos(2*pi*t/T_)) - (Foam::sin(2*pi*t/T_)*Foam::sin(2*pi*t/T_)));
      wDot_[p][2] = 0; 
      */

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
   //scalar t = mesh_.time().value();
    scalar pi = Foam::constant::mathematical::pi, pi2 = pi*pi, T2=T_*T_;
    forAll(motMap_, p)
    {
      //Info << "correcting ALE Motion Mapping\n";
      scalar X = mesh_.points()[p][0];
      scalar Y = mesh_.points()[p][1];
      //scalar Z = mesh_.points()[p][2];
      /*
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
      */
      FatalErrorIn("aleModel.C") << "The coding for the selected motion has to be updated." << abort(FatalError);
      wDot_[p][0] = (X - motMap_[p][0]) * (pi2/T2);
      wDot_[p][1] = (Y - motMap_[p][1]) * (4*pi2/T2);
      wDot_[p][2] = 0; 
      
    }   
  } else if (fictitiousMotionType() == "sinusoid") {
   //scalar t = mesh_.time().value();
    scalar pi = Foam::constant::mathematical::pi, pi2 = pi*pi, T2=T_*T_;
    forAll(motMap_, p)
    {
      //Info << "correcting ALE Motion Mapping\n";
      scalar X = mesh_.points()[p][0];
      scalar Y = mesh_.points()[p][1];
      //scalar Z = mesh_.points()[p][2];
      /*
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
*/
      FatalErrorIn("aleModel.C") << "The coding for the selected motion has to be updated." << abort(FatalError);
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

void aleModel::printProperties()
{
  Info << nl << nl << "Model: " << name_ << nl;
  Info << "aleModel: " << model_ << nl;
  Info << "fictitious Motion Type: " << fictitiousMotionType() << nl;
  Info << "beta: " << beta_ << nl;
  Info << "T: " << T_ << nl;
  Info << nl << nl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

pointTensorField aleModel::piola (pointTensorField& matF, pointTensorField& matH, pointScalarField& matJ)
{
  pointTensorField matP
  (
     IOobject (name_+"P", mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE), pMesh_,
     dimensionedTensor("matP", dimensionSet(1,-1,-2,0,0,0,0), tensor::zero)
  );

  if (model_ == "neoHookean") { 
    forAll(mesh_.points(), n) {
      matP[n] = mu_.value()*pow(matJ[n], -2.0/3.0)*matF[n]
	- ((mu_.value()/3.0)*pow(matJ[n],(-5.0/3.0))*(matF[n] && matF[n])*matH[n])
	+ ( kappa_.value()*(matJ[n]-1.0) )*matH[n];
    }

  } else if (model_ == "mooneyRivlin") {
    forAll(mesh_.points(),n) {
      matP[n] = mu_.value()*(matF[n] - matH[n]/matJ[n]) + lambda_.value() * (matJ[n]-1.0) * matH[n];
    }
  } else {
    FatalErrorIn("aleModel.C") << "ALE Piola Model is not properly defined." << abort(FatalError);
  }
  
  return matP;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*
pointScalarField aleModel::getMaterialPressure (pointTensorField& matF, pointTensorField& matH, pointScalarField& matJ)
{
  pointScalarField matPres
  (
     IOobject ("matPres", mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::AUTO_WRITE),
     pMesh_,
     dimensionedScalar("matPres", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
  );

  if (model_ == "neoHookean") { 
    forAll(mesh_.points(), n) {
      matPres[n] = kappa_.value()*(matJ[n]-1.0) ;
    }
  } else {
    FatalErrorIn("aleModel.C") << "Material Pressure Model is not properly defined." << abort(FatalError);
  }
  
  return matPres;
}
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
}

// ************************************************************************* //
