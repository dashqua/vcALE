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

#include "solidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(solidModel, 0);


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

solidModel::solidModel
(
    const pointTensorField& F,
    const dictionary& dict,
    const fvMesh& vm
)
:
    mesh_(vm),

    P_(
        IOobject(
            "P",
            F.time().timeName(),
            F.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        F.mesh(),
        dimensionedTensor("P", dimensionSet(1,-1,-2,0,0,0,0), tensor::zero)
    ),

    p_(
        IOobject
        (
            "p",
            F.time().timeName(),
            F.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        F.mesh(),
        dimensionedScalar("p", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
    ),

    model_(dict.lookup("solidModel")),

    rho_(dict.lookup("rho")),
    E_(dict.lookup("E")),
    nu_(dict.lookup("nu")),
    mu_(E_/(2.0*(1.0 + nu_))),
    lambda_(nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))),
    kappa_(lambda_ + (2.0/3.0)*mu_),

    Up_ (sqrt((lambda_+2.0*mu_)/rho_)),
    Us_ (sqrt(mu_/rho_)),

    op(mesh_),

    // Von Mises vars
    b_ (F & F.T()),
  
    //Hm_(dimensionedScalar("Hm", dimensionSet(1,-1,-2,0,0,0,0), 0.0)),
    //Ys0_(dimensionedScalar("Ys", dimensionSet(1,-1,-2,0,0,0,0), 0.0)),

    Ys_(
      IOobject(
        "Ys",
	F.time().timeName(),
	F.db(),
	IOobject::NO_READ,
	IOobject::NO_WRITE
	),
      F.mesh(),
      Ys0_
    ),

    strain_p_(
      IOobject(
        "strain_p",
	F.time().timeName(),
	F.db(),
	IOobject::NO_READ,
	IOobject::NO_WRITE
	),
      F.mesh(),
      dimensionedScalar("strain_p", dimless, 0.0)
   ),

   CpInv_ (F),
   tau_ (F),

   vMises_(
     IOobject(
       "vMises",
       F.time().timeName(),
       F.db(),
       IOobject::NO_READ,
       IOobject::NO_WRITE
       ),
     F.mesh(),
     dimensionedScalar("vMises", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
   )
  
{
  if (model_ == "vonMises"){
    Hm_  = dict.subDict("vonMisesDict").lookup("Hm");
    Ys0_ = dict.subDict("vonMisesDict").lookup("Ys");
  }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //
solidModel::~solidModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

  void solidModel::correct(pointTensorField& F, pointTensorField& H, pointScalarField& J)
{
  /*
    F has to be 'trueF' now.
    Because we only need the 'true' tensors here, F is retrieved as an argument from runTime
    and H and J and re-computed from that F (which is, as mentionned, trueF).
 
    EDIT: H obviously needs F to be inversed, which is not possible in vanilla OF
    so H will be passed as well. At this point, it is also worth sending J.

    F, H and J come not const, so a const_cast is necessary.
  */
    //const pointMesh& pMesh_ = P_.mesh();
    //const objectRegistry& db = pMesh_.thisDb();

    if (model_ == "neoHookean"){
      const pointTensorField& F__ = const_cast<pointTensorField&>(F);
      const pointTensorField& H__ = const_cast<pointTensorField&>(H);
      const pointScalarField& J__ = const_cast<pointScalarField&>(J);
      const tensorField& H_ = H__.internalField();
      const tensorField& F_ = F__.internalField();
      const scalarField& J_ = J__.internalField();
      
      forAll(mesh_.points(), nodeID){
	  p_[nodeID] = kappa_.value()*(J_[nodeID]-1.0);

	  P_[nodeID] =
            mu_.value()*pow(J_[nodeID],(-2.0/3.0))*F_[nodeID]
	    - ((mu_.value()/3.0)*pow(J_[nodeID],(-5.0/3.0))*(F_[nodeID] && F_[nodeID])*H_[nodeID])
	    + p_[nodeID]*H_[nodeID];
      }
    } else if (model_ == "vonMises") {
      const pointTensorField& F_ = const_cast<pointTensorField&>(F);
      const pointScalarField& J_ = const_cast<pointScalarField&>(J);
      const pointTensorField Finv_ = inv(F_).ref();
      //p_ = kappa_*(Foam::log(J_)/J_);
      b_ = F_ & CpInv_ & F_.T();

      forAll(mesh_.points(), nodeID){
	p_[nodeID] = kappa_.value()*(log(J_[nodeID]))/J_[nodeID];
	op.eigenStructure(b_[nodeID]);
	const vector& eVal = op.eigenValue();
	const tensor& eVec = op.eigenVector();
	
	// Principle Trial Deviatoric Kirchoff Stress Vector
	vector tauDevT = vector(
		2.0*mu_.value()*log(sqrt(eVal.x())) - (2./3.)*mu_.value()*log(J_[nodeID]),
		2.0*mu_.value()*log(sqrt(eVal.y())) - (2./3.)*mu_.value()*log(J_[nodeID]),
		2.0*mu_.value()*log(sqrt(eVal.z())) - (2./3.)*mu_.value()*log(J_[nodeID])
				);
	// Yield Criterion
	vector directionV = vector::zero;
	scalar plasticM = 0.0;	
	double f = sqrt((3.0/2.0)*(tauDevT&tauDevT)) - (Ys0_.value() + Hm_.value()*strain_p_[nodeID]);

	vector tauDev = tauDevT;
	vector eStretch = vector::zero;
	if (f > 0.0){
	  directionV = tauDevT/(sqrt(2.0/3.0)*sqrt(tauDevT & tauDevT));
	  plasticM = f/(3.0*mu_.value() + Hm_.value());

	  // Elastic Stretch Vector
	  eStretch = vector(
			    exp(log(sqrt(eVal.x())) - plasticM*directionV[0]),
			    exp(log(sqrt(eVal.y())) - plasticM*directionV[1]),
			    exp(log(sqrt(eVal.z())) - plasticM*directionV[2])
			    );
	  // Principle Deviatoric Kirchoff Stress Tensor
	  tauDev = vector(
	     (1.0 - ((2.0*mu_.value()*plasticM) / (sqrt(2./3.)*sqrt(tauDevT&tauDevT)))) * tauDevT[0],
	     (1.0 - ((2.0*mu_.value()*plasticM) / (sqrt(2./3.)*sqrt(tauDevT&tauDevT)))) * tauDevT[1],
	     (1.0 - ((2.0*mu_.value()*plasticM) / (sqrt(2./3.)*sqrt(tauDevT&tauDevT)))) * tauDevT[2]
			  );
	  // Update Left Cauchy Green Strain Tensor
	  b_[nodeID] = tensor::zero;
	  b_[nodeID] += 
	    (eStretch[0]*eStretch[0])*( vector(eVec[3*0], eVec[3*0+1], eVec[3*0+2]) *vector(eVec[3*0], eVec[3*0+1], eVec[3*0+2]) );
	  b_[nodeID] +=
	    (eStretch[1]*eStretch[1])*( vector(eVec[3*1], eVec[3*1+1], eVec[3*1+2]) *vector(eVec[3*1], eVec[3*1+1], eVec[3*1+2]) );
	  b_[nodeID] +=
	    (eStretch[2]*eStretch[2])*( vector(eVec[3*2], eVec[3*2+1], eVec[3*2+2]) *vector(eVec[3*2], eVec[3*2+1], eVec[3*2+2]) );
	  // Update Plastic Strain
	  strain_p_[nodeID] += plasticM;
	}
	
	// Kirchoff Stress Tensor
	tau_[nodeID] = tensor::zero;
	for (int i=0; i<3; i++) {
	  tau_[nodeID] += (tauDev[i] + (J_[nodeID]*p_[nodeID]))
			  *(vector(eVec[3*i], eVec[3*i+1], eVec[3*i+2])
			  *vector(eVec[3*i], eVec[3*i+1], eVec[3*i+2]));
	}
	// Update PK1
	P_[nodeID] = tau_[nodeID] & Finv_[nodeID].T();
	// Update von-Mises stresses
	vMises_[nodeID] = sqrt(1.5*(tauDev && tauDev));
	// Update CpInv
	CpInv_ = Finv_ & b_ & Finv_.T();
      }

    } else {
	FatalErrorIn("solidModel.C") << "Solid Model is not properly defined." << abort(FatalError) ;
    }

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void solidModel::printMaterialProperties()
{
    Info<< "\nPrinting material properties ..." << nl
        << "Constitutive model = " << model_ << nl
        << "Density = " << rho_.value() << " " << rho_.dimensions() << nl
        << "Young's modulus = " << E_.value() << " " << E_.dimensions() << nl
        << "Poisson's ratio = " << nu_.value() << " " << nu_.dimensions() << nl
        << "Lame's first parameter lambda = " << lambda_.value() << " "
        << lambda_.dimensions() << nl
        << "Lame's second parameter mu = " << mu_.value() << " "
        << mu_.dimensions() << nl
        << "Bulk modulus kappa = " << kappa_.value() << " "
        << kappa_.dimensions() << nl
        << "Linear pressure wave speed = " << Up_.value() << " "
        << Up_.dimensions() << nl
        << "Linear shear wave speed = " << Us_.value() << " "
        << Us_.dimensions() << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
