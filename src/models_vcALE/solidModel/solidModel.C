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

    inverseC_(
        IOobject(
	    "inverseC",
	    F.time().timeName(),
	    F.db(),
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	F.mesh(),
	dimensionedTensor("inverseC", dimensionSet(1,-1,-2,0,0,0,0), tensor::zero) // to check
    ),

    epsilon_( dict.lookup("epsilon") ),
    
    one_(
        IOobject
        (
            "one",
            F.time().timeName(),
            F.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        F.mesh(),
        1.0
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

    energyAlgorithm_(
        IOobject
        (
            "energyAlgorithm",
            F.time().timeName(),
            F.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        F.mesh(),
        dimensionedScalar
        (
            "energyAlgorithm",
            dimensionSet(1,-1,-2,0,0,0,0),
            0.0
        )
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

    op(mesh_)
{
    p_.write();
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

    const pointTensorField& F__ = const_cast<pointTensorField&>(F);
    const pointTensorField& H__ = const_cast<pointTensorField&>(H);
    const pointScalarField& J__ = const_cast<pointScalarField&>(J);
    const tensorField& H_ = H__.internalField();
    const tensorField& F_ = F__.internalField();
    const scalarField& J_ = J__.internalField();

    if (model_ == "neoHookean"){
      forAll(mesh_.points(), nodeID){
	  p_[nodeID] = kappa_.value()*(J_[nodeID]-1.0);

	  P_[nodeID] =
            mu_.value()*pow(J_[nodeID],(-2.0/3.0))*F_[nodeID]
	    - ((mu_.value()/3.0)*pow(J_[nodeID],(-5.0/3.0))*(F_[nodeID] && F_[nodeID])*H_[nodeID])
	    + p_[nodeID]*H_[nodeID];
      }
    } else if (model_ == "vonMises") {
      scalar deltaGamma = 0;                // step 2a
      vector nu_np1 = vector::zero;       // step 2b
      pointScalarField J_nplus1 = det(F__); // step 3
      pointTensorField bTrial_enp1 =  (F__ & inverseC()) & F__.T() ; //step 5


      pointTensorField tauPrimeTrial(
		      IOobject("tauPrimeTrial", F.time().timeName(), F.db(), IOobject::NO_READ, IOobject::NO_WRITE), 
		      F.mesh(),
		      dimensionedTensor("tauPrimeTrial", dimensionSet(0,0,0,0,0,0,0), tensor::zero)); // step 8
      pointTensorField tauPrime(
		      IOobject("tauPrime", F.time().timeName(), F.db(), IOobject::NO_READ, IOobject::NO_WRITE), 
		      F.mesh(), 
		      dimensionedTensor("tauPrime", dimensionSet(0,0,0,0,0,0,0), tensor::zero)); // step 12
      pointTensorField tau(
		      IOobject("tau", F.time().timeName(), F.db(), IOobject::NO_READ, IOobject::NO_WRITE), 
		      F.mesh(),
		      dimensionedTensor("tau", dimensionSet(0,0,0,0,0,0,0), tensor::zero)); // step 13
      pointTensorField F_m1T = op.invT(F__);// step 14
      pointTensorField invF  = op.inverse(F__);// step 16
      pointTensorField b_enp1(
		      IOobject("b_enp1", F.time().timeName(), F.db(), IOobject::NO_READ, IOobject::NO_WRITE), 
		      F.mesh(),
		      dimensionedTensor("b_enp1", dimensionSet(0,0,0,0,0,0,0), tensor::zero)); // step 15
      forAll(mesh_.points(), nodeID){
	p_[nodeID] = kappa_.value() * Foam::log(J_nplus1[nodeID]) / J_nplus1[nodeID]; // step 4
	
	op.eigenStructure( bTrial_enp1[nodeID] );    // step 6a
	vector eigVal_ = op.eigenValue();            // step 6b
	tensor eigVec_ = op.eigenVector();           // step 6c

        scalar lambdaTrial_e1 = Foam::sqrt(eigVal_[0]);
	scalar lambdaTrial_e2 = Foam::sqrt(eigVal_[1]);
	scalar lambdaTrial_e3 = Foam::sqrt(eigVal_[2]);
	
	vector n1_np1 = vector(eigVec_.xx(), eigVec_.yx(), eigVec_.zx());//eigVec_.col(0);
	vector n2_np1 = vector(eigVec_.xy(), eigVec_.yy(), eigVec_.zy());
	vector n3_np1 = vector(eigVec_.xz(), eigVec_.yz(), eigVec_.zz());  // step 7

	tauPrimeTrial[nodeID].xx() =
	  (2*mu_.value()*Foam::log(lambdaTrial_e1)) - ((mu_.value()*2./3.) * Foam::log(J_nplus1[nodeID]))  ; // step 8a
	tauPrimeTrial[nodeID].yy() =
	  (2*mu_.value()*Foam::log(lambdaTrial_e2)) - ((mu_.value()*2./3.) * Foam::log(J_nplus1[nodeID]))  ; // step 8b
	tauPrimeTrial[nodeID].zz() =
	  (2*mu_.value()*Foam::log(lambdaTrial_e3)) - ((mu_.value()*2./3.) * Foam::log(J_nplus1[nodeID]))  ; // step 8c

	// IMPORTANT
	// H will have to be given as a parameter
	// tau0 will have to be computed
	// yield Criterion may be reshaped as a function.
	scalar H = 0.1e09 ; // hardening parameter
	scalar tau0 = 0.4e09; // initial yield stress
	scalar yieldCriterion = Foam::sqrt(1.5*(tauPrimeTrial[nodeID]&&tauPrimeTrial[nodeID])) - (tau0 + H*epsilon().value());

	if (yieldCriterion > 0) {
	  scalar normTauPrimeTrialLocal = Foam::mag(tauPrimeTrial[nodeID]);
	  nu_np1[0] = tauPrimeTrial[nodeID].xx()/(Foam::sqrt(2./3.)*normTauPrimeTrialLocal); // step 9a
	  nu_np1[1] = tauPrimeTrial[nodeID].yy()/(Foam::sqrt(2./3.)*normTauPrimeTrialLocal); // step 9b
	  nu_np1[2] = tauPrimeTrial[nodeID].zz()/(Foam::sqrt(2./3.)*normTauPrimeTrialLocal); // step 9c

	  deltaGamma = yieldCriterion/(3*mu_.value() + H); // step 10

	  scalar lambdanp1_e1 = Foam::exp(Foam::log(lambdaTrial_e1)-(deltaGamma*nu_np1[0])); // step 11a
	  scalar lambdanp1_e2 = Foam::exp(Foam::log(lambdaTrial_e2)-(deltaGamma*nu_np1[1])); // step 11a
	  scalar lambdanp1_e3 = Foam::exp(Foam::log(lambdaTrial_e3)-(deltaGamma*nu_np1[2])); // step 11a

	  tauPrime[nodeID].xx() = (1 - (2*mu_.value()*deltaGamma/(Foam::sqrt(2./3.) * normTauPrimeTrialLocal))) * tauPrimeTrial[nodeID].xx();// step 12a
	  tauPrime[nodeID].yy() = (1 - (2*mu_.value()*deltaGamma/(Foam::sqrt(2./3.) * normTauPrimeTrialLocal))) * tauPrimeTrial[nodeID].yy();// step 12b
	  tauPrime[nodeID].zz() = (1 - (2*mu_.value()*deltaGamma/(Foam::sqrt(2./3.) * normTauPrimeTrialLocal))) * tauPrimeTrial[nodeID].zz();// step 12c	  

	  // Here is used J, the argument passed during the call.
	  // However it should be the same than J_np1
	  scalar tau11 = tauPrime[nodeID].xx() + (J[nodeID]*p_[nodeID]) ;  // step 13a
	  scalar tau22 = tauPrime[nodeID].yy() + (J[nodeID]*p_[nodeID]) ;  // step 13b
	  scalar tau33 = tauPrime[nodeID].zz() + (J[nodeID]*p_[nodeID]) ;  // step 13c
	  tau[nodeID] = (tau11 * (n1_np1*n1_np1)) + (tau22 * (n2_np1*n2_np1)) + (tau33 * (n3_np1*n3_np1));
	  P_[nodeID] = tau[nodeID] & F_m1T[nodeID]; // step 14
	  b_enp1[nodeID] = (lambdanp1_e1*(n1_np1*n1_np1)) + (lambdanp1_e2*(n2_np1*n2_np1)) + (lambdanp1_e3*(n3_np1*n3_np1)); // step 15
	  inverseC_[nodeID] = (invF[nodeID] & b_enp1[nodeID]) & F_m1T[nodeID]; // step 16a
	  epsilon_ += deltaGamma; // step 16b
	}
	// Note: from step 14 (at most): possibility to do operations out of loop
	// to gain visibility
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
