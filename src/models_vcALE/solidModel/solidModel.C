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

    P_( IOobject("P", F.time().timeName(), F.db(), IOobject::NO_READ, IOobject::NO_WRITE), F.mesh(),
        dimensionedTensor("P", dimensionSet(1,-1,-2,0,0,0,0), tensor::zero)
    ),

    p_( IOobject("p", F.time().timeName(), F.db(), IOobject::NO_READ, IOobject::NO_WRITE), F.mesh(),
        dimensionedScalar("p", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
    ),

    model_(dict.subDict("solidModel").lookup("solidModel")),

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
  
    Ys_( IOobject("Ys", F.time().timeName(), F.db(), IOobject::NO_READ, IOobject::NO_WRITE), F.mesh(), Ys0_),

    strain_p_( IOobject("strain_p", F.time().timeName(), F.db(), IOobject::NO_READ, IOobject::NO_WRITE), F.mesh(),
      dimensionedScalar("strain_p", dimless, 0.0)
    ),
    
    CpInv_ (F),
    tau_ (F),

    vMises_( IOobject("vMises", F.time().timeName(), F.db(), IOobject::NO_READ, IOobject::NO_WRITE), F.mesh(),
     dimensionedScalar("vMises", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
    ),

    HardeningLaw_("constant"),
    delta_(0.0),
    tauInf_(0.0),
    tolerance_(GREAT)
  
{
  if (model_ == "vonMises"){
    Hm_    = dict.subDict("solidModel").subDict("vonMisesDict").lookup("Hm");
    Ys0_   = dict.subDict("solidModel").subDict("vonMisesDict").lookup("Ys");
    dict.subDict("solidModel").subDict("vonMisesDict").lookup("HardeningLaw") >> HardeningLaw_;
    Info << "Hardening Law: " << HardeningLaw_ << nl;    
    if (HardeningLaw_ == "nonLinear") {
      dict.subDict("solidModel").subDict("vonMisesDict").lookup("delta") >> delta_;
      dict.subDict("solidModel").subDict("vonMisesDict").lookup("tauInf") >> tauInf_;
      dict.subDict("solidModel").subDict("vonMisesDict").lookup("tolerance") >> tolerance_;
      Info << "-> delta: " << delta_ << nl
	   << "-> tauInf: " << tauInf_ << nl
	   << "-> tolerance: " << tolerance_ <<nl;
    }
  }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //
solidModel::~solidModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidModel::correct(
		pointTensorField& F, 
		pointTensorField& H, 
		pointScalarField& J, 
		pointScalarField& spatJ, 
		pointScalarField& matJ,
		//
		pointTensorField& CpInv,      // inverse of right Cauchy
		pointScalarField& strain_p/*,   // plastic strain
		pointTensorField& P,       // Piola
		pointScalarField& pressure */
		){
// This method does like '::correct' but updates 
// 1. CpInv_
// 2. strain_p_
// 3. Piola
// 4. pressure
// from the 'outside'.
// It means they have to be given as extra arguments
//
// IMPORTANT NOTE: this is a polymorphic alternative to the usual '::method'.
//                 As a consequence, it has to have the EXACT SAME CONSEQUENCES
//                 -> The updated pressure has to be the one from this model library, 
//                    and it will be retrieved by the code using the usual 'model.pressure()'
//                 -> The Piola update is done using its reference.
//                 -> The strain_p update is done using its reference
//                 -> The CpInv_ update is done using its reference
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (model_ == "vonMises") {
      // Reconstruction of trueF J from spatJ
      pointScalarField invMatJ = op.inverseScalar(matJ);
      const pointScalarField& spatJ_ = const_cast<pointScalarField&>(spatJ);
      const pointTensorField& F_ = const_cast<pointTensorField&>(F);
      const pointScalarField& J_ = const_cast<pointScalarField&>(J);
      const pointTensorField Finv_ = inv(F_).ref();
     
      b_ = F_ & CpInv & F_.T(); //

      forAll(mesh_.points(), nodeID){
	scalar newJ = spatJ[nodeID] * invMatJ[nodeID];
	p_[nodeID] = kappa_.value() * ( log(newJ) / newJ );   // /!\ the correct pressure is stored in the lib
	op.eigenStructure(b_[nodeID]);
	const vector& eVal = op.eigenValue();    const tensor& eVec = op.eigenVector();
	
	// Principle Trial Deviatoric Kirchoff Stress Vector
	vector tauDevT = vector(
		2.0*mu_.value()*log(sqrt(eVal.x())) - (2./3.)*mu_.value()*log(J_[nodeID]),
		2.0*mu_.value()*log(sqrt(eVal.y())) - (2./3.)*mu_.value()*log(J_[nodeID]),
		2.0*mu_.value()*log(sqrt(eVal.z())) - (2./3.)*mu_.value()*log(J_[nodeID])
				);
	// Yield Criterion
	vector directionV = vector::zero;    scalar plasticM = 0.0;
	double f = sqrt((3.0/2.0)*(tauDevT&tauDevT)) - (Ys0_.value() + Hm_.value()*strain_p[nodeID]);//
	if (HardeningLaw_ == "nonLinear") { 
		f = sqrt((3.0/2.0)*(tauDevT&tauDevT)) - (Ys0_.value() + Hm_.value()*strain_p[nodeID] + (tauInf_-Ys0_.value())*(1-Foam::exp(-delta_*strain_p[nodeID])) ); 
	}

	vector tauDev = tauDevT;
	vector eStretch = vector::zero;
	if (f > 0.0){
	  directionV = tauDevT/(sqrt(2.0/3.0)*sqrt(tauDevT & tauDevT));
	  plasticM = f/(3.0*mu_.value() + Hm_.value());
	  if (HardeningLaw_ == "nonLinear") {
		  // add passing plasticM
	    plasticM = NewtonRaphson(plasticM, F_[nodeID], b_[nodeID], strain_p[nodeID], f);
	  } 

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
	  strain_p[nodeID] += plasticM;//
	}
	
	// Kirchoff Stress Tensor
	tau_[nodeID] = tensor::zero;
	for (int i=0; i<3; i++) {
          scalar pID = p_[nodeID];
	  tau_[nodeID] += (tauDev[i] + (J_[nodeID]*pID))
			  *(vector(eVec[3*i], eVec[3*i+1], eVec[3*i+2])
			  *vector(eVec[3*i], eVec[3*i+1], eVec[3*i+2]));
	}
	// Update PK1
	P_[nodeID] = tau_[nodeID] & Finv_[nodeID].T();
	// Update von-Mises stresses
	vMises_[nodeID] = sqrt(1.5*(tauDev && tauDev));
	// Update CpInv
	CpInv = Finv_ & b_ & Finv_.T();//
      }

    } else { /* If this is not a Von Mises model, do nothing */
   	FatalErrorIn("solidModel.C") << "error in parameters of solid model" << abort(FatalError);
    }


}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void solidModel::correct(pointTensorField& F, pointTensorField& H, pointScalarField& J, pointScalarField& spatJ, pointScalarField& matJ)
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
      pointScalarField invMatJ = op.inverseScalar(matJ);
      
      forAll(mesh_.points(), nodeID){ ///////////////////////////////////////
          scalar newJ = spatJ[nodeID] * invMatJ[nodeID];
	  //p_[nodeID] = kappa_.value() * (newJ -1.0);
	  p_[nodeID] = kappa_.value()*(J_[nodeID]-1.0);

	  P_[nodeID] =
            mu_.value()*pow(J_[nodeID],(-2.0/3.0))*F_[nodeID]
	    - ((mu_.value()/3.0)*pow(J_[nodeID],(-5.0/3.0))*(F_[nodeID] && F_[nodeID])*H_[nodeID])
	    + p_[nodeID]*H_[nodeID];
      }
    } else if (model_ == "MooneyRivlin") { 
       //Compressible Mooney-Rivlin
       // taken from: An upwind vertex centred Finite Volume solver for Lagrangian solid dynamics
       // Miquel Aguirre 1, Antonio J. Gil ∗, Javier Bonet, Chun Hean Lee Zienkiewicz
      const pointTensorField& F__ = const_cast<pointTensorField&>(F);
      const pointTensorField& H__ = const_cast<pointTensorField&>(H);
      const pointScalarField& J__ = const_cast<pointScalarField&>(J);
      const tensorField& H_ = H__.internalField();
      const tensorField& F_ = F__.internalField();
      const scalarField& J_ = J__.internalField();

      forAll(mesh_.points(), nodeID){
	  p_[nodeID] = kappa_.value() * (J_[nodeID]-1.0);
          P_[nodeID] = 
              (mu_.value() * (F_[nodeID] - (H_[nodeID]/J_[nodeID]))) + (lambda_.value()*(J_[nodeID]-1.0)*H_[nodeID]);
      }

    } else if (model_ == "vonMises") {
      // Recconstruction of trueF J from spatJ
      pointScalarField invMatJ = op.inverseScalar(matJ);
      const pointScalarField& spatJ_ = const_cast<pointScalarField&>(spatJ);
      //const pointScalarField  spatJ_clone = spatJ_;
      //const pointScalarField newJ_   = spatJ_clone * invMatJ;

      const pointTensorField& F_ = const_cast<pointTensorField&>(F);
      const pointScalarField& J_ = const_cast<pointScalarField&>(J);


      const pointTensorField Finv_ = inv(F_).ref();
      //p_ = kappa_*(Foam::log(J_)/J_);
     
      b_ = F_ & CpInv_ & F_.T();

      forAll(mesh_.points(), nodeID){
	scalar newJ = spatJ[nodeID] * invMatJ[nodeID];
	//p_[nodeID] = kappa_.value()*( log(J_[nodeID]) )/J_[nodeID];     // pressure from det(F)
	p_[nodeID] = kappa_.value() * ( log(newJ) / newJ );               // pressure from solved spatJ
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
          scalar pID = p_[nodeID];
	  //scalar pID = kappa_.value()*( log(J_[nodeID]) )/J_[nodeID];//p_[nodeID];   // pressure from det
	  //scalar pID = kappa_.value()*( log(newJ) )/newJ;                              // pressure from solved spatJ
	  tau_[nodeID] += (tauDev[i] + (J_[nodeID]*pID))
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

scalar solidModel::NewtonRaphson(scalar plasticM_, tensor Fnode, tensor bnode, scalar strain_pnode, scalar f) {
  scalar plasticM = plasticM_, plasticMold = plasticM_;
  //scalar J = det(Fnode);
  scalar rhs = GREAT;
  // parameters to define: delta_ , tauInf_
  scalar muBar = mu_.value();   //(1./3.)*mu_.value()*Foam::pow(J,-2./3.) * tr(bnode);
  scalar F = f + Ys0_.value() + Hm_.value()*strain_pnode + (tauInf_-Ys0_.value())*(1-Foam::exp(-delta_*strain_pnode));
  while (mag(rhs) > tolerance_){
    scalar T = (3*muBar*plasticMold) + Ys0_.value() + Hm_.value()*(strain_pnode+plasticMold) + (tauInf_-Ys0_.value())*(1-Foam::exp(-delta_*(strain_pnode+plasticMold)));
    rhs = (T-F)/(3*muBar + Hm_.value() + delta_*(tauInf_-Ys0_.value())*Foam::exp(-delta_*(strain_pnode+plasticMold)) );
    plasticMold = plasticM;
    plasticM = plasticMold - rhs;
    //Info << "plasticM: " << plasticM << nl;
  }
  return plasticM;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void solidModel::printProperties() {
    Info<< "\nPrinting material properties ..." << nl
        << "Constitutive model = " << model_ << nl
        << "Density = " << rho_.value() << " " << rho_.dimensions() << nl
        << "Young's modulus = " << E_.value() << " " << E_.dimensions() << nl
        << "Poisson's ratio = " << nu_.value() << " " << nu_.dimensions() << nl
        << "Lame's first parameter lambda = " << lambda_.value() << " " << lambda_.dimensions() << nl
        << "Lame's second parameter mu = " << mu_.value() << " " << mu_.dimensions() << nl
        << "Bulk modulus kappa = " << kappa_.value() << " " << kappa_.dimensions() << nl
	<< "Hardening modulus = " << Hm_.value() << " " << Hm_.dimensions() << nl
	<< "Initial yield stress = " << Ys0_.value() << " " << Ys0_.dimensions() << nl
        << "Linear pressure wave speed = " << Up_.value() << " " << Up_.dimensions() << nl
        << "Linear shear wave speed = " << Us_.value() << " " << Us_.dimensions() << endl;
}

void solidModel::writeVars() {
  if (model_ == "vonMises"){
    //CpInv_.write();
    strain_p_.write();
    vMises_.write();
  }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
