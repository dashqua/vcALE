pR = op.inverseScalar(matJ) * lm ;

// Compute interior terms for right hand sides
rhsF.primitiveFieldRef() = op.surfaceSum(vC*magSf);
rhsMatF.primitiveFieldRef() = op.surfaceSum(wC*magSf);

rhsLm4.primitiveFieldRef() = op.surfaceSum(stabC*magSf);
if (PHatTerm == "withoutInvolution") {
    rhsLm1.primitiveFieldRef() = op.surfaceSum(tC*magSf);
} else if (PHatTerm == "withInvolution") {
  rhsLm1.primitiveFieldRef() = 
      ((aleH && op.surfaceSum(P1C*magSf))*vector_100)
    + ((aleH && op.surfaceSum(P2C*magSf))*vector_010)
    + ((aleH && op.surfaceSum(P3C*magSf))*vector_001);   
}

if (GCLtechnique == "withInvolution") {
    rhsJ.primitiveFieldRef() = aleH && op.surfaceSum(wC*magSf);
    rhsLm2.primitiveFieldRef() = op.surfaceSum(tC2*magSf) & (aleH.T() & solvedW);
    rhsLm3.primitiveFieldRef() = pR * rhsJ.primitiveField();
    //rhsLm.primitiveFieldRef() = rhsLm1 + rhsLm2 + rhsLm3;
} else if (GCLtechnique == "withoutInvolution") {
    rhsJ.primitiveFieldRef() = op.surfaceSum(HTwC*magSf);
    rhsLm2.primitiveFieldRef() = op.surfaceSum(pRwHC*magSf);
    //rhsLm.primitiveFieldRef() = rhsLm1 + rhsLm2;
}

rhsW.primitiveFieldRef() = op.surfaceSum(PstarHC*magSf) / rho.value() ; // materialModel_.wDot().primitiveField()

//rhsE.primitiveFieldRef() = op.surfaceSum(Pmat_pmat*magSf);
rhsE1.primitiveFieldRef() = op.surfaceSum( inter.pointToEdge( PHat.T()&(lm/rho.value()))                 & Sf );
rhsE2.primitiveFieldRef() = op.surfaceSum( inter.pointToEdge( (aleH.T()&solvedW) * E/matJ) & Sf );
rhsE = rhsE1 + rhsE2;

// Compute boundary terms for right hand sides
forAll(mesh.boundary(), patch) {
    forAll(mesh.boundaryMesh()[patch], facei) {
        const label& face = mesh.boundaryMesh()[patch].start() + facei;
        forAll(mesh.faces()[face], nodei) {
            const label& node = mesh.faces()[face][nodei];
            label nodeB = -1, nodeC = -1;
            if (nodei == 0) {
                nodeB = mesh.faces()[face][1];
                nodeC = mesh.faces()[face][2];
            } else if (nodei == 1) {
                nodeB = mesh.faces()[face][2];
                nodeC = mesh.faces()[face][0];
            } else if (nodei == 2) {
                nodeB = mesh.faces()[face][0];
                nodeC = mesh.faces()[face][1];
            }

	    vector v1 = v_[node];  vector vH1 = vHat[node];  vector w1 = solvedW[node];
	    vector v2 = v_[nodeB]; vector vH2 = vHat[nodeB]; vector w2 = solvedW[nodeB];
	    vector v3 = v_[nodeC]; vector vH3 = vHat[nodeC]; vector w3 = solvedW[nodeC];
 
            tensor PH1 = matP[node]  & aleH[node];    tensor PstarH1 = starP[node]  & aleH[node];
	    tensor PH2 = matP[nodeB] & aleH[nodeB];   tensor PstarH2 = starP[nodeB] & aleH[nodeB];
	    tensor PH3 = matP[nodeC] & aleH[nodeC];   tensor PstarH3 = starP[nodeC] & aleH[nodeC];

            scalar jC_1(0.0), jC_2(0.0), jC_3(0.0);
	      
	    vector tC_1(vector::zero), 
		   tC_2(vector::zero), 
	      	   tC_3(vector::zero),
		   vHatClocal(vector::zero),
	           wClocal(vector::zero),
	           rhsLmSum(vector::zero),
		   rhsLm1Clocal(vector::zero),
		   rhsLm2Clocal(vector::zero);

	    tensor tC2_1(tensor::zero), 
		   tC2_2(tensor::zero), 
	           tC2_3(tensor::zero),
	           tC3_1(tensor::zero), 
		   tC3_2(tensor::zero), 
	           tC3_3(tensor::zero), 
	           PmatHClocal(tensor::zero),
	           PstarHClocal(tensor::zero);

	    vector PxlocA = vector(P[node].xx(),P[node].yx(),P[node].zx()),
	           PylocA = vector(P[node].xy(),P[node].yy(),P[node].zy()),
	           PzlocA = vector(P[node].xz(),P[node].yz(),P[node].zz());
	    vector PxlocB = vector(P[nodeB].xx(),P[nodeB].yx(),P[nodeB].zx()),
	           PylocB = vector(P[nodeB].xy(),P[nodeB].yy(),P[nodeB].zy()),
	           PzlocB = vector(P[nodeB].xz(),P[nodeB].yz(),P[nodeB].zz());
 	    vector PxlocC = vector(P[nodeC].xx(),P[nodeC].yx(),P[nodeC].zx()),
	           PylocC = vector(P[nodeC].xy(),P[nodeC].yy(),P[nodeC].zy()),
	           PzlocC = vector(P[nodeC].xz(),P[nodeC].yz(),P[nodeC].zz());

	    vector v100(1,0,0), v010(0,1,0), v001(0,0,1);
	    pointTensorField HT = aleH.T();
	    const fvPatch& patx = mesh.boundary()[patch];
            vector Npf = N_p.boundaryField()[patch][facei];
	    
	    if (patx.name().find("free") != string::npos) {
	      //tC_1 = vector::zero;
	      //tC_2 = vector::zero;
	      //tC_3 = vector::zero;  //no need to rewrite vector::zero
            } else if (patx.name().find("traction") != string::npos) {
                tC_1 = traction;
                tC_2 = traction;
                tC_3 = traction;
            } else if (patx.name().find("fixed") != string::npos) {
	      if (PHatTerm == "withoutInvolution") {
		tC_1 = PHat[node]  & Npf;   
	        tC_2 = PHat[nodeB] & Npf;
		tC_3 = PHat[nodeC] & Npf;		
	      } else if (PHatTerm == "withInvolution") {
		tC_1 = ((aleH[node] && (PxlocA*Npf))*v100)
		     + ((aleH[node] && (PylocA*Npf))*v010)
	             + ((aleH[node] && (PzlocA*Npf))*v001);
		tC_2 = ((aleH[nodeB] && (PxlocB*Npf))*v100)
		     + ((aleH[nodeB] && (PylocB*Npf))*v010)
	             + ((aleH[nodeB] && (PzlocB*Npf))*v001);
		tC_2 = ((aleH[nodeC] && (PxlocC*Npf))*v100)
		     + ((aleH[nodeC] && (PylocC*Npf))*v010)
	             + ((aleH[nodeC] && (PzlocC*Npf))*v001);		
	      }
            } else if (patx.name().find("roller") != string::npos) { // symmetric
              if (PHatTerm == "withoutInvolution") {     
		      // Note at the moment, traction and ramp are null
                tC_1 = (Npf*Npf) & (PHat[node]   & Npf);
		tC_2 = (Npf*Npf) & (PHat[nodeB]  & Npf);
		tC_3 = (Npf*Npf) & (PHat[nodeC]  & Npf);
	      } else if (PHatTerm == "withInvolution") {
                FatalErrorIn("gEqns.H") << "'roller BC' and 'PHatTern' == 'withInvolution' not coded yet." << abort(FatalError);
	      }

/*

              //scalar ramp = 0.0;
	       //    /!\   THIS PART IS NOT VALID
	      tC_1 = ((Npf*Npf) & ((PHat[node]&Npf) - (aleUp___[node]*lm[node])))
                   ;//+ ((tensor::I - (Npf*Npf)) & (ramp*traction));
	      tC_2 = ((Npf*Npf) & ((PHat[nodeB]&Npf) - (aleUp___[nodeB]*lm[nodeB])))
                   ;//+ ((tensor::I - (Npf*Npf)) & (ramp*traction));
	      tC_3 = ((Npf*Npf) & ((PHat[nodeC]&Npf) - (aleUp___[nodeC]*lm[nodeC])))
                   ;//+ ((tensor::I - (Npf*Npf)) & (ramp*traction));
		//above is cancelled ATM

	      tC_1 = vector::zero, tC_2 = vector::zero, tC_3 = vector::zero;

	      lm[node] = ((Npf*Npf) & ((PHat[node]&Npf) - (aleUp___[node]*lm[node]))); // + ((tensor::I - (Npf*Npf)) & (ramp*traction));
*/
            } else if (patx.name().find("procBoundary") != string::npos) {
                // this case is for parallel simulation : do nothing
		// at the moment, this BC does not work, and produces
		// something very wrong !
	    } else {	    
	      FatalErrorIn("gEqns.H") << "Problem with BC " << patx.name() << abort(FatalError);
	    }
	    
      	    vHatClocal   = 0.125*(6*vH1 + vH2 + vH3)/3.;
	    wClocal      = 0.125*(6*w1 + w2 + w3)/3.;
	    PmatHClocal  = 0.125*(6*PH1 + PH2 + PH3)/3.;
	    PstarHClocal = 0.125*(6*PstarH1 + PstarH2 + PstarH3)/3.;
            ///////////rhsLmSum     = 0.125*(6*tC_1 + tC_2 + tC_3)/3.;
	    rhsLm1Clocal =  0.125*(6*tC_1 + tC_2 + tC_3)/3.;
	    
	    if (GCLtechnique == "withInvolution") {
	        tC2_1 = (pR[node]  * Npf);
	        tC2_2 = (pR[nodeB] * Npf);
	        tC2_3 = (pR[nodeC] * Npf);
	        //
	        tC3_1 = w1 * Npf;
	        tC3_2 = w2 * Npf;
	        tC3_3 = w3 * Npf;
	        //
	        //rhsLmSum += 0.125*( (6*tC2_1 + tC2_2 + tC2_3 ) & (HT[node] & solvedW[node]) )/3.;
		// no BC for rhsLm3
		rhsLm2Clocal = 0.125*( (6*tC2_1 + tC2_2 + tC2_3 ) & (HT[node] & solvedW[node]) )/3.;
	    } else if (GCLtechnique == "withoutInvolution") {
   	        jC_1 = (HT[node]  & w1) & Npf;
	        jC_2 = (HT[nodeB] & w2) & Npf;
	        jC_3 = (HT[nodeC] & w3) & Npf;
                //rhsLmSum += 0.125*( pR[node] * (6*jC_1 + jC_2 + jC_3))/3.;
		rhsLm2Clocal = 0.125*( pR[node] * (6*jC_1 + jC_2 + jC_3))/3.;
	    }

	    //rhsLm[node]   += rhsLmSum   * (magSf_p.boundaryField()[patch][facei]);
	    rhsLm1[node]  += rhsLm1Clocal  * (magSf_p.boundaryField()[patch][facei]);
	    rhsLm2[node]  += rhsLm2Clocal  * (magSf_p.boundaryField()[patch][facei]);
            rhsF[node]    += vHatClocal    * (Sf_p.boundaryField()[patch][facei]);
	    rhsMatF[node] += wClocal       * (Sf_p.boundaryField()[patch][facei]);
	    
	    //rhsW[node]    += PmatHClocal & (Sf_p.boundaryField()[patch][facei]);	    
	    //rhsW[node]    += PstarHClocal & (Sf_p.boundaryField()[patch][facei]);
	    
	    rhsW[node] = vector::zero;  // At the moment: fix BCs to 0 for w	    	    
	    
	    //solvedW[node]    -= (solvedW[node] & Npf) * Npf; // removing the normal component
	    // check BC for rhsE
	    // J has no BC


	    /*
	    vector rhsEClocal1 = 0.125*( 6*(PH1.T()&v1) + (PH2.T()&v2) + (PH3.T()&v3) )/3.;
            scalar rhsEClocal2 = 0.125*( 6*(jC_1*E[node]/matJ[node]) + (jC_2*E[nodeB]/matJ[nodeB]) + (jC_3*E[nodeC]/matJ[nodeC]) )/3.;
	    rhsE[node] += (rhsEClocal1 & Sf_p.boundaryField()[patch][facei]) + (rhsEClocal2 * magSf_p.boundaryField()[patch][facei]);
	    */
        }
    }
}

// Compute right hand sides
//op.volumeIntegrate(rhsLm,   V);
op.volumeIntegrate(rhsLm1,   V);
op.volumeIntegrate(rhsLm2,   V);
op.volumeIntegrate(rhsLm3,   V);
op.volumeIntegrate(rhsLm4,   V);
op.volumeIntegrate(rhsF,    V);
op.volumeIntegrate(rhsJ,    V);
op.volumeIntegrate(rhsMatF, V);
op.volumeIntegrate(rhsW,    V);
op.volumeIntegrate(rhsE, V);

// Compute solvedW
solvedW += deltaT * materialModel_.wDot();//rhsW;
//solvedW.correctBoundaryConditions();

// Compute node coordinates
x += deltaT * vHat;

// same for xw
xw += deltaT * solvedW;

// Compute deformation gradient tensor
F += deltaT * rhsF;

// Compute material Jacobian
matJ += deltaT * rhsJ;

//Compute material deformation gradient
matF += deltaT * rhsMatF;

// Compute linear momentum
//trac = rhsLm1;
//trac.correctBoundaryConditions();
lm += deltaT * (rhsLm1 + rhsLm2 + rhsLm3 + rhsLm4);//rhsLm;
lm.correctBoundaryConditions();

// Compute materiel energy
E += deltaT * rhsE;
E.correctBoundaryConditions();
