// update pR 
//pR = op.inverseScalar(matJ) * lm ;

// Compute interior terms for right hand sides
//rhsF.primitiveFieldRef() = op.surfaceSum(vC*magSf);
//rhsMatF.primitiveFieldRef() = op.surfaceSum(wC*magSf); 
rhsF.primitiveFieldRef()    = op.surfaceSum( inter.pointToEdge(vHat)    * Sf );
rhsMatF.primitiveFieldRef() = op.surfaceSum( inter.pointToEdge(solvedW) * Sf ); 

rhsLm4.primitiveFieldRef() = op.surfaceSum(stabC*magSf);

if (PHatTerm == "withInvolution") {
    rhsLm1.primitiveFieldRef() = 
	  ((aleH && op.surfaceSum(P1C*magSf))*vector_100) + ((aleH && op.surfaceSum(P2C*magSf))*vector_010) + ((aleH && op.surfaceSum(P3C*magSf))*vector_001);   
} else if (PHatTerm == "withoutInvolution") {
    rhsLm1.primitiveFieldRef() = op.surfaceSum(tC*magSf);
}

if (GCLtechnique == "withInvolution") {
    rhsJ.primitiveFieldRef() = aleH && op.surfaceSum( inter.pointToEdge(solvedW) * Sf);
    rhsLm2.primitiveFieldRef() = op.surfaceSum(tC2*magSf) & (aleH.T() & solvedW);
    rhsLm3.primitiveFieldRef() = pR * rhsJ.primitiveField();
} else if (GCLtechnique == "withoutInvolution") {
    rhsJ.primitiveFieldRef() = op.surfaceSum(HTwC*magSf);
    rhsLm2.primitiveFieldRef() = op.surfaceSum(pRwHC*magSf);
    //rhsLm3.primitiveFieldRef() -= rhsLm3.primitiveFieldRef();
}

rhsW.primitiveFieldRef() = op.surfaceSum( inter.pointToEdge(starP) & Sf ) / rho.value();
//op.surfaceSum(PstarHC*magSf) / rho.value();

rhsSpatJ1.primitiveFieldRef() = spatH && op.surfaceSum( inter.pointToEdge( v_ ) * Sf );
rhsSpatJ2.primitiveFieldRef() = op.surfaceSum( inter.pointToEdge( (spatH.T()&trueF)&solvedW ) & Sf );
rhsSpatJ3.primitiveFieldRef() = op.surfaceSum( (stabJacC & (inter.pointToEdge(spatH)&N) ) & (inter.pointToEdge(spatH) & Sf) ) ;

rhsE1.primitiveFieldRef() = op.surfaceSum( inter.pointToEdge( PHat.T()&(lm/rho.value()))      & Sf );
rhsE2.primitiveFieldRef() = op.surfaceSum( inter.pointToEdge( (aleH.T()&solvedW) * E/matJ)    & Sf );
rhsE = rhsE1 + rhsE2;

rhsD.primitiveFieldRef()     = //op.surfaceSum( inter.pointToEdge((D&aleH.T()) & solvedW) * Sf ); 
             ( op.surfaceSum( inter.pointToEdge(D1*(aleH.T()&solvedW)) & Sf ) * vector_100) 
	   + ( op.surfaceSum( inter.pointToEdge(D2*(aleH.T()&solvedW)) & Sf ) * vector_010)
	   + ( op.surfaceSum( inter.pointToEdge(D3*(aleH.T()&solvedW)) & Sf ) * vector_001);
rhsAlpha.primitiveFieldRef() = op.surfaceSum( inter.pointToEdge(alpha*(aleH.T()&solvedW)) & Sf );


// predeclare useful vars for BCs treatment
pointTensorField HT = aleH.T();
pointTensorField spatHT = spatH.T();

// Compute boundary terms for right hand sides
forAll(mesh.boundary(), patch) {
    forAll(mesh.boundaryMesh()[patch], facei) {
        const label& face = mesh.boundaryMesh()[patch].start() + facei;
        vector Npf = N_p.boundaryField()[patch][facei];
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

	    vector pR1    = pR[node],       pR2 = pR[nodeB],       pR3 = pR[nodeC]; 
	    tensor D1     = D[node],         D2 = D[nodeB],         D3 = D[nodeC];
	    scalar alpha1 = alpha[node], alpha2 = alpha[nodeB], alpha3 = alpha[nodeC];

	    vector v1 = v_[node];  vector vH1 = vHat[node];  vector w1 = solvedW[node];
	    vector v2 = v_[nodeB]; vector vH2 = vHat[nodeB]; vector w2 = solvedW[nodeB];
	    vector v3 = v_[nodeC]; vector vH3 = vHat[nodeC]; vector w3 = solvedW[nodeC];

            scalar jC_1(0.0), jC_2(0.0), jC_3(0.0), rhsSpat1Clocal(0.0), HTwClocal(0.0);
	     
	    vector tC_1(vector::zero), tC_2(vector::zero), tC_3(vector::zero),
		   vHatClocal(vector::zero),
	           wClocal(vector::zero),
	           rhsLmSum(vector::zero), rhsLm1Clocal(vector::zero), rhsLm2Clocal(vector::zero),
		   rhsSpat2Clocal(vector::zero),
		   rhsAlphaClocal(vector::zero),
		   rhsWClocal(vector::zero);

	    tensor tC2_1(tensor::zero), tC2_2(tensor::zero), tC2_3(tensor::zero),
	           tC3_1(tensor::zero), tC3_2(tensor::zero), tC3_3(tensor::zero),
		   rhsDClocal(tensor::zero); 

	    vector PxlocA = vector(P[node].xx(),P[node].yx(),P[node].zx()),
	           PylocA = vector(P[node].xy(),P[node].yy(),P[node].zy()),
	           PzlocA = vector(P[node].xz(),P[node].yz(),P[node].zz());
	    vector PxlocB = vector(P[nodeB].xx(),P[nodeB].yx(),P[nodeB].zx()),
	           PylocB = vector(P[nodeB].xy(),P[nodeB].yy(),P[nodeB].zy()),
	           PzlocB = vector(P[nodeB].xz(),P[nodeB].yz(),P[nodeB].zz());
 	    vector PxlocC = vector(P[nodeC].xx(),P[nodeC].yx(),P[nodeC].zx()),
	           PylocC = vector(P[nodeC].xy(),P[nodeC].yy(),P[nodeC].zy()),
	           PzlocC = vector(P[nodeC].xz(),P[nodeC].yz(),P[nodeC].zz());

	    //tensor HT1 = aleH[node].T(), HT2 = aleH[nodeB].T(), HT3 = aleH[nodeC].T();
	    //tensor spatHT1 = spatH[node].T(), spatHT2 = spatH[nodeB].T(), spatHT3 = spatH[nodeC].T();
            // important: all quantities operating on available variables (ex: aleH.T(), spatH.T())
	    //            have to be declared carefully here.

            tensor gvC_1 = v1 * Npf, gvC_2 = v2 * Npf, gvC_3 = v3 * Npf;
/*
	    tensor     HT1 = aleH[node].T(),      HT2 = aleH[nodeB].T(),      HT3 = aleH[nodeC].T();
            tensor spatHT1 = spatH[node].T(), spatHT2 = spatH[nodeB].T(), spatHT3 = spatH[nodeC].T();
*/
	    tensor     HT1 = HT[node],      HT2 = HT[nodeB],      HT3 = HT[nodeC];
            tensor spatHT1 = spatHT[node], spatHT2 = spatHT[nodeB], spatHT3 = spatHT[nodeC];
	    
	    vector HTFw1 = (spatHT1 & trueF[node])  & w1, 
		   HTFw2 = (spatHT2 & trueF[nodeB]) & w2, 
		   HTFw3 = (spatHT3 & trueF[nodeC]) & w3;
	    

	    const fvPatch& patx = mesh.boundary()[patch];
	    if (patx.name().find("free") != string::npos) {                 // OK
	      //tC_1 = vector::zero; //tC_2 = vector::zero; //tC_3 = vector::zero; 
            } else if (patx.name().find("traction") != string::npos) {      // OK
                tC_1 = traction;    tC_2 = traction;    tC_3 = traction;
            } else if (patx.name().find("fixed") != string::npos) {         // OK
	      if (PHatTerm == "withoutInvolution") {
		tC_1 = PHat[node]  & Npf;    tC_2 = PHat[nodeB] & Npf;    tC_3 = PHat[nodeC] & Npf;		
	      } else if (PHatTerm == "withInvolution") {
		tC_1 = ((aleH[node] && (PxlocA*Npf))*v100)  + ((aleH[node] && (PylocA*Npf))*v010)  + ((aleH[node] && (PzlocA*Npf))*v001);
		tC_2 = ((aleH[nodeB] && (PxlocB*Npf))*v100) + ((aleH[nodeB] && (PylocB*Npf))*v010) + ((aleH[nodeB] && (PzlocB*Npf))*v001);
		tC_2 = ((aleH[nodeC] && (PxlocC*Npf))*v100) + ((aleH[nodeC] && (PylocC*Npf))*v010) + ((aleH[nodeC] && (PzlocC*Npf))*v001);		
	      }
            } else if (patx.name().find("roller") != string::npos) { // symmetric
              if (PHatTerm == "withoutInvolution") {     
		// Note at the moment, traction and ramp are null
                tC_1 = ((Npf*Npf) & (PHat[node]  & Npf)) ;
		tC_2 = ((Npf*Npf) & (PHat[nodeB] & Npf)) ;
		tC_3 = ((Npf*Npf) & (PHat[nodeC] & Npf)) ;
		//- ((Npf*Npf) & (aleUp___[nodeC]*lm[nodeC] ));//+ ((I - Npf*Npf) & traction);
		
		// WEAK TREATMENT	
		vH1 = (I-(Npf*Npf)) & vH1;  vH2 = (I-(Npf*Npf)) & vH2;  vH3 = (I-(Npf*Npf)) & vH3;
		pR1 = (I-(Npf*Npf)) & pR1;  pR2 = (I-(Npf*Npf)) & pR2;	pR3 = (I-(Npf*Npf)) & pR3;
		
	      } else if (PHatTerm == "withInvolution") {
                FatalErrorIn("gEqns.H") << "roller: PHatTerm=withInv." << abort(FatalError);
	      }
	    } else if (patx.name().find("skew") != string::npos) {  // skew-symmetric
              if (PHatTerm == "withoutInvolution") {
		// Note at the moment, traction and ramp are null
                tC_1 = ((I-(Npf*Npf)) & (PHat[node]  & Npf)) ;
                tC_2 = ((I-(Npf*Npf)) & (PHat[nodeB] & Npf)) ;
 		tC_3 = ((I-(Npf*Npf)) & (PHat[nodeC] & Npf)) ;
		//+ ((Npf*Npf) & traction);
		
		// WEAK TREATMENT
		vH1 = (Npf*Npf) & vH1;    vH2 = (Npf*Npf) & vH2;    vH3 = (Npf*Npf) & vH3;
                pR1 = (Npf*Npf) & pR1;    pR2 = (Npf*Npf) & pR2;    pR3 = (Npf*Npf) & pR3;

	      } else if (PHatTerm == "withInvolution") {
		FatalErrorIn("gEqns.H") << "skew: PHatTerm=withInv." << abort(FatalError);
	      }
            } else if (patx.name().find("procBoundary") != string::npos) {
                // this case is for parallel simulation : do nothing
		// at the moment, this BC does not work, and produces
		// something very wrong !
	    }	    

      	    vHatClocal     = 0.125*(6*vH1 + vH2 + vH3)/3.;
	    wClocal        = 0.125*(6*w1 + w2 + w3)/3.;
	    HTwClocal      = 0.125*(6*((HT1&w1)&Npf) + ((HT2&w2)&Npf) + ((HT3&w3)&Npf))/3.;// here: future term for rhsMatJ
	    rhsSpat1Clocal = 0.125*( spatH[node] && (6*gvC_1 + gvC_2 + gvC_3))/3.;
	    rhsSpat2Clocal = 0.125*(6*HTFw1 + HTFw2 + HTFw3)/3.;
	    rhsWClocal     = 0.125*(6*(starP[node]&Npf) + (starP[nodeB]&Npf) + (starP[nodeC]&Npf))/3.;// here: future term for rhsW
	    rhsDClocal     = 0.125*( 6*(D1*((HT1&w1)&Npf)) + (D2*((HT2&w2)&Npf)) + (D3*((HT3&w3)&Npf)) )/3.;
	    rhsAlphaClocal = 0.125*( 6*(alpha1*(HT1&w1)) + (alpha2*(HT2&w2)) + (alpha3*(HT3&w3)))/3.;
	    rhsLm1Clocal   = 0.125*(6*tC_1 + tC_2 + tC_3)/3.;

	    if (GCLtechnique == "withInvolution") {
	        tC2_1 = (pR1  * Npf);    tC2_2 = (pR2 * Npf);    tC2_3 = (pR3 * Npf);
	        //tC3_1 = w1 * Npf;        tC3_2 = w2 * Npf;       tC3_3 = w3 * Npf;
		rhsLm2Clocal = 0.125*( (6*tC2_1 + tC2_2 + tC2_3 ) & (aleH[node].T() & solvedW[node]) )/3.;
	    } else if (GCLtechnique == "withoutInvolution") {
   	        //jC_1 = (aleH[node].T() & w1) & Npf;    jC_2 = (aleH[nodeB].T() & w2) & Npf;    jC_3 = (aleH[nodeC].T() & w3) & Npf;
		jC_1 = (HT1 & w1) & Npf;    jC_2 = (HT2 & w2) & Npf;    jC_3 = (HT3 & w3) & Npf;
		rhsLm2Clocal = 0.125*( pR[node] * (6*jC_1 + jC_2 + jC_3))/3.;
		//rhsLm2Clocal = 0.125*( 6*(((pR1*w1)&aleH[node]) & Npf) + (((pR2*w2)&aleH[nodeB]) & Npf) + (((pR3*w3)&aleH[nodeC]) & Npf) )/3.;
	    }
            
	    rhsF[node]      += vHatClocal     * (Sf_p.boundaryField()[patch][facei]);
	    rhsMatF[node]   += wClocal        * (Sf_p.boundaryField()[patch][facei]);
	    rhsJ[node]   += HTwClocal         * (magSf_p.boundaryField()[patch][facei]);  //vector::zero;   // At the moment: w is 0 on boundaries
            rhsSpatJ1[node] += rhsSpat1Clocal * (magSf_p.boundaryField()[patch][facei]);
	    rhsSpatJ2[node] += rhsSpat2Clocal & (Sf_p.boundaryField()[patch][facei]);
	    rhsW[node] = vector::zero;//rhsW[node]      += rhsWClocal     * (magSf_p.boundaryField()[patch][facei]); //rhsW[node]       = vector::zero;  // At the moment: fix BCs to 0 for w	    	   	   
            rhsD[node]      += rhsDClocal     * (magSf_p.boundaryField()[patch][facei]);
	    rhsAlpha[node]  += rhsAlphaClocal & (Sf_p.boundaryField()[patch][facei]);
	    rhsLm1[node]    += rhsLm1Clocal   * (magSf_p.boundaryField()[patch][facei]);
	    rhsLm2[node]    += rhsLm2Clocal   * (magSf_p.boundaryField()[patch][facei]);
        }
    }
}

// Compute right hand sides
/*
op.volumeIntegrateVector5(rhsLm1, rhsLm2, rhsLm3, rhsLm4, rhsW, V);
op.volumeIntegrateTensor2(rhsF, rhsMatF, V);
op.volumeIntegrateScalar2(rhsJ, rhsE, V);
//op.volumeIntegrate(rhsLm1,   V);   //vect
//op.volumeIntegrate(rhsLm2,   V);   //vect
//op.volumeIntegrate(rhsLm3,   V);   //vect
//op.volumeIntegrate(rhsLm4,   V);   //vect
//op.volumeIntegrate(rhsF,    V);    //tensor
//op.volumeIntegrate(rhsJ,    V);    //scalar
//op.volumeIntegrate(rhsMatF, V);    //tensor
//op.volumeIntegrate(rhsW,    V);    //vect
//op.volumeIntegrate(rhsE, V);       //scalar
op.volumeIntegrate(rhsSpatJ1, V);       //scalar
op.volumeIntegrate(rhsSpatJ2, V);       //scalar
op.volumeIntegrate(rhsSpatJ3, V);       //scalar
op.volumeIntegrate(rhsD, V);            //tensor
op.volumeIntegrate(rhsAlpha, V);        //scalar
//op.volumeIntegrate(externalForces, V);  //vector
//externalForces.correctBoundaryConditions();
*/
op.volumeIntegrateGlobal(
    rhsF,      rhsMatF,   rhsD,
    rhsLm1,    rhsLm2,    rhsLm3,    rhsLm4, rhsW,
    rhsSpatJ1, rhsSpatJ2, rhsSpatJ3, rhsJ,   rhsE, rhsAlpha,
    V
);

// Update
forAll(x, n) {
  x[n] += deltaT.value() * vHat[n];                                    // Compute node for spatial motion
  xw[n] += deltaT.value() * solvedW[n];                                // Compute nodes for material motion
  spatJ[n] += deltaT.value() * (rhsSpatJ1[n] + rhsSpatJ2[n] + rhsSpatJ3[n]); // Compute spatial Jacobian
  F[n] += deltaT.value() * rhsF[n];                                    // Compute spatial Deformation Gradient
  matJ[n] += deltaT.value() * rhsJ[n];                                 // Compute material Jacobian
  matF[n] += deltaT.value() * rhsMatF[n];                              // Compute material Deformation Gradient
  E[n] += deltaT.value() * rhsE[n];                                    // Compute energy
  DH[n]     += deltaT.value() * rhsD[n];                               // Compute DHat
  alphaH[n] += deltaT.value() * rhsAlpha[n];                           // Compute alphaHat
  lm[n] += deltaT.value() * (rhsLm1[n] + rhsLm2[n] + rhsLm3[n] + rhsLm4[n] + externalForces[n]); // Compute LM 
}

// update solvedW
if (true) { //if (curTimeStep == updateEvery){
  if (usePstar) {
    solvedW += deltaT * rhsW;
  } else {
    solvedW += deltaT * materialModel_.wDot();
  } 
} 
//solvedW.correctBoundaryConditions();

/*
x += deltaT * vHat;                                    // Compute node for spatial motion
xw += deltaT * solvedW;                                // Compute nodes for material motion
spatJ += deltaT * (rhsSpatJ1 + rhsSpatJ2 + rhsSpatJ3); // Compute spatial Jacobian
F += deltaT * rhsF;                                    // Compute spatial Deformation Gradient
matJ += deltaT * rhsJ;                                 // Compute material Jacobian
matF += deltaT * rhsMatF;                              // Compute material Deformation Gradient
E += deltaT * rhsE;                                    // Compute energy
DH     += deltaT * rhsD;                               // Compute DHat
alphaH += deltaT * rhsAlpha;                           // Compute alphaHat
lm += deltaT * (rhsLm1 + rhsLm2 + rhsLm3 + rhsLm4 + externalForces); // Compute LM 
*/
lm.correctBoundaryConditions();

// Strong Roller/SKEW-SYM BC
forAll(mesh.boundary(), patch) {
  const fvPatch& patx = mesh.boundary()[patch];
  if (patx.name().find("roller") != string::npos) {
    forAll(mesh.boundaryMesh()[patch], facei) {
      const label& face = mesh.boundaryMesh()[patch].start() + facei;
      vector Npf = N_p.boundaryField()[patch][facei];
      forAll(mesh.faces()[face], nodei) {
        const label& node = mesh.faces()[face][nodei];
        lm[node]      = (I - (Npf*Npf)) & lm[node];
	//solvedW[node] = (I - (Npf*Npf)) & solvedW[node];
      }
    }
  } else if (patx.name().find("skew") != string::npos) {
    forAll(mesh.boundaryMesh()[patch], facei) {
      const label& face = mesh.boundaryMesh()[patch].start() + facei;
      vector Npf = N_p.boundaryField()[patch][facei];
      forAll(mesh.faces()[face], nodei) {
        const label& node = mesh.faces()[face][nodei];
        lm[node]      = (Npf*Npf) & lm[node];
	//solvedW[node] =  (I - (Npf*Npf)) & solvedW[node];
      }
    }
    
  } else { }
}

