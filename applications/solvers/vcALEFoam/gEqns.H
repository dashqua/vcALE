pR = op.inverseScalar(matJ) * lm ;

// Compute interior terms for right hand sides
rhsF.primitiveFieldRef() = op.surfaceSum(vC*magSf);
rhsMatF.primitiveFieldRef() = op.surfaceSum(wC*magSf);

rhsLm1.primitiveFieldRef() = op.surfaceSum(tC*magSf);
if (GCLtechnique == "withInvolution") {
    rhsJ.primitiveFieldRef() = aleH && op.surfaceSum(wC*magSf);
    rhsLm2.primitiveFieldRef() = op.surfaceSum(tC2*magSf) & (aleH.T() & w_);
    rhsLm3.primitiveFieldRef() = pR * (aleH && (op.surfaceSum(wC*magSf)));  // use rhsJ to be faster
    rhsLm.primitiveFieldRef() = rhsLm1 + rhsLm2 + rhsLm3; //rhsLm1 + rhsLm2 + rhsLm3;
} else if (GCLtechnique == "withoutInvolution") {
    rhsJ.primitiveFieldRef() = op.surfaceSum(HTwC*magSf);
    rhsLm2.primitiveFieldRef() = pR.primitiveFieldRef() * rhsJ.primitiveFieldRef();
    rhsLm.primitiveFieldRef() = rhsLm1 + rhsLm2;
} else {
    FatalErrorIn("gEqns.H") << "GCL technique is not properly defined." << abort(FatalError);
}

// Compute boundary terms for right hand sides
forAll(mesh.boundary(), patch)
{
    forAll(mesh.boundaryMesh()[patch], facei)
    {
        const label& face = mesh.boundaryMesh()[patch].start() + facei;

        forAll(mesh.faces()[face], nodei)
        {
            const label& node = mesh.faces()[face][nodei];
            label nodeB = -1;
            label nodeC = -1;

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

	    vector vH1 = vHat[node];
	    vector vH2 = vHat[nodeB];
	    vector vH3 = vHat[nodeC];
 
	    vector w1 = w_[node];
	    vector w2 = w_[nodeB];
	    vector w3 = w_[nodeC];

            scalar jC_1(0.0),
	           jC_2(0.0),
   	           jC_3(0.0);
	      
	    vector tC_1(vector::zero), 
		   tC_2(vector::zero), 
	      	   tC_3(vector::zero),
		   vHatClocal(vector::zero),
	           wClocal(vector::zero),
	      rhsLmSum(vector::zero);

	    tensor tC2_1(tensor::zero), 
		   tC2_2(tensor::zero), 
	           tC2_3(tensor::zero),
	           tC3_1(tensor::zero), 
		   tC3_2(tensor::zero), 
	           tC3_3(tensor::zero);
	    
	    pointTensorField HT = aleH.T();

	    	
	    const fvPatch& patx = mesh.boundary()[patch];
	    if (patx.name().find("free") != string::npos) {
	      //tC_1 = vector::zero;
	      //tC_2 = vector::zero;
	      //tC_3 = vector::zero;  //no need to rewrite vector::zero
            } else if (patx.name().find("traction") != string::npos) {
                tC_1 = traction;
                tC_2 = traction;
                tC_3 = traction;           
            } else if (patx.name().find("fixed") != string::npos) {
		tC_1 = PHat[node]  & N_p.boundaryField()[patch][facei];   
	        tC_2 = PHat[nodeB] & N_p.boundaryField()[patch][facei];
		tC_3 = PHat[nodeC] & N_p.boundaryField()[patch][facei];
	    } else {	    
		FatalErrorIn("gEqns.H") << "A boundary is not properly defined for tC" << abort(FatalError);
	    }
	    
      	    vHatClocal   = 0.125*(6.0*vH1 + vH2 + vH3);
	    wClocal      = 0.125*(6.0*w1 + w2 + w3);

            rhsLmSum = 0.125*(6*tC_1 + tC_2 + tC_3);
	    if (GCLtechnique == "withInvolution") {
	        tC2_1 = (pR[node]  * N_p.boundaryField()[patch][facei]);
	        tC2_2 = (pR[nodeB] * N_p.boundaryField()[patch][facei]);
	        tC2_3 = (pR[nodeC] * N_p.boundaryField()[patch][facei]);
	        //
	        tC3_1 = w1 * N_p.boundaryField()[patch][facei];
	        tC3_2 = w2 * N_p.boundaryField()[patch][facei];
	        tC3_3 = w3 * N_p.boundaryField()[patch][facei];
	        //
	        rhsLmSum += 0.125*( (6*tC2_1 + tC2_2 + tC2_3 ) & (HT[node] & w_[node]) );
	        //rhsLmSum += 0.125*( pR[node] * (aleH[node] && (6*tC3_1 + tC3_2 + tC3_3)));  //always
	    } else if (GCLtechnique == "withoutInvolution") {
   	        jC_1 = (HT[node]  & w1) & N_p.boundaryField()[patch][facei];
	        jC_2 = (HT[nodeB] & w2) & N_p.boundaryField()[patch][facei];
	        jC_3 = (HT[nodeC] & w3) & N_p.boundaryField()[patch][facei];
                rhsLmSum += 0.125*( pR[node] * (6*jC_1 + jC_2 + jC_3));
	    } else {
	        FatalErrorIn("gEqns.H") << "GCL technique is not properly defined." << abort(FatalError);
	    }
	    rhsLm[node] += rhsLmSum * (magSf_p.boundaryField()[patch][facei]/3.0);
            rhsF[node] +=  vHatClocal*(Sf_p.boundaryField()[patch][facei]/3.0);
	    rhsMatF[node] += wClocal*(Sf_p.boundaryField()[patch][facei]/3.0);
        }
    }
}

// Compute right hand sides
op.volumeIntegrate(rhsLm, V);
op.volumeIntegrate(rhsF, V);
op.volumeIntegrate(rhsJ, V);
op.volumeIntegrate(rhsMatF, V);

// Compute node coordinates
x += deltaT * vHat;

// same for xw
xw += deltaT * w_;

// Compute deformation gradient tensor
F += deltaT*rhsF;

// Compute material Jacobian
matJ += deltaT*rhsJ;

//Compute material deformation gradient
matF += deltaT*rhsMatF;

// Compute linear momentum
// lm += deltaT*rhsLm;

// pTilde receives rhs, then passes to lm
//pTilde += deltaT*rhsLm;
//pR = op.inverseScalar(matJ) * pTilde;
//lm = aleJ * pR;
pointScalarField coef = matJ ; //aleJ;
coef *= op.inverseScalar(matJ);
coef *= deltaT;
pointVectorField ooo = coef * rhsLm;
lm += ooo ;
lm.correctBoundaryConditions();
  //

 /*
pTilde += deltaT*rhsLm;
lm = (aleJ * ( op.inverseScalar(matJ)*pTilde ) ); 
lm.correctBoundaryConditions();
pTilde = op.inverseScalar(aleJ) * ( matJ*lm );
*/
