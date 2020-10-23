// Compute interior terms for right hand sides
rhsLm.primitiveFieldRef() = op.surfaceSum(tC*magSf);
rhsF.primitiveFieldRef() = op.surfaceSum(vC*Sf);
rhsJ.primitiveFieldRef() = aleH && op.surfaceSum(wC*Sf);
rhsMatF.primitiveFieldRef() =  op.surfaceSum(wC*Sf);

if (angularMomentumConservation == "yes")
{
    xe = inter.pointToEdge(x);
    rhsAm.primitiveFieldRef() = op.surfaceSum((xe ^ tC)*magSf);
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

            if (nodei == 0)
            {
                nodeB = mesh.faces()[face][1];
                nodeC = mesh.faces()[face][2];
            }
            else if (nodei == 1)
            {
                nodeB = mesh.faces()[face][2];
                nodeC = mesh.faces()[face][0];
            }
            else if (nodei == 2)
            {
                nodeB = mesh.faces()[face][0];
                nodeC = mesh.faces()[face][1];
            }

	    vector vH1 = vHat[node];
	    vector vH2 = vHat[nodeB];
	    vector vH3 = vHat[nodeC];
 
	    vector w1 = w_[node];
	    vector w2 = w_[nodeB];
	    vector w3 = w_[nodeC];

            vector x1 = x[node];
            vector x2 = x[nodeB];
            vector x3 = x[nodeC];

	    vector t1(vector::zero), 
		   t2(vector::zero), 
		   t3(vector::zero),
		   vHatClocal(vector::zero),
		   wClocal(vector::zero),
		   tC(vector::zero);
	    	
	    const fvPatch& patx = mesh.boundary()[patch];
	    if (patx.name().find("free") != string::npos)
            {
                t1 = vector::zero;
                t2 = vector::zero;
                t3 = vector::zero;
            }
            else if (patx.name().find("traction") != string::npos)
            {
                t1 = traction;
                t2 = traction;
                t3 = traction;
            }
	    else if (patx.name().find("fixed") != string::npos) 
	    {
		t1 = PHat[node]  & N_p.boundaryField()[patch][facei];   
	        t2 = PHat[nodeB] & N_p.boundaryField()[patch][facei];
		t3 = PHat[nodeC] & N_p.boundaryField()[patch][facei];
		// make sure vC is 0
		vHatClocal = vector::zero;
	    }
	    else
	    {	    
		Info << "Error: A boundary is not properly defined for tC" << nl;
		return 0;
	    }
			//
            tC = 0.125*(6.0*t1 + t2 + t3);
      	    vHatClocal = 0.125*(6.0*vH1 + vH2 + vH3);
	    wClocal = 0.125*(6.0*w1 + w2 + w3);

            rhsLm[node] += tC*(magSf_p.boundaryField()[patch][facei]/3.0);
            rhsF[node] +=
                vHatClocal*(Sf_p.boundaryField()[patch][facei]/3.0);
	    rhsMatF[node] +=
		 wClocal*(Sf_p.boundaryField()[patch][facei]/3.0);

            if (angularMomentumConservation == "yes")
            {
                vector xC = 0.125*(6.0*x1 + x2 +x3);
                rhsAm[node] += (xC ^ tC)*(magSf_p.boundaryField()[patch][facei]/3.0);
            }
        }
    }
}

// Compute right hand sides
op.volumeIntegrate(rhsLm, V);
op.volumeIntegrate(rhsF, V);
op.volumeIntegrate(rhsJ, V);
op.volumeIntegrate(rhsMatF, V);

// Conservation of angular momentum
if (angularMomentumConservation == "yes")
{
    op.volumeIntegrate(rhsAm, V);
    am.AMconservation(rhsLm, rhsLm1, rhsAm, i, V);
}


// Compute node coordinates
//v_ = op.inverseScalar(aleJ) * lm/rho;
//vHat = v_ + ((F & op.inverse(matF)) & w_);
x += deltaT * vHat;

// same for xw
xw += deltaT * w_;

// Compute linear momentum
lm += deltaT*rhsLm;
lm.correctBoundaryConditions();

// Compute deformation gradient tensor
F += deltaT*rhsF;

// Compute material Jacobian
matJ += deltaT*rhsJ;

//Compute material deformation gradient
matF += deltaT*rhsMatF;
