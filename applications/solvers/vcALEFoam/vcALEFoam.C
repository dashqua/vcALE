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

Application
    solidVertexFoam

Description
    A large strain solid mechanics solver based on a linear momentum/strains
    mixed formualtion. An explicit Total Lagrangian formulation utilisiing
    a monolithic Total Variation Diminishing Runge-Kutta time integrator.
    A discrete angular momentum projection algorithm based on two global
    Lagrange Multipliers is added for angular momentum conservation.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pointFields.H"
#include "gradientSchemes.H"
#include "interpolationSchemes.H"
#include "angularMomentum.H"
#include "solidModel.H"
#include "mechanics.H"
#include "dualMesh.H"
#include "operations.H"
#include "aleModel.H"
#include "paperFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readControls.H"
    #include "createFields.H"
    #include "checks.H"
    #include "postPro.H" 
	
    while (runTime.loop())
    {
        if (timeStepping == "variable")
        {	 
	    scalar minStretch = GREAT;
	    forAll(stretch, n) { if (stretch[n] < minStretch) { minStretch = stretch[n]; } }
            deltaT = (cfl*h)/(model.Up()/minStretch);	
            runTime.setDeltaT(deltaT);
        }
        t += deltaT; tstep++;
	
        Info << "\nTime Step =" << tstep << "\n deltaT = " << deltaT.value() << " s"
             << "\n Time = " << t.value() << " s" << endl;

	x.oldTime(); 
	xw.oldTime(); 
        spatJ.oldTime();	
        F.oldTime();        
	matJ.oldTime();
	matF.oldTime();
	E.oldTime(); 
	DH.oldTime();
	alphaH.oldTime();
        lm.oldTime(); 
	solvedW.oldTime();
	
	for(int n=0; n<nbRKstages; n++) {
            #include "gEqns.H"
            //
            #include "updateVariables.H"
            //
	    x       = (RKcoefTable[n][0] * x.oldTime()) + (RKcoefTable[n][1] * x);
	    xw      = (RKcoefTable[n][0] * xw.oldTime()) + (RKcoefTable[n][1] * xw);
	    spatJ   = (RKcoefTable[n][0] * spatJ.oldTime()) + (RKcoefTable[n][1] * spatJ);
	    F       = (RKcoefTable[n][0] * F.oldTime()) + (RKcoefTable[n][1] * F);
	    matJ    = (RKcoefTable[n][0] * matJ.oldTime()) + (RKcoefTable[n][1] * matJ);
	    matF    = (RKcoefTable[n][0] * matF.oldTime()) + (RKcoefTable[n][1] * matF);
	    E       = (RKcoefTable[n][0] * E.oldTime()) + (RKcoefTable[n][1] * E);
	    DH      = (RKcoefTable[n][0] * DH.oldTime()) + (RKcoefTable[n][1] * DH);
	    alphaH  = (RKcoefTable[n][0] * alphaH.oldTime()) + (RKcoefTable[n][1] * alphaH);
	    lm      = (RKcoefTable[n][0] * lm.oldTime()) + (RKcoefTable[n][1] * lm);
	    solvedW = (RKcoefTable[n][0] * solvedW.oldTime()) + (RKcoefTable[n][1] * solvedW);
	}


/*	
	forAll(RKstage, i) {
            #include "gEqns.H"
            if (RKstage[i] == 0) {
                #include "updateVariables.H"
            }
        }
	
	x        = 0.5*(x.oldTime() + x);
	xw       = 0.5*(xw.oldTime() + xw);
	spatJ    = 0.5*(spatJ.oldTime() + spatJ);
	F        = 0.5*(F.oldTime()     + F);
	matJ     = 0.5*(matJ.oldTime()  + matJ);
	matF     = 0.5*(matF.oldTime()  + matF);
	E        = 0.5*(E.oldTime() + E);
	DH       = 0.5*(DH.oldTime() + DH);
	alphaH   = 0.5*(alphaH.oldTime() + alphaH);
        lm       = 0.5*(lm.oldTime()    + lm);   
	solvedW  = 0.5*(solvedW.oldTime() + solvedW);
	
        #include "updateVariables.H"
*/

        if (runTime.outputTime())
        {
            #include "postPro.H"
  	}

	/*
        if (curTimeStep < updateEvery){
          Info << "Updating material motion ..." << nl;
	  curTimeStep += 1;
	} else if (curTimeStep == updateEvery) {
          curTimeStep = 0;
	} else {
          FatalErrorIn("vcALEFoam.C") <<"Error in curTimeStep "<< abort(FatalError);
	}
	*/

        Info << " Simulation completed = "
             << (t.value()/runTime.endTime().value())*100 << "%" 
	     << " - Clock = "
	     << runTime.elapsedClockTime() << "s"
	     << nl;

    }

    Info<< "\nExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
