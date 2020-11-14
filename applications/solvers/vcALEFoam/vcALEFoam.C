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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readControls.H"
    #include "createFields.H"
    #include "postPro.H" 
	
    while (runTime.loop())
    {
	//    Info << "before tstep.............." << nl;
        if (timeStepping == "variable")
        {
            deltaT = (cfl*h)/model.Up();	
            runTime.setDeltaT(deltaT);
        }
	//Info << "after tstep" << nl;
        t += deltaT; tstep++;
	
        Info << "\nTime Step =" << tstep << "\n deltaT = " << deltaT.value() << " s"
             << "\n Time = " << t.value() << " s" << endl;
	     
	//Info << "oldTimes" << nl;
	
        F.oldTime();        
	matJ.oldTime();
	matF.oldTime();
        lm.oldTime(); pR.oldTime(); pTilde.oldTime();
	solvedW.oldTime();
	x.oldTime(); 
	xw.oldTime(); 
	E.oldTime(); 

	forAll(RKstage, i) {
		//Info << "gEqns.H" << nl;
            #include "gEqns.H"
            if (RKstage[i] == 0) {
		    //Info << "updateVariables.H" << nl;
                #include "updateVariables.H"
            }
        }

	F        = 0.5*(F.oldTime()    + F);
	matJ     = 0.5*(matJ.oldTime() + matJ);
	matF     = 0.5*(matF.oldTime() + matF);
        lm       = 0.5*(lm.oldTime()   + lm);
	solvedW  = 0.5*(solvedW.oldTime() + solvedW);
	x        = 0.5*(x.oldTime() + x);
	xw       = 0.5*(xw.oldTime() + xw);
	E        = 0.5*(E.oldTime() + E);
	
        #include "updateVariables.H"

        if (runTime.outputTime())
        {
            #include "postPro.H"
  	}

	//Info << "tvalue: " << t.value() << nl;
	//Info << "endTime: " << runTime.endTime().value() << nl;
	//Info << "testent:" << (t.value()/runTime.endTime().value())*100 << nl;
        Info << " Simulation completed = "
             << (t.value()/runTime.endTime().value())*100 << "%" << endl;

	//Info << "end tstep..................." << nl;
    }

    Info<< "\nExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
