/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    eofFoam

Description
    Transient solver for coupled system of equations.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fixedGradientFvPatchFields.H"
#include "zeroIonFluxFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    lduMatrix::debug = 0; // Linear solver residuals display disabled
    
    Info<< "\nCalculating concentration distribution\n" << endl;

int ite=0;

   while (runTime.loop()) 
    {
#       include "readSolverControls.H"
        Info << "Time = " << runTime.timeName() << nl << endl;
        Info << "C1 Residual" << tab << tab << "C2 Residual" << tab << tab << "Phi Residual" << endl;

scalar maxResidual_U = 0.0;

#       include "readPISOControls.H"
#       include "CourantNo.H"

		scalar maxResidual = 0.0;
		int iCorr = 0;

if(ite ==0) {
    solve
    (
         fvm::laplacian(PhiInitial) ==  -alpha * Z1*C1-alpha*Z2*C2
    );

Phi = PhiInitial;
electricField=-fvc::grad(Phi);
}

      do {
      do {
        	maxResidual = 0.0;
	Phi.storePrevIter();
	C1.storePrevIter();
	C2.storePrevIter();

	# include "PhiEqn.H"
	# include "C1Eqn.H"
	# include "C2Eqn.H"   
            
Info << eqnResidualC1 << tab << tab << eqnResidualC2 << tab << tab <<  eqnResidualPhi << endl;
	} while (eqnResidualC1 > convergenceTolerance || eqnResidualC2 > convergenceTolerance || eqnResidualPhi > convergenceTolerance);

// Solving Transient Navier Stokes equations


        // Pressure-velocity SIMPLE corrector loop
        for (int corr = 0; corr < nCorr; corr++)
        {
            // Momentum predictor

            tmp<fvVectorMatrix> UEqn
            (
             fvm::ddt(U)
	+ fvm::div(phi, U)
	+ (((beta * Z1 * C1) + (beta * Z2 * C2))* (fvc::grad(Phi)))
     - fvm::laplacian(nu, U)
              
            );

             UEqn().relax();

            eqnResidualU = solve
    (
        UEqn() == -fvc::grad(p)
    ).initialResidual();

    maxResidual_U = max(eqnResidualU, maxResidual_U);

            p.boundaryField().updateCoeffs();
            volScalarField rUA = 1.0/UEqn().A();
            U = rUA*UEqn().H();
            UEqn.clear();
            phi = fvc::interpolate(U) & mesh.Sf();
            adjustPhi(phi, U, p);

            // Store pressure for under-relaxation
            p.storePrevIter();

            // Non-orthogonal pressure corrector loop
            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rUA, p) == fvc::div(phi)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                if (nonOrth == nNonOrthCorr)
                {
                    phi -= pEqn.flux();
                }
            }

#           include "continuityErrs.H"

            // Explicitly relax pressure for momentum corrector
            p.relax();

            // Momentum corrector
            U -= rUA*fvc::grad(p);
            U.correctBoundaryConditions();
Info << eqnResidualU << endl;
        }
}while (eqnResidualU > convergenceTolerance2);
         //turbulence->correct();

Flux1=-D1*fvc::grad(C1)-D1*Z1*C1*fvc::grad(Phi) + (C1*U) ;
Flux1_diff = -D1*fvc::grad(C1);
Flux1_em = -D1*Z1*C1*fvc::grad(Phi);
Flux1_conv = (C1*U) ;
Flux2= -D2*fvc::grad(C2)-D2*Z2*C2*fvc::grad(Phi) + (C2*U) ;
Flux2_diff =  -D2*fvc::grad(C2);
Flux2_em = -D2*Z2*C2*fvc::grad(Phi);
Flux2_conv = (C2*U) ;
rho_e = ((Z1*C1)+(Z2*C2))*96485.3;
rho_E = (((Z1*C1)+(Z2*C2))*96485.3)* electricField;
I_1 = (Z1*Flux1)*96485.3;
I_2 = (Z2*Flux2)*96485.3;
I = I_1 + I_2 ;
I_1_diff = (Z1*Flux1_diff)*96485.3;
I_2_diff = (Z2*Flux2_diff)*96485.3;
I_1_em = (Z1*Flux1_em)*96485.3;
I_2_em = (Z2*Flux2_em)*96485.3;
I_1_conv = (Z1*Flux1_conv)*96485.3;
I_2_conv = (Z2*Flux2_conv)*96485.3;

			Info << eqnResidualC1 << tab << tab << eqnResidualC2 << tab << tab << eqnResidualPhi << tab << eqnResidualU << endl;


		runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

//#       include "convergenceCheck.H"  
ite++;

    }



    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
