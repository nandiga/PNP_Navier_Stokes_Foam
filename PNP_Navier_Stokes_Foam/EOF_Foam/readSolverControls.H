const dictionary& solverControls = mesh.solutionDict().subDict("solverControls");

int nCorr1(readInt(solverControls.lookup("nCorrectorsMax")));
scalar convergenceTolerance(readScalar(solverControls.lookup("globalTolerance")));
scalar convergenceTolerance2(readScalar(solverControls.lookup("globalTolerance2")));

