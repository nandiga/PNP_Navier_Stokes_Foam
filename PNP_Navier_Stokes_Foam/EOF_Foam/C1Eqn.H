{
//	C1.storePrevIter();

fvScalarMatrix C1Eqn
    (
        fvm::ddt(C1)
      
     == fvm::laplacian(D1,C1) - fvm::div(Z1*D1*dphi, C1) - fvm::div(phi,C1)
    );

   //C1Eqn.relax();
 
    eqnResidualC1 = C1Eqn.solve().initialResidual();
C1.relax();
    maxResidual = max(eqnResidualC1, maxResidual);

}
