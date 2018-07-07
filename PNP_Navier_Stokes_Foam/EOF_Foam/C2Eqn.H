{
//	C2.storePrevIter();

    fvScalarMatrix C2Eqn
    (
        fvm::ddt(C2)

       == fvm::laplacian(D2,C2)-fvm::div(Z2*D2*dphi, C2) - fvm::div(phi,C2)

    );

//  C2Eqn.relax();
    
    eqnResidualC2 = C2Eqn.solve().initialResidual();
 C2.relax(); 
    maxResidual = max(eqnResidualC2, maxResidual);  


}
