{
//	Phi.storePrevIter();
    fvScalarMatrix PhiEqn
    (
        fvm::laplacian(Phi) ==  -alpha * Z1*C1-alpha*Z2*C2 
    );

//   PhiEqn.relax();
    
    eqnResidualPhi = PhiEqn.solve().initialResidual(); 
 Phi.relax(); 
    maxResidual = max(eqnResidualPhi, maxResidual);  
electricField=-fvc::grad(Phi);
dphi=linearInterpolate(electricField) & mesh.Sf();

}
