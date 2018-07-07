    // Solve the Momentum equation

    tmp<fvVectorMatrix> UEqn
    (
      fvm::div(phi, U)
	+ (((beta * Z1 * C1) + (beta * Z2 * C2))* (fvc::grad(Phi)))
     - fvm::laplacian(nu, U)
    );

    UEqn().relax();

    eqnResidualU = solve
    (
        UEqn() == -fvc::grad(p)
    ).initialResidual();

    maxResidual = max(eqnResidualU, maxResidual);

