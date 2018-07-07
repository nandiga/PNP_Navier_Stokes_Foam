// initialize values for convergence checks

    scalar eqnResidual = 1;
    scalar convergenceCriterion = 0;

    simple.readIfPresent("convergence", convergenceCriterion);

