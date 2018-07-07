Info<< "Reading field C1 \n" << endl;

    volScalarField C1
    (
        IOobject
        (
            "C1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

Info<< "Reading field C2 \n" << endl;

    volScalarField C2
    (
        IOobject
        (
            "C2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


Info<< "Reading field Phi \n" << endl; 

    volScalarField Phi 
    (
        IOobject
        (
            "Phi", 
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

Info<< "Reading initial field Phi \n" << endl; 

    volScalarField PhiInitial 
    (
        IOobject
        (
            "PhiInitial", 
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

Info<< "Reading field Flux\n" << endl;
volVectorField electricField
(
        IOobject
        (
            "electricField",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        -fvc::grad(Phi)//mesh
    );

    Info << "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info << "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );



   Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


   /* Info<< "Reading diffusivity D1\n" << endl;*/

    dimensionedScalar D1
    (
        transportProperties.lookup("D1")
    );

/* Info<< "Reading diffusivity D2\n" << endl;*/

    dimensionedScalar D2
    (
        transportProperties.lookup("D2")
    );

/* Info<< "Reading parameter Z1\n" << endl;*/

    dimensionedScalar Z1
    (
        transportProperties.lookup("Z1")
    );

/* Info<< "Reading parameter Z2\n" << endl;*/

    dimensionedScalar Z2
    (
        transportProperties.lookup("Z2")
    );

/* Info<< "Reading parameter K1\n" << endl;*/

    dimensionedScalar K1
    (
        transportProperties.lookup("K1")
    );

/* Info<< "Reading parameter K2\n" << endl;*/

    dimensionedScalar K2
    (
        transportProperties.lookup("K2")
    );

/* Info<< "Reading parameter R\n" << endl;*/

    dimensionedScalar R
    (
        transportProperties.lookup("R")
    );

/* Info<< "Reading parameter F\n" << endl;*/

  dimensionedScalar F
    (
        transportProperties.lookup("F")
    );

/* Info<< "Reading parameter T\n" << endl;*/

  dimensionedScalar T
    (
        transportProperties.lookup("T")
    );

/* Info<< "Reading parameter alpha\n" << endl;*/

    dimensionedScalar alpha
    (
        transportProperties.lookup("alpha")
    );

/* Info<< "Reading parameter beta\n" << endl;*/

    dimensionedScalar beta
    (
        transportProperties.lookup("beta")
    );


/* Info<< "Reading parameter nu\n" << endl;*/

    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );

	scalar eqnResidualC1;
	scalar eqnResidualC2;
	scalar eqnResidualPhi;
        scalar eqnResidualU;

 Info<< "Reading/calculating face flux field phi\n" << endl;
 
surfaceScalarField dphi
 (
    IOobject
    (
         "dphi",
         runTime.timeName(),
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::AUTO_WRITE
     ),
  linearInterpolate(electricField) & mesh.Sf()
 );

Info<< "Reading field space charge \n" << endl;

    volScalarField rho_e
    (
        IOobject
        (
            "rho_e",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        ((Z1*C1)+(Z2*C2))*96485.3
    );


Info<< "Reading body force rho_E \n" << endl;
    volVectorField rho_E
    (
        IOobject
        (
            "rho_E",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        (((Z1*C1)+(Z2*C2))*F.value())* electricField
    );


Info<< "Reading field Flux1 due to C1\n" << endl;
    volVectorField Flux1
    (
        IOobject
        (
            "Flux1",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        -D1*fvc::grad(C1)-D1*Z1*C1*fvc::grad(Phi) + (C1*U)

    );

Info<< "Reading field Flux1_diff due to diffusion \n" << endl;
    volVectorField Flux1_diff
    (
        IOobject
        (
            "Flux1_diff",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        -D1*fvc::grad(C1)
    );

Info<< "Reading field Flux1_em due to electromigration\n" << endl;
    volVectorField Flux1_em
    (
        IOobject
        (
            "Flux1_em",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        -D1*Z1*C1*fvc::grad(Phi) 
    );

Info<< "Reading field Flux1_conv due to convective velocity \n" << endl;
    volVectorField Flux1_conv
    (
        IOobject
        (
            "Flux1_conv",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        (C1*U)
    );

 Info<< "Reading field Flux2 due to C2 \n" << endl;
    volVectorField Flux2
    (
        IOobject
        (
            "Flux2",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       -D2*fvc::grad(C2)-D2*Z2*C2*fvc::grad(Phi) + (C2*U)
    );

 Info<< "Reading field Flux2_diff due to diffusion \n" << endl;
    volVectorField Flux2_diff
    (
        IOobject
        (
            "Flux2_diff",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       -D2*fvc::grad(C2)
    );

Info<< "Reading field Flux2_em due to electric migration \n" << endl;
    volVectorField Flux2_em
    (
        IOobject
        (
            "Flux2_em",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
      -D2*Z2*C2*fvc::grad(Phi)
    );

Info<< "Reading field Flux2_conv due to convective velocity \n" << endl;
    volVectorField Flux2_conv
    (
        IOobject
        (
            "Flux2_conv",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        (C2*U)
    );
Info<< "Reading field current I_1 \n" << endl;
    volVectorField I_1
    (
        IOobject
        (
            "I_1",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       (Z1*Flux1)*F.value()
    );

Info<< "Reading field current I_2 \n" << endl;
    volVectorField I_2
    (
        IOobject
        (
            "I_2",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       (Z2*Flux2)*F.value()
    );

Info<< "Reading field current I \n" << endl;
    volVectorField I
    (
        IOobject
        (
            "I",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       I_1 + I_2 
    );

Info<< "Reading field current I_1_diff \n" << endl;
    volVectorField I_1_diff
    (
        IOobject
        (
            "I_1_diff",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       (Z1*Flux1_diff)*F.value()
    );

Info<< "Reading field current I_2_diff \n" << endl;
    volVectorField I_2_diff
    (
        IOobject
        (
            "I_2_diff",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       (Z2*Flux2_diff)*F.value()
    );

Info<< "Reading field current I_1_em \n" << endl;
    volVectorField I_1_em
    (
        IOobject
        (
            "I_1_em",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       (Z1*Flux1_em)*F.value()
    );

Info<< "Reading field current I_2_em \n" << endl;
    volVectorField I_2_em
    (
        IOobject
        (
            "I_2_em",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       (Z2*Flux2_em)*F.value()
    );


Info<< "Reading field current I_1_conv \n" << endl;
    volVectorField I_1_conv
    (
        IOobject
        (
            "I_1_conv",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       (Z1*Flux1_conv)*F.value()
    );

Info<< "Reading field current I_2_conv \n" << endl;
    volVectorField I_2_conv
    (
        IOobject
        (
            "I_2_conv",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
       (Z2*Flux2_conv)*F.value()
    );

#   include "createPhi.H"


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    //setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);
 setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);


   
