    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
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

    #include "createPhi.H"

    threePhaseMixture threePhaseProperties(U, phi);

    volScalarField& alpha1(threePhaseProperties.alpha1());
    volScalarField& alpha2(threePhaseProperties.alpha2());
    volScalarField& alpha3(threePhaseProperties.alpha3());
    volScalarField& alpha4(threePhaseProperties.alpha4()); //new

    const dimensionedScalar& rho1 = threePhaseProperties.rho1();
    const dimensionedScalar& rho2 = threePhaseProperties.rho2();
    const dimensionedScalar& rho3 = threePhaseProperties.rho3();
    const dimensionedScalar& rho4 = threePhaseProperties.rho4();

    //dimensionedScalar D23(threePhaseProperties.lookup("D23"));
	

// New ------------------------------------------------------------------

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
		
	dimensionedScalar D23
	(
		transportProperties.lookup("D23")
	);
	dimensionedScalar D24
	(
		transportProperties.lookup("D24")
	);
	dimensionedScalar D34
	(
		transportProperties.lookup("D34")
	);
	//dimensionedScalar rho4
	//(
	//	transportProperties.lookup("rho4")
	//);

//	volScalarField alpha4	removed to put into threephasemixture
//	(
//	   IOobject
//	   (
//		"alpha.alpha4",
//		runTime.timeName(),
//		mesh, 
//		IOobject::MUST_READ,
//		IOobject::AUTO_WRITE
//	   ),
//	   mesh
//	);
// End of new -----------------------------------------------------------

    // Need to store rho for ddt(rho, U)
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        alpha1*rho1 + alpha2*rho2 + alpha3*rho3 +alpha4*rho4,
        alpha1.boundaryField().types()
    );
    rho.oldTime();


    // Mass flux
    // Initialisation does not matter because rhoPhi is reset after the
    // alpha solution before it is used in the U equation.
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rhoPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho1*phi
    );


    // Construct interface from alpha distribution
    threePhaseInterfaceProperties interface(threePhaseProperties);


    // Construct incompressible turbulence model
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, threePhaseProperties)
    );


    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*gh
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }

    fv::IOoptionList fvOptions(mesh);
