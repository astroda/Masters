{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");
    
    //word alpha3Scheme("div(phi,alpha3)");
    //word alpha4Scheme("div(phi,alpha4)");


    surfaceScalarField phir
    (
        IOobject
        (
            "phir",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        interface.cAlpha()*mag(phi/mesh.magSf())*interface.nHatf()
    );

    for (int gCorr=0; gCorr<nAlphaCorr; gCorr++)
    {
        // Create the limiter to be used for all phase-fractions
        scalarField allLambda(mesh.nFaces(), 1.0);

        // Split the limiter into a surfaceScalarField
        slicedSurfaceScalarField lambda
        (
            IOobject
            (
                "lambda",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimless,
            allLambda,
            false   // Use slices for the couples
        );


        // Create the complete convection flux for alpha1
        surfaceScalarField phiAlpha1
        (
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha1,
                alpharScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, alpha3, alpharScheme),
                alpha1,
                alpharScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, alpha4, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        // Create the bounded (upwind) flux for alpha1
        surfaceScalarField phiAlpha1BD
        (
            upwind<scalar>(mesh, phi).flux(alpha1)
        );

        // Calculate the flux correction for alpha1
        phiAlpha1 -= phiAlpha1BD;

        // Calculate the limiter for alpha1
        MULES::limiter
        (
            allLambda,
            1.0/runTime.deltaT().value(),
            geometricOneField(),
            alpha1,
            phiAlpha1BD,
            phiAlpha1,
            zeroField(),
            zeroField(),
            1,
            0,
            3
        );

        // Create the complete flux for alpha2
        surfaceScalarField phiAlpha2
        (
            fvc::flux
            (
                phi,
                alpha2,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(phir, alpha1, alpharScheme),
                alpha2,
                alpharScheme
            )
        );

//New ----------------------------------------------------------------------------------
	surfaceScalarField phiAlpha3
        (
            fvc::flux
            (
                phi,
                alpha3,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(phir, alpha1, alpharScheme),
                alpha3,
                alpharScheme
            )
        );
	
	//UW bounded flux for alpha 3
	surfaceScalarField phiAlpha3BD
        (
            upwind<scalar>(mesh, phi).flux(alpha3)
        );
        // Calculate the flux correction for alpha3
        phiAlpha3 -= phiAlpha3BD;

	MULES::limiter
        (
            allLambda,
            1.0/runTime.deltaT().value(),
            geometricOneField(),
            alpha3,
            phiAlpha3BD,
            phiAlpha3,
            zeroField(),
            zeroField(),
            1,
            0,
            3
        );
//--------------------------------------------------------------------------------------

        // Create the bounded (upwind) flux for alpha2
        surfaceScalarField phiAlpha2BD
        (
            upwind<scalar>(mesh, phi).flux(alpha2)
        );

        // Calculate the flux correction for alpha2
        phiAlpha2 -= phiAlpha2BD;

        // Further limit the limiter for alpha2
        MULES::limiter
        (
            allLambda,
            1.0/runTime.deltaT().value(),
            geometricOneField(),
            alpha2,
            phiAlpha2BD,
            phiAlpha2,
            zeroField(),
            zeroField(),
            1,
            0,
            3
        );

        // Construct the limited fluxes
        phiAlpha1 = phiAlpha1BD + lambda*phiAlpha1;
        phiAlpha2 = phiAlpha2BD + lambda*phiAlpha2;
        phiAlpha3 = phiAlpha3BD + lambda*phiAlpha3; //new

        // Solve for alpha1
        solve(fvm::ddt(alpha1) + fvc::div(phiAlpha1));

        // Create the diffusion coefficients for alpha2<->alpha3
        //volScalarField Dc23(D23*max(alpha3, scalar(0))*pos(alpha2));
        //volScalarField Dc32(D23*max(alpha2, scalar(0))*pos(alpha3));

        //volScalarField Dc24(D24*max(alpha4, scalar(0))*pos(alpha2));
        //volScalarField Dc42(D24*max(alpha2, scalar(0))*pos(alpha4));

	//here pos means that if alpha2,3 is present then it is 1, if not then it is 0

        // Add the diffusive flux for alpha3->alpha2
        // replce with source terms in alpha2Eqn:
	//phiAlpha2 -= (fvc::snGrad(alpha2)*D2 - fvc::snGrad(alpha3)*D3 - fvc::snGrad(alpha4)*D4)*mesh.magSf()*fvc::snGrad(alpha1);
	//phiAlpha3 -= (-fvc::snGrad(alpha2)*D2 + fvc::snGrad(alpha3)*D3 - fvc::snGrad(alpha4)*D4)*mesh.magSf()*fvc::snGrad(alpha1); 
        //phiAlpha3 -= D3*mesh.magSf()*fvc::snGrad(alpha1); //new
	//Using alpha here works because its equal to the gradient	
	//phiAlpha2 -= fvc::interpolate(alpha2*D2 - alpha3*D3 - alpha4*D4)*mesh.magSf()*fvc::snGrad(alpha1); *******************************************************
	//phiAlpha3 -= fvc::interpolate(-alpha2*D2 + alpha3*D3 - alpha4*D4)*mesh.magSf()*fvc::snGrad(alpha1); ******************************************************

	//Source term definition used to couple/connect the phases:
        /*volScalarField source2
        (
	  fvc::laplacian(D2, alpha2)
	//- fvc::div(phiAlpha2) 
        );        

	volScalarField source3 //volScalarField source3
        (
	  fvc::laplacian(D3, alpha3)
	//- fvc::div(phiAlpha3)     //was div(phi, alpha3)
        );

        volScalarField source4
        (
	  fvc::laplacian(D4, alpha4)
	//- fvc::div(phi, alpha4) 
        );

	volVectorField flux2
	(
	 - ((alpha3/D34 + alpha2/D23 +alpha4/D34)*fvc::grad(alpha2) + alpha2*(1/D23 - 1/D24)*fvc::grad(alpha3))/((alpha2/D24 + alpha3/D23 + alpha4/D24)*(alpha3/D34 + alpha2/D23 + alpha4/D34) - (alpha2*alpha3*(1/D23-1/D24)*(1/D23 - 1/D34)))
	);


	volVectorField flux3
	(
	 - ((alpha2/D24 + alpha3/D23 +alpha4/D24)*fvc::grad(alpha3) + alpha3*(1/D23 - 1/D34)*fvc::grad(alpha2))/((alpha2/D24 + alpha3/D23 + alpha4/D24)*(alpha3/D34 + alpha2/D23 + alpha4/D34) - (alpha2*alpha3*(1/D23-1/D24)*(1/D23 - 1/D34)))

	);*/
	
	volScalarField S
	(

	  D34*pos(alpha2) + D24*pos(alpha3) + D23*pos(alpha4) + D23/1000

	);

	volVectorField flux2
	(
	  (fvc::grad(alpha2)*D24*(alpha2*D34+(1-alpha2)*D23)/S + fvc::grad(alpha3)*alpha2*D34*(D24-D23)/S)*(1-alpha1)
	);

	volVectorField flux3
	(
	  (fvc::grad(alpha2)*alpha3*D24*(D34-D23)/S + fvc::grad(alpha3)*D34*(alpha3*D24+(1-alpha3)*D23)/S)*(1-alpha1) //error in the second grad term
	);


        // Solve for alpha2
        fvScalarMatrix alpha2Eqn
        (
           fvm::ddt(alpha2)
         + fvc::div(phiAlpha2)
	 - fvc::div(flux2)
         //- fvm::laplacian(D2, alpha2)
	 //+ source3 
	 //+ source4
        );
        alpha2Eqn.solve();



	// Solve for alpha3 ------------------------------
	fvScalarMatrix alpha3Eqn
	(
	   fvm::ddt(alpha3)
	 + fvc::div(phiAlpha3)
	 - fvc::div(flux3)
	 //- fvm::laplacian(D3, alpha3)
	 //+ source2
	 //+ source4

	);
	alpha3Eqn.solve();

	// Solve for alpha4 ------------------------------
	/*fvScalarMatrix alpha4Eqn
	(
	   fvm::ddt(alpha4)
	 + fvc::div(phi, alpha4)
	 - fvm::laplacian(D4, alpha4)
	 - source2
	 - source3

	);
	alpha4Eqn.solve();
	*/

	
	
	alpha4 = 1.0 - alpha1 - alpha2 - alpha3;

        // Construct the complete mass flux
        rhoPhi =			//surfaceScalarField
              rho1*fvc::flux(phi, alpha1, alphaScheme)		//dimScalaraField*surfaceScalarField*volScalarField
              + rho2*fvc::flux(phi, alpha2, alphaScheme)
              + rho3*fvc::flux(phi, alpha3, alphaScheme)
              + rho4*fvc::flux(phi, alpha4, alphaScheme);	//+ rho4*phi*alpha4;
	// End of New ------------------------------------
	
	//constrain the alpha fields and cut off small numerical errors
	alpha1 = max(scalar(0), min(scalar(1), alpha1));
	alpha2 = max(scalar(0), min(scalar(1), alpha2));
	alpha3 = max(scalar(0), min(scalar(1), alpha3));
	alpha4 = max(scalar(0), min(scalar(1), alpha4));
	

    }

    Info<< "Air phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;

    Info<< "Liquid phase volume fraction = "
        << alpha2.weightedAverage(mesh.V()).value()
        << "  Min(alpha2) = " << min(alpha2).value()
        << "  Max(alpha2) = " << max(alpha2).value()
        << endl;

    Info<< "Other liquid phase volume fraction = "
	<< alpha3.weightedAverage(mesh.V()).value()
        << "  Min(alpha3) = " << min(alpha3).value()
        << "  Max(alpha3) = " << max(alpha3).value()
        << endl;

    Info<< "Other liquid phase volume fraction = "
	<< alpha4.weightedAverage(mesh.V()).value()
        << "  Min(alpha4) = " << min(alpha4).value()
        << "  Max(alpha4) = " << max(alpha4).value()
        << endl;

}
