/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    threePhaseMixture

Description

SourceFiles
    threePhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef threePhaseMixture_H
#define threePhaseMixture_H

#include "incompressible/transportModel/transportModel.H"
#include "IOdictionary.H"
#include "incompressible/viscosityModels/viscosityModel/viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class threePhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class threePhaseMixture
:
    public IOdictionary,
    public transportModel
{
    // Private data

        word phase1Name_;
        word phase2Name_;
        word phase3Name_;
	word phase4Name_;

        volScalarField alpha1_;
        volScalarField alpha2_;
        volScalarField alpha3_;
	volScalarField alpha4_; //new

        const volVectorField& U_;
        const surfaceScalarField& phi_;

        volScalarField nu_;

        autoPtr<viscosityModel> nuModel1_;
        autoPtr<viscosityModel> nuModel2_;
        autoPtr<viscosityModel> nuModel3_;
        autoPtr<viscosityModel> nuModel4_;

        dimensionedScalar rho1_;
        dimensionedScalar rho2_;
        dimensionedScalar rho3_;
        dimensionedScalar rho4_;


    // Private Member Functions

        //- Calculate and return the laminar viscosity
        void calcNu();


public:

    // Constructors

        //- Construct from components
        threePhaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    ~threePhaseMixture()
    {}


    // Member Functions

        const word phase1Name() const
        {
            return phase1Name_;
        }

        const word phase2Name() const
        {
            return phase2Name_;
        }

        const word phase3Name() const
        {
            return phase3Name_;
        }

        const word phase4Name() const
        {
            return phase4Name_;
        }

        const volScalarField& alpha1() const
        {
            return alpha1_;
        }

        volScalarField& alpha1()
        {
            return alpha1_;
        }

        const volScalarField& alpha2() const
        {
            return alpha2_;
        }

        volScalarField& alpha2()
        {
            return alpha2_;
        }

        const volScalarField& alpha3() const
        {
            return alpha3_;
        }

        volScalarField& alpha3()
        {
            return alpha3_;
        }

        const volScalarField& alpha4() const //new
        {
            return alpha4_;
        }

        volScalarField& alpha4()
        {
            return alpha4_;
        }

        //- Return const-access to phase1 density
        const dimensionedScalar& rho1() const
        {
            return rho1_;
        }

        //- Return const-access to phase2 density
        const dimensionedScalar& rho2() const
        {
            return rho2_;
        };

        //- Return const-access to phase3 density
        const dimensionedScalar& rho3() const
        {
            return rho3_;
        };

        //- Return const-access to phase4 density
        const dimensionedScalar& rho4() const
        {
            return rho4_;
        };

        //- Return the velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the flux
        const surfaceScalarField& phi() const
        {
            return phi_;
        }

        //- Return const-access to phase1 viscosityModel
        const viscosityModel& nuModel1() const
        {
            return nuModel1_();
        }

        //- Return const-access to phase2 viscosityModel
        const viscosityModel& nuModel2() const
        {
            return nuModel2_();
        }

        //- Return const-access to phase3 viscosityModel
        const viscosityModel& nuModel3() const
        {
            return nuModel3_();
        }

        //- Return const-access to phase4 viscosityModel
        const viscosityModel& nuModel4() const
        {
            return nuModel4_();
        }

        //- Return the dynamic laminar viscosity
        tmp<volScalarField> mu() const;

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> muf() const;

        //- Return the kinematic laminar viscosity
        tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Return the face-interpolated dynamic laminar viscosity
        tmp<surfaceScalarField> nuf() const;

        //- Correct the laminar viscosity
        void correct()
        {
            calcNu();
        }

        //- Read base transportProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
