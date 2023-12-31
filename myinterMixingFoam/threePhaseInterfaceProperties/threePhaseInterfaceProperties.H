/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    threePhaseInterfaceProperties

Description
    Properties to aid interFoam :
    1. Correct the alpha boundary condition for dynamic contact angle.
    2. Calculate interface curvature.

SourceFiles
    threePhaseInterfaceProperties.C

\*---------------------------------------------------------------------------*/

#ifndef threePhaseInterfaceProperties_H
#define threePhaseInterfaceProperties_H

#include "threePhaseMixture.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class threePhaseInterfaceProperties Declaration
\*---------------------------------------------------------------------------*/

class threePhaseInterfaceProperties
{
    // Private data

        const threePhaseMixture& mixture_;

        //- Compression coefficient
        scalar cAlpha_;

        //- Surface tension 1-2
        dimensionedScalar sigma12_;

        //- Surface tension 1-3
        dimensionedScalar sigma13_;

        //- Surface tension 1-4
        dimensionedScalar sigma14_;

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;

        surfaceScalarField nHatf_;
        volScalarField K_;


    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        threePhaseInterfaceProperties(const threePhaseInterfaceProperties&);
        void operator=(const threePhaseInterfaceProperties&);

        //- Correction for the boundary condition on the unit normal nHat on
        //  walls to produce the correct contact dynamic angle
        //  calculated from the component of U parallel to the wall
        void correctContactAngle
        (
            surfaceVectorField::GeometricBoundaryField& nHat
        ) const;

        //- Re-calculate the interface curvature
        void calculateK();


public:

    //- Conversion factor for degrees into radians
    static const scalar convertToRad;


    // Constructors

        //- Construct from volume fraction field alpha and IOdictionary
        threePhaseInterfaceProperties(const threePhaseMixture& mixture);


    // Member Functions

        scalar cAlpha() const
        {
            return cAlpha_;
        }

        const dimensionedScalar& deltaN() const
        {
            return deltaN_;
        }

        const surfaceScalarField& nHatf() const
        {
            return nHatf_;
        }

        const volScalarField& K() const
        {
            return K_;
        }

        tmp<volScalarField> sigma() const
        {
            volScalarField limitedAlpha2(max(mixture_.alpha2(), scalar(0)));
            volScalarField limitedAlpha3(max(mixture_.alpha3(), scalar(0)));
            volScalarField limitedAlpha4(max(mixture_.alpha4(), scalar(0)));

            return
                (limitedAlpha2*sigma12_ + limitedAlpha3*sigma13_ + limitedAlpha4*sigma14_)
               /(limitedAlpha2 + limitedAlpha3 + limitedAlpha4 + SMALL);
        }

        tmp<volScalarField> sigmaK() const
        {
            return sigma()*K_;
        }

        void correct()
        {
            calculateK();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
