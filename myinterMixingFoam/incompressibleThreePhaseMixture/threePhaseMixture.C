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

\*---------------------------------------------------------------------------*/

#include "threePhaseMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

//- Calculate and return the laminar viscosity
void Foam::threePhaseMixture::calcNu()
{
    nuModel1_->correct();
    nuModel2_->correct();
    nuModel3_->correct();
    nuModel3_->correct(); //new

    // Average kinematic viscosity calculated from dynamic viscosity
    nu_ = mu()/(alpha1_*rho1_ + alpha2_*rho2_ + alpha3_*rho3_ + alpha4_*rho4_); //added kast term
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::threePhaseMixture::threePhaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    phase1Name_(wordList(lookup("phases"))[0]),
    phase2Name_(wordList(lookup("phases"))[1]),
    phase3Name_(wordList(lookup("phases"))[2]),
    phase4Name_(wordList(lookup("phases"))[3]),

    alpha1_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase1Name_),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    alpha2_
    (
        IOobject
        (
            IOobject::groupName("alpha", phase2Name_),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    alpha3_ 
    (
        IOobject
        (
            IOobject::groupName("alpha", phase3Name_),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    alpha4_ //new
    (
        IOobject
        (
            IOobject::groupName("alpha", phase4Name_),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    U_(U),
    phi_(phi),

    nu_
    (
        IOobject
        (
            "nu",
            U.time().timeName(),
            U.db()
        ),
        U.mesh(),
        dimensionedScalar("nu", dimensionSet(0, 2, -1, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    ),

    nuModel1_
    (
        viscosityModel::New
        (
            "nu1",
            subDict(phase1Name_),
            U,
            phi
        )
    ),
    nuModel2_
    (
        viscosityModel::New
        (
            "nu2",
            subDict(phase2Name_),
            U,
            phi
        )
    ),
    nuModel3_
    (
        viscosityModel::New
        (
            "nu3",
            subDict(phase3Name_),
            U,
            phi
        )
    ),
    nuModel4_				//new
    (
        viscosityModel::New
        (
            "nu4",
            subDict(phase4Name_),
            U,
            phi
        )
    ),

    rho1_(nuModel1_->viscosityProperties().lookup("rho")),
    rho2_(nuModel2_->viscosityProperties().lookup("rho")),
    rho3_(nuModel3_->viscosityProperties().lookup("rho")),
    rho4_(nuModel4_->viscosityProperties().lookup("rho")) //new
{
    alpha4_ == 1.0 - alpha1_ - alpha2_ - alpha3_;
    calcNu();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::threePhaseMixture::mu() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            "mu",
            alpha1_*rho1_*nuModel1_->nu()
          + alpha2_*rho2_*nuModel2_->nu()
          + alpha3_*rho3_*nuModel3_->nu()
          + alpha4_*rho4_*nuModel4_->nu()
        )
    );
}


Foam::tmp<Foam::surfaceScalarField> Foam::threePhaseMixture::muf() const
{
    surfaceScalarField alpha1f(fvc::interpolate(alpha1_));
    surfaceScalarField alpha2f(fvc::interpolate(alpha2_));
    surfaceScalarField alpha3f(fvc::interpolate(alpha3_));
    surfaceScalarField alpha4f(fvc::interpolate(alpha4_));

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "mu",
            alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
          + alpha2f*rho2_*fvc::interpolate(nuModel2_->nu())
          + alpha3f*rho3_*fvc::interpolate(nuModel3_->nu())
          + alpha4f*rho4_*fvc::interpolate(nuModel4_->nu()) //new new new
        )
    );
}


Foam::tmp<Foam::surfaceScalarField> Foam::threePhaseMixture::nuf() const
{
    surfaceScalarField alpha1f(fvc::interpolate(alpha1_));
    surfaceScalarField alpha2f(fvc::interpolate(alpha2_));
    surfaceScalarField alpha3f(fvc::interpolate(alpha3_));
    surfaceScalarField alpha4f(fvc::interpolate(alpha4_));

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "nu",
            (
                alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
              + alpha2f*rho2_*fvc::interpolate(nuModel2_->nu())
              + alpha3f*rho3_*fvc::interpolate(nuModel3_->nu())
              + alpha4f*rho4_*fvc::interpolate(nuModel4_->nu())
            )/(alpha1f*rho1_ + alpha2f*rho2_ + alpha3f*rho3_ + alpha4f*rho4_)
        )
    );
}


bool Foam::threePhaseMixture::read()
{
    if (transportModel::read())
    {
        if
        (
            nuModel1_().read(*this)
         && nuModel2_().read(*this)
         && nuModel3_().read(*this)
         && nuModel4_().read(*this)
        )
        {
            nuModel1_->viscosityProperties().lookup("rho") >> rho1_;
            nuModel2_->viscosityProperties().lookup("rho") >> rho2_;
            nuModel3_->viscosityProperties().lookup("rho") >> rho3_;
            nuModel4_->viscosityProperties().lookup("rho") >> rho4_;

            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
