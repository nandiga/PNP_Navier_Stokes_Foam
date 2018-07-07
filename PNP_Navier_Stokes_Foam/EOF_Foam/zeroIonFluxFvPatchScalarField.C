/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "zeroIonFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvCFD.H"
#include "fvcSnGrad.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zeroIonFluxFvPatchScalarField::
zeroIonFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    ElectricFieldName_("undefined-ElectricFieldName"),
    Zi_(1.0)
{
gradient() = 0.0;
}

Foam::zeroIonFluxFvPatchScalarField::
zeroIonFluxFvPatchScalarField
(
    const zeroIonFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    ElectricFieldName_(ptf.ElectricFieldName_),
    Zi_(ptf.Zi_)
{}


Foam::zeroIonFluxFvPatchScalarField::
zeroIonFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict),
    ElectricFieldName_(dict.lookup("electricField"))
{
Zi_=readScalar(dict.lookup("Zi"));
gradient()=0.0;

   if (dict.found("gradient"))
    {
        fvPatchScalarField::operator=
        (
         scalarField("gradient", dict, p.size())
        );
     }
    else
    {
       fvPatchScalarField::operator=(gradient());
    }

}

Foam::zeroIonFluxFvPatchScalarField::
zeroIonFluxFvPatchScalarField
(
    const zeroIonFluxFvPatchScalarField& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf),
    ElectricFieldName_(tppsf.ElectricFieldName_),
    Zi_(tppsf.Zi_)
{}


Foam::zeroIonFluxFvPatchScalarField::
zeroIonFluxFvPatchScalarField
(
    const zeroIonFluxFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF),
    ElectricFieldName_(tppsf.ElectricFieldName_),
    Zi_(tppsf.Zi_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zeroIonFluxFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
}


void Foam::zeroIonFluxFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);
}


void Foam::zeroIonFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

   const vectorField& electricField = 
        patch().lookupPatchField<volVectorField, scalar>(ElectricFieldName_);

   const fvPatchScalarField& Cw = *this;

//    *this[patchI].type() == "fixedGradient";

    gradient()= Cw*Zi_*electricField & patch().nf(); 

    fixedGradientFvPatchScalarField::updateCoeffs(); 
}


void Foam::zeroIonFluxFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedGradientFvPatchScalarField::write(os);
    os.writeKeyword("Zi") << Zi_ << token::END_STATEMENT << nl;
    os.writeKeyword("ElectricFieldName") << ElectricFieldName_ << token::END_STATEMENT << nl;
writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        zeroIonFluxFvPatchScalarField
    );
}

// ************************************************************************* //
