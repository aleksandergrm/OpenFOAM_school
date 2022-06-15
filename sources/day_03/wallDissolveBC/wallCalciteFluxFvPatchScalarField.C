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

#include "wallCalciteFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

wallCalciteFluxFvPatchScalarField::wallCalciteFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    diffusion_(0.0),
    k1_(0.0),
    alpha_(0.0),
    Ceq_(0.0)
{ }


wallCalciteFluxFvPatchScalarField::wallCalciteFluxFvPatchScalarField
(
    const wallCalciteFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    diffusion_(ptf.diffusion_),
    k1_(ptf.k1_),
    alpha_(0.0),
    Ceq_(ptf.Ceq_)
{ }


wallCalciteFluxFvPatchScalarField::wallCalciteFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    diffusion_(readScalar(dict.lookup("diffusion"))),
    k1_(readScalar(dict.lookup("k1"))),
    alpha_(readScalar(dict.lookup("alpha"))),
    Ceq_(readScalar(dict.lookup("Ceq")))
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    if (dict.found("gradient"))
    {
        gradient() = scalarField("gradient", dict, p.size());
        fixedGradientFvPatchScalarField::updateCoeffs();
        fixedGradientFvPatchScalarField::evaluate();
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


wallCalciteFluxFvPatchScalarField::wallCalciteFluxFvPatchScalarField
(
    const wallCalciteFluxFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf),
    diffusion_(wbppsf.diffusion_),
    k1_(wbppsf.k1_),
    alpha_(wbppsf.alpha_),
    Ceq_(wbppsf.Ceq_)
{ }


wallCalciteFluxFvPatchScalarField::wallCalciteFluxFvPatchScalarField
(
    const wallCalciteFluxFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    diffusion_(wbppsf.diffusion_),
    k1_(wbppsf.k1_),
    alpha_(wbppsf.alpha_),
    Ceq_(wbppsf.Ceq_)
{ }


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void wallCalciteFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<scalar>& C = patch().lookupPatchField<volScalarField, scalar>("T");

    scalarField Ceq = scalarField(C.size(), Ceq_);

    gradient() = k1_/diffusion_*(1 - C/Ceq_);

    fixedGradientFvPatchScalarField::updateCoeffs();
}


// Write
void wallCalciteFluxFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("diffusion") << diffusion_ << token::END_STATEMENT << nl;
    os.writeKeyword("k1") << k1_ << token::END_STATEMENT << nl;
    os.writeKeyword("alpha") << alpha_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ceq") << Ceq_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, wallCalciteFluxFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
