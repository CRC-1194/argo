/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | Copyright 2011 Tomislav Maric 
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

Description
    Base model for computing divergence free fields: 
        * Cell-centered velocity field. 
        * Face-centered flux field. 
        * Point-centered velocity field.

Author
    Tomislav Maric maric@mma.tu-darmstadt.de 
    Mathematical Modeling and Analysis
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/

#include "divFreeFieldModel.H"
#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(divFreeFieldModel, 0);
defineRunTimeSelectionTable(divFreeFieldModel, dictionary)
defineRunTimeSelectionTable(divFreeFieldModel, time)

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

divFreeFieldModel::divFreeFieldModel
(
    const Time& time,
    const dictionary& dict
)
:
    time_(time),
    cellVelocityName_(dict.lookupOrDefault<word>("cellVelocityName", "U")),
    volFluxName_(dict.lookupOrDefault<word>("volFluxName", "phi")),
    pointVelocityName_(dict.lookupOrDefault<word>("pointVelocityName", ""))
{
    if (cellVelocityName_.empty() && volFluxName_.empty() && pointVelocityName_.empty())
        WarningIn("divFreeModel::divFreeModel") 
            << "No fields are selected to be used by the divFreeModel." << endl;  
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<divFreeFieldModel>
divFreeFieldModel::New
(
    const word& name,
    const Time& time,
    const dictionary& dict
)
{
    auto* ctorPtr = dictionaryConstructorTable(name);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "divFreeModel",
            name,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<divFreeFieldModel>(ctorPtr(time, dict));
}

autoPtr<divFreeFieldModel>
divFreeFieldModel::New
(
    const Time& time,
    const dictionary& dict
)
{
    const word name = dict.get<word>("type"); 

    return divFreeFieldModel::New(name, time, dict); 
}

autoPtr<divFreeFieldModel>
divFreeFieldModel::New
(
    const Time& time
)
{
    const dictionary& dict = time.controlDict().subDict("divFree");

    return divFreeFieldModel::New(time, dict); 
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Calculate and set the cell centered velocity field.
void divFreeFieldModel::setCellVelocity(const Time&)
{
    if (!cellVelocityName_.empty())
    {
        const objectRegistry& reg
        (
            time_.lookupObject<objectRegistry>(polyMesh::defaultRegion)
        );

        if (! reg.foundObject<volVectorField>(cellVelocityName_))
        {
            FatalErrorIn("divFreeFieldModel::setCellVelocity(const Time&)")
                << "Cell velocity field " 
                << cellVelocityName_ << " is not registered." << endl;
        }

        // Make the velocity field modifiable by casting away constness.
        volVectorField & U = const_cast<volVectorField&>
        (
            reg.lookupObject<volVectorField>(cellVelocityName_)
        );

        const fvMesh& mesh = U.mesh();
        // Cell centers.
        const volVectorField& C = mesh.C();
        scalar t = time_.timeOutputValue();

        // Set the internal cell centered velocity field.
        forAll (U, I)
        {
            U[I] = velocity(C[I], t);
        }

        // TODO: Overwrite the boundary field values for parallel process patches.
        //       Important for interpolation. TM. 
        // Boundary field update.
        U.correctBoundaryConditions(); 
    }
}

// Calculate and set the scalar flux field.
void divFreeFieldModel::setFlux(const Time&)
{
    if (!volFluxName_.empty())
    {
        // Get non-const access to the objects.
        const objectRegistry& reg
        (
            time_.lookupObject<objectRegistry>(polyMesh::defaultRegion)
        );

        if (! reg.foundObject<surfaceScalarField>(volFluxName_))
        {
            FatalErrorIn("divFreeFieldModel::setFlux(const Time&)")
                << "Volumetric flux field " 
                << volFluxName_ << " is not registered." << endl;
        }

        surfaceScalarField& phi = const_cast<surfaceScalarField&>
        (
            reg.lookupObject<surfaceScalarField>(volFluxName_)
        );

        const fvMesh& mesh = phi.mesh();
        // Face area normal vectors.
        const surfaceVectorField& Sf = mesh.Sf();
        // Face centers.
        const surfaceVectorField& Cf = mesh.Cf();
        // Get the time.
        scalar t = time_.timeOutputValue();

        // Set the internal volumetric flux.
        forAll (phi, I)
        {
            phi[I] = Sf[I] & velocity(Cf[I], t);
        }

        // Set the boundary volumetric flux.
        forAll(mesh.boundary(), patchi)
        {
            fvsPatchScalarField& phiBoundaryField  = phi.boundaryFieldRef()[patchi];
            const fvsPatchVectorField& SfBoundaryField  = Sf.boundaryField()[patchi];
            const fvsPatchVectorField& CfBoundaryField  = Cf.boundaryField()[patchi];

            forAll(phiBoundaryField, faceI)
            {
                phiBoundaryField[faceI] = SfBoundaryField[faceI] & 
                    velocity(CfBoundaryField[faceI], t);
            }
        }
    }
}

void divFreeFieldModel::setPointVelocity(pointVectorField& Up)
{
    const pointMesh& mesh = Up.mesh();
    const pointField& points = mesh.mesh().points();
    scalar t = time_.timeOutputValue();

    // Set the internal point velocity field.
    forAll (Up, I)
        Up[I] = velocity(points[I], t);

    // Set the boundary point velocity field.
    const pointVectorField::Boundary&  UpBoundary = Up.boundaryFieldRef(); 
    forAll (UpBoundary, patchI)
    {
        const auto& UpBpatch = UpBoundary[patchI].patch();
        const auto& meshPoints = UpBpatch.meshPoints(); 

        forAll(meshPoints, pointI)
            Up[pointI] = velocity(points[pointI], t);  
    }
}

// Calculate and set the point velocity field.
void divFreeFieldModel::setPointVelocity(const Time&)
{
    if (!pointVelocityName_.empty())
    {
        // Get access to the registered objects.
        const objectRegistry& reg
        (
            time_.lookupObject<objectRegistry>(polyMesh::defaultRegion)
        );

        if (! reg.foundObject<pointVectorField>(pointVelocityName_))
        {
            FatalErrorIn("divFreeFieldModel::setPointVelocity(const Time&)")
                << "Point velocity field " << pointVelocityName_ << " is not registered." << endl;
        }

        // Get non-const point velocity access.
        pointVectorField& Up = const_cast<pointVectorField&>
        (
            reg.lookupObject<pointVectorField>(pointVelocityName_)
        );

        setPointVelocity(Up); 
    }
}

vector divFreeFieldModel::velocity(point, scalar) const
{
    notImplemented("divFreeFieldModel::velocity(point, scalar) const");

    return vector(0,0,0);
}

bool divFreeFieldModel::execute() 
{
    setCellVelocity(time_);

    setFlux(time_);

    setPointVelocity(time_);

    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
