/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 AUTHOR,AFFILIATION
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

#include "volumeFractionCalculator.hpp"

#include "IOobject.H"
#include "addToRunTimeSelectionTable.H" 
#include "dictionary.H"
#include "fvc.H"
#include "pointMesh.H"

namespace Foam::TriSurfaceImmersion {

    defineTypeNameAndDebug(volumeFractionCalculator, 0);
    defineRunTimeSelectionTable(volumeFractionCalculator, Dictionary)


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
volumeFractionCalculator::volumeFractionCalculator
(
    const dictionary& configDict,
    const fvMesh& mesh,
    const triSurface& surface
)
:
    mesh_{mesh},
    runTime_{mesh.time()},
    surface_{surface},
    // init fields here
    cellSqrSearchDist_
    {
        IOobject
        (
            "cellSqrSearchDist", 
            runTime_.timeName(), 
            mesh_, 
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::average(pow(mesh.deltaCoeffs(), -2))
    },
    cellSignedDist_
    {
        IOobject
        (
            "cellSignedDist", 
            runTime_.timeName(), 
            mesh, 
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("cellSignedDist", dimLength,0),
        "zeroGradient"
    },
    cellSignedDist0_{"cellSignedDist0", cellSignedDist_}, 
    faceSignedDist_{
        IOobject
        (
            "faceSignedDist", 
            runTime_.timeName(), 
            mesh, 
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("faceSignedDist", dimLength,0)
    },
    cellsToPointsInterp_{mesh},
    pMesh_{mesh_},
    // init fields here
    pointSqrSearchDist_ 
    {
        IOobject
        (
            "pointSqrSearchDist", 
            runTime_.timeName(), 
            mesh_, 
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),

        pMesh_,
        dimensionedScalar("pointSqrSearchDist", dimLength,0),
        "zeroGradient"
    },
    pointSignedDist_ 
    {
        IOobject
        (
            "pointSignedDist", 
            runTime_.timeName(), 
            mesh_, 
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh_,
        dimensionedScalar("pointSignedDist", dimLength,0),
        "zeroGradient"
    }
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
autoPtr<volumeFractionCalculator>
volumeFractionCalculator::New
(
    const dictionary& configDict,
    const fvMesh& mesh,
    const triSurface& surface
)
{
    const word name = configDict.get<word>("type");

    DictionaryConstructorTable::iterator cstrIter =
        DictionaryConstructorTablePtr_->find(name);

    if (cstrIter == DictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn (
            "volumeFractionCalculator::New(const word& name)"
        )   << "Unknown volumeFractionCalculator type "
            << name << nl << nl
            << "Valid volumeFractionCalculators are : " << endl
            << DictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<volumeFractionCalculator>(cstrIter()(configDict, mesh, surface));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void volumeFractionCalculator::printTypeName() const
{
    Info << "VoF calculator type: " << this->type() << endl;
}

void volumeFractionCalculator::writeFields() const
{
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


} // End namespace Foam::TriSurfaceImmersion

// ************************************************************************* //