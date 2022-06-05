/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 Tomislav Maric, Tobias Tolle, TU Darmstadt
                       Anja Lippert, BOSCH CR 
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

#include "pandoraCurvatureAverageExtension.hpp"

#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "fvc.H"
#include "processorFvPatchField.H"

namespace Foam {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pandoraCurvatureAverageExtension, false);
addToRunTimeSelectionTable(pandoraCurvatureExtension, pandoraCurvatureAverageExtension, Dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
pandoraCurvatureAverageExtension::pandoraCurvatureAverageExtension()
:
    pandoraCurvatureExtension{dictionary{}},
    nExtensionIterations_{3}
{}


pandoraCurvatureAverageExtension::pandoraCurvatureAverageExtension(const dictionary& dict)
:
    pandoraCurvatureExtension{dict},
    nExtensionIterations_{dict.get<label>("nExtensionIterations")}
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void pandoraCurvatureAverageExtension::extend(volScalarField& curvature, const volScalarField& isInterfaceCell)
{
    const auto& mesh = curvature.mesh();
    const auto& owner = mesh.owner();
    const auto& neighbour = mesh.neighbour();

    // Check for compatible size
    if (isInterfaceCell.size() != mesh.nCells())
    {
        // Error
    }

    forAll(curvature, cid)
    {
        if (isInterfaceCell[cid] != 1.0)
        {
            curvature[cid] = tagValue_; 
        }
    }
    curvature.correctBoundaryConditions();

    auto faceCurvatureTmp = fvc::interpolate(curvature);
    auto& faceCurvature = faceCurvatureTmp.ref();
    labelField faceCount{curvature.size(), 0};
    scalarField surfaceSum{curvature.size(), 0.0};

    for (label iter = 0; iter != nExtensionIterations_; ++iter)
    {
        // Propagate curvature values to faces
        forAll(neighbour, fid)
        {
            if (curvature[owner[fid]] == tagValue_ && curvature[neighbour[fid]] == tagValue_)
            {
                faceCurvature[fid] = tagValue_;   
            }
            else if (curvature[owner[fid]] != tagValue_ && curvature[neighbour[fid]] == tagValue_)
            {
                faceCurvature[fid] = curvature[owner[fid]];
            }
            else if (curvature[owner[fid]] == tagValue_ && curvature[neighbour[fid]] != tagValue_)
            {
                faceCurvature[fid] = curvature[neighbour[fid]];
            }
        }

        auto& meshBoundaries = curvature.boundaryFieldRef();

        for (auto& mBoundary : meshBoundaries)
        {
            if (isType<processorFvPatch>(mBoundary.patch()))
            {
                // Ensure processor neighbour fields are up-to-date
                mBoundary.initEvaluate();
                mBoundary.evaluate();

                const auto& pPatch = mBoundary.patch();
                const auto& faceToCell = pPatch.faceCells();
                auto neighbourValuesTmp = mBoundary.patchNeighbourField();
                const auto& neighbourValues = neighbourValuesTmp.cref();

                forAll(mBoundary, I)
                {
                    if (curvature[faceToCell[I]] == tagValue_ && neighbourValues[I] == tagValue_)
                    {
                        mBoundary[I] = tagValue_;
                    }
                    else if (curvature[faceToCell[I]] != tagValue_ && neighbourValues[I] == tagValue_)
                    {
                        mBoundary[I] = curvature[faceToCell[I]];
                    }
                    else if (curvature[faceToCell[I]] == tagValue_ && neighbourValues[I] != tagValue_)
                    {
                        mBoundary[I] = neighbourValues[I];
                    }
                }
            }
        }

        // Compute surface sums and face count for cells
        faceCount = 0; 
        surfaceSum = 0.0;

        forAll(neighbour, fid)
        {
            if (faceCurvature[fid] != tagValue_)
            {
                faceCount[owner[fid]] += 1;
                surfaceSum[owner[fid]] += faceCurvature[fid];
                faceCount[neighbour[fid]] += 1;
                surfaceSum[neighbour[fid]] += faceCurvature[fid];
            }
        }

        // TODO: iterate processor boundaries and repeat steps above
        for (const auto& mBoundary : meshBoundaries)
        {
            if (isType<processorFvPatch>(mBoundary.patch()))
            {
                const auto& pPatch = mBoundary.patch();
                const auto& faceToCell = pPatch.faceCells();

                forAll(mBoundary, I)
                {
                    if (mBoundary[I] != tagValue_)
                    {
                        faceCount[faceToCell[I]] += 1;
                        surfaceSum[faceToCell[I]] += mBoundary[I];
                    }
                }
            }
        }

        forAll(curvature, cid)
        {
            if ((curvature[cid] == tagValue_) && (faceCount[cid] > 0))
            {
                curvature[cid] = surfaceSum[cid] / faceCount[cid];
            }
        }
        curvature.correctBoundaryConditions();
    }

    // Reset cells with tag value to zero to avoid explosions
    forAll(curvature, cid)
    {
        if (curvature[cid] == tagValue_)
        {
            curvature[cid] = 0.0;
        }
    }
    curvature.correctBoundaryConditions();

    auto& meshBoundaries = curvature.boundaryFieldRef();
    for (auto& mBoundary : meshBoundaries)
    {
        if (isType<processorFvPatch>(mBoundary.patch()))
        {
            forAll(mBoundary, I)
            {
                if (mBoundary[I] == tagValue_)
                {
                    mBoundary[I] = 0.0;
                }
            }
        }
    }
    curvature.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
