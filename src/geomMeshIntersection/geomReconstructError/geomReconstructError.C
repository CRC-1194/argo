/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::geomReconstructError

Description
    A class for computing the standard errors of the geometrical two phase
    flow methods (VoF, MoF).

Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de
    Mathematical Modelling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany

\*---------------------------------------------------------------------------*/

#include "geomReconstructError.H"
#include "polyhedron.H"

namespace Foam {
    namespace GeometricalTransport {

template<typename MeshIntersection, typename PhaseInterface>
tmp<volScalarField> volSymmDiff(
    MeshIntersection const& meshIntersection, 
    PhaseInterface const& interface,
    const volScalarField& alpha1
)
{
    const auto& runTime = alpha1.time(); 
    const auto& mesh = alpha1.mesh(); 


    tmp<volScalarField> volSymmDiff(
        new volScalarField(
            IOobject(
                "volSymmDiff", 
                runTime.timeName(), 
                mesh, 
                IOobject::NO_READ, 
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimVolume, 0)
        )
    );

    Info << "Computing the volume of symmetric difference error ..." << endl;
    runTime.cpuTimeIncrement();

    typedef pointVectorVector polyhedron; 
    typedef PointSequenceIntersection<polyhedron> polyhedronIntersection;
    
    // Loop through interface cells.
    forAll (alpha1, cellI)
    {
        // TODO: parametrize the alpha1_ tolerance.
        // Initialize geomControl, empty constructor, reads default dictionary.  
        //if ((alpha1[cellI] > 1e-09) && (alpha1[cellI] < (1-1e-09)))
        //{
            //// Calculate the polyhedron.
            
            //// Build an outward oriented halfspace from the interface element.
            //auto halfspace = build<halfspace>(interface_[cellI])
            //// Point it towards the first phase.
            //halfspace.flip(); 

            //// Build the polyhedron from the cell. 
            //auto cellPolyhedron = build<polyhedron>(cellI, mesh); 

            //// Intersect the cell with a halfspace.  
            //auto cellIntersection = 
                //intersect<polyhedronIntersection>(cellPolyhedron,halfspace);

            //// Get a copy of the list of polyhedrons from the tool mesh.
            //const auto& intersectionPolyhedra = meshIntersection.polyhedra();

            //// Positive symmetric difference error.
            //scalar VsdPos = 0;
            //// Negative symmetric difference error.
            //scalar VsdNeg = 0;

            //// Compute the volume of symmetric difference error.  
            ////std::list<polyhedron>::const_iterator pcellIt;
            ////pcellIt = intersectionPolys[cellI].begin();
            ////for (; pcellIt != intersectionPolys[cellI].end(); ++pcellIt)
            //for (const auto& interPolyhedron : intersectionPolyhedra)
            //{
                //// Clip the polyhedron with the cell.
                //polyhedron clippedPoly = pcellIt->intersect(cellPoly);

                //// Compute the volume of intersection with the interface
                //// polyhedron.

                //scalar Vintersect = clippedPoly.intersect(interfacePoly).mag();

                //// Add the volume of intersection to the positive error.
                //VsdPos += Vintersect; 

                //if (clippedPoly.intersects(interfacePlane))
                //{
                    //// Add the difference volume to the negative error.
                    ////VsdNeg += pcellIt->mag() - Vintersect; 
                    //VsdNeg += clippedPoly.mag() - Vintersect; 
                //}

            //} // End loop through tool mesh polyhedra.

            //VsdPos = interfacePoly.mag() - VsdPos;

            //volSymmDiff[cellI] = VsdPos + VsdNeg;

        //} // End loop through interface cells.
    } // End loop through all cells.

    Info << "Done in " << intersectedMesh_.baseTime().cpuTimeIncrement() 
        << " seconds. " << endl;

    return volSymmDiff; 
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace GeometricalTransport 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
