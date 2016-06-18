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

Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de
    Mathematical Modelling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany

\*---------------------------------------------------------------------------*/

namespace Foam {
    namespace GeometricalTransport {

template<typename MeshIntersection, typename GeomInterface>
tmp<volScalarField> volSymmDiff(
    MeshIntersection const& meshIntersection, 
    GeomInterface const& interface,
    const volScalarField& alpha1
)
{
    const auto& runTime = alpha1.time(); 
    const auto& mesh = alpha1.mesh(); 

    tmp<volScalarField> volSymmDiffTmp(
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
    volScalarField& volSymmDiff = volSymmDiffTmp(); 

    Info << "Computing the volume of symmetric difference..." << endl;
    runTime.cpuTimeIncrement();

    typedef pointVectorVector polyhedron; 
    typedef PolygonSequenceIntersection<polyhedron> polyhedronIntersection;

    // Get a copy of the list of polyhedrons from the tool mesh.
    const auto& cellPolyhedra = meshIntersection.cellPolyhedra();
    const auto& V = mesh.V(); 
    
    // Loop through interface cells.
    forAll (alpha1, cellI)
    {
        // TODO: Initialize geomControl and get tolerance.
        if ((alpha1[cellI] > 1e-09) && (alpha1[cellI] < (1-1e-09)))
        {
            // Calculate the polyhedron.
            // Build an outward oriented halfspace from the interface element.
            auto hspace = build<halfspace>(interface[cellI]);

            // Positive symmetric difference error.
            scalar Vpos = 0;
            // Negative symmetric difference error.
            scalar VsdNeg = 0;

            // Compute the volume of symmetric difference error.  
            for (const auto& cellPoly : cellPolyhedra[cellI])
            {
                hspace.flip();
                auto intersection = intersect<polyhedronIntersection>(cellPoly,hspace);
                VsdNeg += volume(intersection.polyhedron()); 

                hspace.flip();
                intersection = intersect<polyhedronIntersection>(cellPoly,hspace);
                Vpos += volume(intersection.polyhedron()); 
            } // End loop through tool mesh polyhedra.
            volSymmDiff[cellI] += VsdNeg + ((alpha1[cellI] * V[cellI]) - Vpos); 

        } // End loop through interface cells.
    } // End loop through all cells.
    Info << "Done in " << runTime.cpuTimeIncrement() << " seconds. " << endl;

    return volSymmDiff; 
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace GeometricalTransport 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
