/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) Tomislav Maric and TU Darmstadt 
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

Description
    This application is only used for testing of the CCI mesh intersection
    algorithm. 

    Intersect two meshes and compute the resulting volume fraction field.
    
    Test description: 
    
    Use bounding boxes of both meshes to randomly place the tool mesh within
    the base mesh such that the tool mesh bounding box is contained within the
    base mesh bounding box.  
    
Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de, maric@mma.tu-darmstadt.de, tomislav.maric@gmx.com
    Mathematical Modeling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "Random.H"
#include "OFstream.H"

#include "geomMeshIntersection.H"

using namespace Foam;
using namespace GeometricalTransport;

#include <chrono>
#include <fstream>

using namespace std::chrono;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Intersect the base mesh with the tool mesh and store the volume fraction given by the intersection on the base mesh."
    );

    argList::addOption
    (
        "fieldName", 
        "alpha.water",
        "Name of the volume fraction field." 
    ); 

    argList::addOption
    (
        "toolCase", 
        "pathName",
        "Path to the tool case. Default '../toolCase'" 
    ); 

    argList::addOption
    (
        "dataFile", 
        "geomIntersectMeshes.dat",
        "Name of file used to save the reported errors." 
    ); 

    argList::addOption
    (
        "nIterations", 
        100,
        "Number of times the tool mesh is randomly positioned and intersected with the base mesh." 
    ); 

    argList args(argc, argv);

    // Read user-defined options.
    const word fieldName = args.optionLookupOrDefault<word>("fieldName", "alpha.water"); 
    const word dataFileName = args.optionLookupOrDefault<word>("dataFile", "cci.csv"); 
    const label nIterations = args.optionLookupOrDefault<label>("nIterations", 100);  

    high_resolution_clock::time_point p0 = high_resolution_clock::now();
    #include "createMeshes.H"

    Info<< "Reading field alpha\n" << endl;
    volScalarField alpha
    (
        IOobject
        (
            fieldName,
            runTimeBase.timeName(),
            baseMesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        baseMesh, 
        dimensionedScalar ("zero", dimless, 0)
    );
    high_resolution_clock::time_point p1 = high_resolution_clock::now();

    // RANDOM tool mesh positioning. 
    // Position the tool mesh centroid at the base mesh centroid. 
    label nToolCells = toolMesh.nCells(); 
    label nBaseCells = baseMesh.nCells();  
    Pstream::gather(nToolCells, sumOp<label>()); 
    Pstream::gather(nBaseCells, sumOp<label>()); 
    vector toolCenter = gSum(toolMesh.C()) / nToolCells;  
    const vector baseCenter = gSum(baseMesh.C()) / nBaseCells;  
    toolMesh.movePoints(toolMesh.points() + baseCenter - toolCenter); 

    boundBox baseBox = baseMesh.bounds(); 
    const boundBox toolBox = toolMesh.bounds(); 
    const Vector<label> solutionVector = baseMesh.solutionD(); 
    vector toolBoxDelta = toolBox.max() - toolBox.min(); 
    // Null the toolBoxDelta in a dimension that is not used in pseudo2D. 
    forAll(toolBoxDelta, I)
    {
        if (solutionVector[I] == -1)
            toolBoxDelta[I] *= 0;
    }
    point& baseMin = baseBox.min(); 
    point& baseMax = baseBox.max(); 
    // Shrink the base mesh box, so that the random position of the tool mesh
    // centroid is placed in a way that the base mesh bounding box contains  
    // the tool mesh bounding box. 
    baseMin += toolBoxDelta; 
    baseMax -= toolBoxDelta;  
    
    // Generate a random tool mesh centroid position within the shrunk base
    // mesh bounding box. 
    Random r(1e04); 

    // Open the error file for measurement output..
    OFstream errorFile(dataFileName); 
    // Nt : number of cells in the tool mesh 
    // Nb : number of cells in the tool mesh 
    // Ev : volume conservation error: 
    // |tool mesh volume from volume fraction - tool mesh volume | / tool mesh volume. 
    // Ti : initialization time, loading the mesh and fields, 
    // Te : execution time of the CCI mesh intersection operation.
    // Nx : total number of cells that are intersected: larger than the 
    //      number of interface cells, as bulk cells are intersected as well.
    // Ax : average number of intersections per intersected cell: average sum of 
    //      tool cell sizes per intersected cell.
    // Ni : number of interface cells, 
    // Nb : number of bulk cells.
    errorFile << "Nt,Nb,Ev,Ti,Te,Nx,Ax,Ni,Nb\n"; 


    for (int testI = 0; testI < nIterations; ++testI)
    {
        const vector randomPosition = r.position(baseMin,baseMax); 
        // Compute the displacement of the tool mesh with respect to the base mesh
        // centroid.
        vector randomDisplacement = randomPosition - baseCenter; 
        // Null the displacement in a dimension that is not used in pseudo2D. 
        forAll(randomDisplacement, I)
        {
            if (solutionVector[I] == -1)
                randomDisplacement[I] *= 0;
        }

        if (! baseBox.contains(randomPosition))
        {
            Info << baseMin << ":" << baseMax << ":" << randomPosition << nl;
            FatalErrorIn("voFoamTestCellCellIntersectMeshes")                                                 
                << "Random position outside of the base mesh bounding box." 
                << Foam::abort(FatalError);
        }

        // The toolMesh centroid now overlaps the base mesh centroid.
        toolMesh.movePoints(toolMesh.points() + randomDisplacement); 

        high_resolution_clock::time_point t0 = high_resolution_clock::now();
        geomMeshIntersection meshIntersection(baseMesh, toolMesh);
        meshIntersection.setVolFraction(alpha);
        high_resolution_clock::time_point t1 = high_resolution_clock::now();

        alpha.write(); 
        meshIntersection.report(Info, alpha); 

        const scalar Valpha = gSum((baseMesh.V() * alpha)()); 
        const scalar Vt = gSum(toolMesh.V());
        const scalar Ev = mag(Vt - Valpha) / Vt;  
        const scalar Ti = duration_cast<nanoseconds>(p1 - p0).count() / 1e09;
        const scalar Te = duration_cast<nanoseconds>(t1 - t0).count() / 1e09;

        label Ni = 0; 
        label Nb = 0; 

        forAll(alpha, cellI)
        {
            if ((alpha[cellI] > 0) && (alpha[cellI] < 1))
                ++Ni;  
            else 
                ++Nb; 
        }

        errorFile << toolMesh.nCells() << "," 
            << baseMesh.nCells() << ","
            << Ev << ","  
            << Ti << "," 
            << Te << ","
            << meshIntersection.Nx() << ","
            << meshIntersection.Ax() << ","
            << Ni << "," << Nb << nl;


        // Return the tool mesh to the original position.
        toolMesh.movePoints(toolMesh.points() - randomDisplacement); 
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
