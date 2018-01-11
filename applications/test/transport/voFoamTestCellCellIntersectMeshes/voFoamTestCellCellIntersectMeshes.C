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

Description
    Intersect two meshes and compute the resulting volume fraction field.

    This application is only used for testing. 
    
    It moves the tool mesh randomly before intersection to accurately measure
    the execution time. 

    Test description: 
        - The base mesh is either a 3D unit-domain [1x1x1] or a 2D unit-domain
          [1x1]. 
        - The tool mesh is initialized in the center of the base mesh.
        - The tool mesh is moved randomly. 
        - Built for cylinder and sphere tool meshes: random placement of other
          tool mesh shapes might cause the tool mesh not to be within the base
          mesh!

     Possible generalization:
        - Use bounding boxes of both meshes to randomly place the tool mesh
          within the base mesh.  
        - Works only for convex meshes anyway: performance can be measured 
          accurately also with unit-domains, so there is no need for this. 

Author

    Tomislav Maric
    maric@csi.tu-darmstadt.de
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

#include <omp.h>
#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    setMeshIntersectionArgs(argc,argv);

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

    #include "createMeshes.H"

    // RANDOM tool mesh positioning. 
    // Position the tool mesh centroid at the base mesh centroid. 
    // TODO: Parallelize, start with commented lines. TM. 
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
    Random r(3231271); 

    // Open the error file for measurement output..
    const word fieldName = args.optionLookupOrDefault<word>("fieldName", "alpha.water"); 
    const word dataFileName = args.optionLookupOrDefault<word>("dataFile", "geomIntersectMeshes.csv"); 
    const label nIterations = args.optionLookupOrDefault<label>("nIterations", 100);  
    OFstream errorFile(dataFileName); 
    errorFile << "Nt,Nb,Ev,Te\n"; 

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

        Info << "Random displacement = " << randomDisplacement << nl;
        // The toolMesh centroid now overlaps the base mesh centroid.
        toolMesh.movePoints(toolMesh.points() + randomDisplacement); 

        // TODO: Remove, used for testing.
        toolMesh.write(); 
        

        Info<< "Reading field alpha1\n" << endl;
        volScalarField alpha1
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

        #pragma omp single 
        {
            // Precalculate mesh geometry. Lazy evaluation triggers race conditions. TM.
            baseMesh.Sf(); 
            baseMesh.Cf();
            baseMesh.C(); 
        }

        const double t0 = omp_get_wtime(); 
        geomMeshIntersection intersect(baseMesh, toolMesh);
        intersect.setVolFraction(alpha1);
        const double t1 = omp_get_wtime(); 

        alpha1.write(); 
        intersect.report(Info, alpha1); 

        const scalar Vb = gSum((baseMesh.V() * alpha1)()); 
        const scalar Vt = sum(toolMesh.V()).value();  
        const scalar Ev = mag(Vt - Vb);  
        const double Te = t1 - t0;

        errorFile << toolMesh.nCells() << "," 
            << baseMesh.nCells() << ","
            << Ev << ","  
            << Te << "\n";  


        // Return the tool mesh to the original position.
        toolMesh.movePoints(toolMesh.points() - randomDisplacement); 
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
