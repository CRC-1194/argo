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
    Evaluates standard advection errors for interface capturing methods.

Author
    Tomislav Maric maric@mma.tu-darmstadt.de 
    Mathematical Modeling and Analysis
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/


#include "fieldErrorFunctionObject.H"
#include "addToRunTimeSelectionTable.H"
#include "error.H"
#include "fvc.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(fieldErrorFunctionObject, 0);
addToRunTimeSelectionTable(functionObject, fieldErrorFunctionObject, dictionary);

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void fieldErrorFunctionObject::calcErrorVolume()
{
    const volScalarField& field = getCurrentField();

    errorVolume_[time_.timeIndex()] =
        mag(sum(initialFieldPtr_() - mag((1.0 - alphaDropletPhase_) - field))).value() /
        initialFieldSum_;
}

void fieldErrorFunctionObject::calcErrorL1()
{
    const volScalarField& field = getCurrentField();

    const fvMesh& mesh = field.mesh(); 

    const scalarField& cellVolumes  = mesh.V();

    errorL1_[time_.timeIndex()] =
        sum(cellVolumes * mag(initialFieldPtr_() - field));

    // Correct when solving in 2D for the cell layer height, the L1 error
    // should then be SfMag  * (mag(initialFieldPtr_() - field)) and not
    // cellVolumes * (-//-) because the volume size scales with the height
    // of the cell layer. TM.  
    
    scalar cellLayerHeight = -GREAT; 

    const labelList& faceOwner = mesh.faceOwner(); 
    const cellList& cells = mesh.cells(); 
    const faceList& faces = mesh.faces(); 
    const pointField& points = mesh.points(); 

    forAll(mesh.boundary(), patchI)
    {
        const fvPatch& fvp = mesh.boundary()[patchI]; 
        // If there is a patch of type empty, we have a 2D case. 
        if (fvp.type() == "empty")
        {
            // Get the first face of the patch. 
            const face& firstFace = faces[fvp.start()]; 
            // Compute the *reversed* unit normal vector of the first patch face. 
            vector firstFaceNormal = -1 * firstFace.unitNormal(points); 
            // Compute the center of the first patch face. 
            const vector firstFaceCenter = firstFace.centre(points);  
            // Get the corresponding cell. 
            const cell& firstCell = cells[faceOwner[fvp.start()]]; 
            // Get the cell points. 
            const pointField firstCellPoints = firstCell.points(faces, points);  

            // Find the cell point with the maximal distance from the face center
            // in the direction of the *reversed* face normal vector. 
            forAll (firstCellPoints, pointI)
            {
                const point& cellPoint = firstCellPoints[pointI]; 
                scalar pointHeight = firstFaceNormal & (cellPoint - firstFaceCenter); 
                if (pointHeight > cellLayerHeight)
                    cellLayerHeight = pointHeight;
            }

            // Correct the L1 error by dividing it with the height of the cell layer. 
            // to make it proportional to the SfMag - face area. 
            errorL1_[time_.timeIndex()] /= cellLayerHeight; 
            break;
        }
    }

}

void fieldErrorFunctionObject::calcErrorL1norm()
{
    const volScalarField& field = getCurrentField();

    const scalarField& cellVolumes = field.mesh().V();

    errorL1norm_[time_.timeIndex()] =
        sum(mag(initialFieldPtr_() - field)).value() /
        sum(cellVolumes * mag(initialFieldPtr_()));
}

const volScalarField& fieldErrorFunctionObject::getCurrentField() const
{
    // Get the reference to the object registry.
    const objectRegistry& reg
    (
        time_.lookupObject<objectRegistry>(polyMesh::defaultRegion)
    );

    if (! reg.foundObject<volScalarField>(fieldName_))
    {
        FatalErrorIn("fieldErrorFunctionObject::start()")
            << "Volume fraction field " << fieldName_ << " is not registered." << endl;
    }

    // Get the reference to the volume fraction field.
    const volScalarField & field = reg.lookupObject<volScalarField>(fieldName_);

    return field;
}

void fieldErrorFunctionObject::adjustFieldSizes()
{
    errorVolume_.resize(time_.timeIndex() + 1);
    errorL1_.resize(time_.timeIndex() + 1);
    errorL1norm_.resize(time_.timeIndex() + 1);
    elapsedTime_.resize(time_.timeIndex() + 1);
}

void fieldErrorFunctionObject::setElapsedTime()
{
    Info << "Setting elapsed time to : " << time_.timeOutputValue() << endl;
    elapsedTime_[time_.timeIndex()] = time_.timeOutputValue();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fieldErrorFunctionObject::fieldErrorFunctionObject
(
    const word& name,
    const Time& time,
    const dictionary& dict
)
:
    functionObject{name},
    started_{false},
    alphaDropletPhase_{readScalar(dict.lookup("alphaDroplet"))},
    initialFieldPtr_{nullptr},
    time_{time},
    fieldName_{dict.lookup("fieldName")},
    errorVolume_{},
    errorL1_{},
    errorL1norm_{},
    elapsedTime_{},
    initialFieldSum_{0}
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

fieldErrorFunctionObject::~fieldErrorFunctionObject()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool fieldErrorFunctionObject::read(const dictionary& dict)
{
    fieldName_ = word(dict.lookup("fieldName"));
    alphaDropletPhase_ = readScalar(dict.lookup("alphaDroplet"));

    return execute();
}

bool fieldErrorFunctionObject::execute()
{
    if (!started_)
    {
        started_ = start();
    }

    adjustFieldSizes();

    setElapsedTime();

    calcErrorVolume();
    calcErrorL1();
    calcErrorL1norm();

    return true;
}

bool fieldErrorFunctionObject::start()
{
    // Get the reference to the initial volume fraction field.
    const volScalarField & field = getCurrentField();

    // Create the copy of the initial volume fraction field and store it for error calculation.
    initialFieldPtr_ = autoPtr<volScalarField>(
       new volScalarField(
           mag((1.0 - alphaDropletPhase_) - field)
       )
    );

    initialFieldSum_ = sum(mag((1.0 - alphaDropletPhase_) - field)).value();

    if (initialFieldSum_ <= SMALL)
    {
        FatalErrorIn ("fieldFunctionObject::start()")
            << "Field for which the errors are computed (" << fieldName_
            << ") has not been initialized." << abort(FatalError);
    }

    started_ = true;

    return true;
}

bool fieldErrorFunctionObject::end()
{
    // Write the error data.
    write();

    execute();

    return true;
}

bool fieldErrorFunctionObject::write()
{
    // Write the errors in the prefix/fieldErrors.dat file.
    fileName errorFileName (time_.rootPath() + "/" + time_.globalCaseName() 
                            + "/advectionErrors.dat");
    OFstream errorFile (errorFileName);

    errorFile << "# elapsed time | volume conservation errror | L1 advection error | "
        << "L1 advection error norm \n";

    forAll (elapsedTime_, I)
    {
        errorFile << elapsedTime_[I] << " " << errorVolume_[I] << " "
            << errorL1_[I] << " " << errorL1norm_[I] << "\n";
    }

    return true;
}

void fieldErrorFunctionObject::updateMesh(const mapPolyMesh&)
{
    notImplemented("fieldErrorFunctionObject::updateMesh(const mapPolyMesh&)");
}

void fieldErrorFunctionObject::movePoints(const polyMesh&)
{
    notImplemented("fieldErrorFunctionObject::movePoints(const polyMesh&)");
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
