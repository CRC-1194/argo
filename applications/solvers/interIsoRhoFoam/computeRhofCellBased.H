// FIXME: This shouldn't be boilerplate, refactor. TM. 
// Compute rhoPhi using the old approach since it is valid for all bulk cells
rhof.oldTime() == rhof;
rhof == fvc::interpolate(rho);

// #include "SortableList.H"
// #include <assert.h>

// Correct faces which are part of an interface cell (0 < alpha < 1)
const auto& faceOwner = mesh.owner();
//const auto& faceNeighbour = mesh.neighbour();
const auto& cellList = mesh.cells();
const auto& phiIn = phi.primitiveField();
// const auto& meshVertices = mesh.points();
// const auto& meshFaces = mesh.faces();

const reconstructionSchemes& surf = advector.surf(); //Jun: reference to reconstructionSchemes
const DynamicField<label>& interfaceLabels = surf.interfaceLabels(); // Jun: list of interface cells

forAll(interfaceLabels, i) // only interface cells are operated
{
    const label celliLabel = interfaceLabels[i];
   
    if (mag(surf.normal()[celliLabel]) != 0) // cell faces area vector cannot be zero
    { 
        const point x0 = surf.centre()[celliLabel]; // cell interface centre
        vector n0 = -surf.normal()[celliLabel]; // cell interface unit normal
        n0 /= (mag(n0));

        const cell& celli = cellList[celliLabel];
        forAll(celli, fi)
        {
            const label fiLabel = celli[fi];
            const scalar& sfArea = mesh.magSf()[fiLabel];

            if (mesh.isInternalFace(fiLabel))
            {
                bool isDownwindFace = false;

                if (celliLabel == faceOwner[fiLabel])
                {
                    if (phiIn[fiLabel] >= 0)
                    {
                        isDownwindFace = true;
                    }
                }
                else
                {
                     if (phiIn[fiLabel] < 0)
                     {
                         isDownwindFace = true;
                     }
                }

                if (isDownwindFace)
                {

                     cutFace.calcSubFace(fiLabel, n0, x0);
                     scalar alphaFace = mag(cutFace.subFaceArea()) / sfArea;

//                     assert(alphaFace >= 0.0 && "Negative face area fraction");
//                     assert(alphaFace <= 1.0 && "Face area fraction greater than 1");

                     rhof[fiLabel] = alphaFace*rho1.value() + (1 - alphaFace)*rho2.value();
                     muf[fiLabel] = alphaFace*rho1.value()*nu1.value() +(1 - alphaFace)*rho2.value()*nu2.value();
                }
            }
            else
            {
                      // Jun: how to handle boundary patches?
            }

        }
    }
}

/*
*/
