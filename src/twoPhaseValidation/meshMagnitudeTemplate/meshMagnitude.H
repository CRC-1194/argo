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
    Structure template that gets the mesh magnitude field from a
    GeometricField.  Uses partial class template specialization.

Author
    Tomislav Maric maric@mma.tu-darmstadt.de 
    Mathematical Modeling and Analysis
    Technische Universität Darmstadt

Contributors

\*---------------------------------------------------------------------------*/

#ifndef meshMagnitude_H
#define meshMagnitude_H

#include "fvPatchField.H"
#include "fvsPatchField.H"
#include "volMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Structure meshMagnitude Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class PatchField, class GeoMesh>
struct meshMagnitude
{
    static const
    scalarField& get(const GeometricField<Type, PatchField, GeoMesh>& f);
};

template<typename Type>
struct meshMagnitude<Type, fvPatchField, volMesh>
{
    static const scalarField& get(const GeometricField<Type, fvPatchField, volMesh>& vf)
    {
        return vf.mesh().V().field();
    }
};

template<typename Type>
struct meshMagnitude<Type, fvsPatchField, surfaceMesh>
{
    static const scalarField& get(const GeometricField<Type, fvsPatchField, surfaceMesh>& sf)
    {
        return sf.mesh().magSf();
    }
};

/*---------------------------------------------------------------------------*\
                    Overloading To Simplify User Interface
\*---------------------------------------------------------------------------*/

template<typename Type, template <typename> class PatchField, typename GeoMesh>
const scalarField&
getMeshMagnitude(const GeometricField<Type, PatchField, GeoMesh> & f)
{
    return meshMagnitude<Type, PatchField, GeoMesh>::get(f);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
