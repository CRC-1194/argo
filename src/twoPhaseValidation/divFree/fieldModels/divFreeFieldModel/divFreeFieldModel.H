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


#ifndef divFreeFieldModel_H
#define divFreeFieldModel_H

#include "stringList.H"
#include "Time.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         A templated functionObject
\*---------------------------------------------------------------------------*/

class divFreeFieldModel
{
        const Time& time_;

        // Used for field lookup when field references are not provided.
        word cellVelocityName_;

        word volFluxName_;

        word pointVelocityName_;

public:

    //- Runtime type information
    TypeName("divFreeFieldModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        divFreeFieldModel,
        dictionary,
        (const Time& time, const dictionary& dict),
        (time, dict)
    )

    declareRunTimeSelectionTable
    (
        autoPtr,
        divFreeFieldModel,
        time,
        (const Time& time),
        (time)
    )

    // Constructors
    divFreeFieldModel
    (
        const Time& time, 
        const dictionary& dict
    );

    // Selectors
    
    //- Select with a name 
    static autoPtr<divFreeFieldModel> New
    (
        const word& name, 
        const Time& time,
        const dictionary& dict

    );

    //- Select with a name from dictionary
    static autoPtr<divFreeFieldModel> New
    (
        const Time& time,
        const dictionary& dict

    );

    //- Select from controlDict::functions::divFree 
    static autoPtr<divFreeFieldModel> New
    (
        const Time& time
    );

    virtual ~divFreeFieldModel() = default;

    // Member Functions
    
    // Set the cell centered velocity field.
    void setCellVelocity(const Time& time);

    // Set the volumetric flux field.
    void setFlux(const Time& time);

    // Set the point velocity field.
    void setPointVelocity(const Time& time);

    // Return the point velocity field.
    void setPointVelocity(pointVectorField& pf);

    // Compute velocity at point X in time t.
    virtual vector velocity(point X, scalar t) const;

    virtual bool execute(); 

    //- Access time
    inline const Time& getTime() const
    {
        return time_;
    };

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
