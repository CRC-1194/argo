/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Tomislav Maric, TU Darmstadt
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
    Global function templates for file output of geometrical data using the 
    legacy ASCII VTK format. 

\*---------------------------------------------------------------------------*/

#ifndef legacyVtkOutput_H
#define legacyVtkOutput_H

#include "foamVtkLegacyAsciiFormatter.H"
#include "foamVtkOutput.H"
#include <iomanip>
#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/**
 * @brief Write points as legacy VTK polydata.
 *
 * @tparam Points: Iterable collection of points. 
 * @param vtks: VTK stream
 * @param points: Points to be written. 
 * @param fileName: Name of the file. 
 * @param index: Optional index for file sorting in ParaView.
 */
template<typename Points> 
void write_vtk_points(
    std::ostream& vtks,
    Points const& points, 
    const std::string& fileName
)
{
    // Initialize the legacy ascii VTK formatter: 15 decimal point precision.
    Foam::vtk::legacyAsciiFormatter legacyFormat(vtks, 15); 

    // Write the file header based on the chosen format (legacy VTK).
    Foam::vtk::legacy::fileHeader(legacyFormat, fileName, Foam::vtk::fileTag::POLY_DATA);

    // Write point header
    Foam::vtk::legacy::beginPoints(vtks, points.size());

    // Write point coordinates. 
    for(const auto& point : points) 
    {
        for(char cI = 0; cI < 3; ++cI)
            vtks << point[cI] << " "; 
        vtks << "\n"; 
    }
    vtks << "\nPOINT_DATA " << points.size() << "\n";
}

/**
 * @brief  Appends points to the VTK stream and returns the stream.
 *
 * @tparam Stream
 * @tparam Points
 * @param fileName
 * @param points
 * @param fileIndex
 *
 * @return 
 */
template<typename Stream, typename Points> 
std::enable_if_t<std::is_base_of_v<std::ostream, Stream>, Stream>
vtk_stream( 
    std::string fileName,
    Points const& points, 
    unsigned int fileIndex = 0
)
{
    std::stringstream fileNameStream; 
    fileNameStream << fileName << "-" 
        << std::setw(9) << std::setfill('0') <<  fileIndex;  

    std::ofstream vtks (fileNameStream.str() + ".vtk"); 
    write_vtk_points(vtks, points, fileName);
    return vtks;
}

/**
 * @brief Appends scalar values to the legacy VTK stream.
 *
 * @tparam Stream
 * @tparam Scalars
 * @param vtks
 * @param values
 * @param fieldName
 * @param index
 *
 * @return 
 */
template<typename Stream, typename Scalars> 
std::enable_if_t<std::is_arithmetic_v<typename Scalars::value_type>>
write_to_vtk_stream(
    Stream& vtks,
    Scalars const& values, 
    const std::string& fieldName
)
{
    // Format for legacy VTK scalar data output:
    // SCALARS volume_scalars float nCmpts
    // LOOKUP_TABLE default
    vtks << "SCALARS " << fieldName << " float " << "1\n"; 
    vtks << "LOOKUP_TABLE default\n";

    // Write point values.
    for (const auto val : values)
        vtks << val << "\n"; 
}


/**
 * @brief Appends vector values to the legacy VTK stream. 
 *
 * @tparam Stream
 * @tparam Vectors
 * @param vtks
 * @param vecs
 * @param fieldName
 *
 * @return 
 */
template<typename Stream> 
void write_to_vtk_stream(
    Stream& vtks,
    const std::vector<Foam::vector>& vecs, 
    const string& fieldName
)
{
    // Format for legacy VTK vector data output:
    // VECTORS pnorm float
    vtks << "VECTORS " << fieldName << " float " << "\n"; 

    // Write point values.
    for (const auto& vec : vecs)
    {
        for (const auto cmpt : vec)
            vtks << cmpt << " "; 
        vtks << "\n";
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
