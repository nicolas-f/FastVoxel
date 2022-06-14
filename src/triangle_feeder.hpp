/*
 *     This file is part of FastVoxel.
 *
 *     FastVoxel is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     FastVoxel is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *     along with FastVoxel.  If not, see <http://www.gnu.org/licenses/>.
 * FastVoxel is a voxelisation library of polygonal 3d model and do volumes identifications.
 * It is dedicated to finite element solvers
 * @author Nicolas Fortin , Judicaël Picaut judicael.picaut (home) ifsttar.fr
 * Official repository is https://github.com/nicolas-f/FastVoxel
 */

#include "scalar_field_creator.hpp"

#ifndef __TRIANGLE_SCALARFIELDBUILDERS__
#define __TRIANGLE_SCALARFIELDBUILDERS__

namespace ScalarFieldBuilders
{


class TriangleScalarFieldCreator : public ScalarFieldCreator
{
#ifdef _DEBUG
	ivec2 inoutbox;
#endif
public:
 TriangleScalarFieldCreator(const decimal& _resolution);

 /**
  * Append a triangle to the scalar field
  * @param A Coordinate of the vertex A
  * @param B Coordinate of the vertex B
  * @param C Coordinate of the vertex C
  * @param marker Marker of the triangle. [0-32768]
  */
 void SecondStep_PushTri(const dvec3& A,const dvec3& B,const dvec3& C,const SpatialDiscretization::weight_t& marker=1);
 bool LoadPlyModel(const std::string& fileInput);
};

}

#endif
