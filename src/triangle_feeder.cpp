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
 * @author Nicolas Fortin
 * This project is the production of IFSTTAR (www.ifsttar.fr)
 * @copyright GNU Public License.V3
 * Official repository is https://github.com/nicolas-f/FastVoxel
 */
#include "triangle_feeder.hpp"
#include <tools/octree44_triangleElement.hpp>
#include <input_output/ply/rply_interface.hpp>
#include <cstring>
#include <limits>
#include <iomanip>
#include <iostream>

#ifndef MINREF
	#define MINREF(a, b)  if(a>b) a=b;
#endif
#ifndef MAXREF
	#define MAXREF(a, b)  if(a<b) a=b;
#endif
#ifndef MINVEC
	#define MINVEC(av, bv)  MINREF(av[0],bv[0]);MINREF(av[1],bv[1]);MINREF(av[2],bv[2]);
#endif
#ifndef MAXVEC
	#define MAXVEC(av, bv)  MAXREF(av[0],bv[0]);MAXREF(av[1],bv[1]);MAXREF(av[2],bv[2]);
#endif
#ifdef __USE_MULTITHREAD__
	#include <boost/thread/thread.hpp>
	#include <boost/bind.hpp>
#endif

#ifndef MIN
	#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif
#ifndef MAX
	#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif

/**
 * Calcul la liste des coordonn�es de cubes en intersection avec ce triangle
 * @param[in] boxCellCount Nombre de cellules dans les axes x,y,z
 * @param[in] boxCenter Position x,y,z de la boite englobante
 * @param[in] cellSize Taille d'une cellule
 * @param[in] triA Position du sommet A du triangle
 * @param[in] triB Position du sommet B du triangle
 * @param[in] triC Position du sommet C du triangle
 * @param[out] minRange Indice de d�but de l'intervalle
 * @param[out] maxRange Indice de fin de l'intervalle
 */

void GetRangeIntersectedBoundingCubeByTri(
    const ivec3& boxCellCount,
    const dvec3& zeroCellCenter,
    double cellSize,
    const dvec3& triA,
    const dvec3& triB,
    const dvec3& triC,
    ivec3& minRange,
    ivec3& maxRange)
{
    const dvec3 bmin(
        MIN(MIN(triA.x, triB.x), triC.x),
        MIN(MIN(triA.y, triB.y), triC.y),
        MIN(MIN(triA.z, triB.z), triC.z));

    const dvec3 bmax(
        MAX(MAX(triA.x, triB.x), triC.x),
        MAX(MAX(triA.y, triB.y), triC.y),
        MAX(MAX(triA.z, triB.z), triC.z));

    const dvec3 halfCell(cellSize / 2.0,
                         cellSize / 2.0,
                         cellSize / 2.0);

    const dvec3 gridMin = zeroCellCenter - halfCell;

    const dvec3 minOffset = (bmin - gridMin) / cellSize;
    const dvec3 maxOffset = (bmax - gridMin) / cellSize;

    // Include the voxel preceding the minimum boundary.
    // The SAT test removes any additional non-intersecting candidates.
    const ivec3 rawMin(
        static_cast<int>(floor(minOffset.x)) - 1,
        static_cast<int>(floor(minOffset.y)) - 1,
        static_cast<int>(floor(minOffset.z)) - 1);

    const ivec3 rawMax(
        static_cast<int>(floor(maxOffset.x)),
        static_cast<int>(floor(maxOffset.y)),
        static_cast<int>(floor(maxOffset.z)));

    minRange = ivec3(
        MAX(0, rawMin.x),
        MAX(0, rawMin.y),
        MAX(0, rawMin.z));

    maxRange = ivec3(
        MIN(boxCellCount.x - 1, rawMax.x),
        MIN(boxCellCount.y - 1, rawMax.y),
        MIN(boxCellCount.z - 1, rawMax.z));
}

namespace ScalarFieldBuilders
{
	TriangleScalarFieldCreator::TriangleScalarFieldCreator(
        const decimal& _resolution)
	:ScalarFieldCreator(_resolution)
	{
	}
    bool TriangleScalarFieldCreator::LoadPlyModel(
        const std::string& fileInput)
    {
        formatRPLY::t_model model3D;
        if(!formatRPLY::CPly::ImportPly(model3D,fileInput) || model3D.modelVertices.size()==0)
            return false;

        dvec3 minBoundingBox = *model3D.modelVertices.begin();
        dvec3 maxBoundingBox = *model3D.modelVertices.begin();
        std::vector<dvec3> vertices_vec;
        vertices_vec.reserve(model3D.modelVertices.size());

        for (std::list<dvec3>::const_iterator itvert = model3D.modelVertices.begin();
             itvert != model3D.modelVertices.end();
             ++itvert)
        {
            MAXVEC(maxBoundingBox, (*itvert));
            MINVEC(minBoundingBox, (*itvert));
            vertices_vec.push_back(*itvert);
        }

        model3D.modelVertices.clear();
        this->FirstStep_Params(minBoundingBox,maxBoundingBox);
        std::cout << std::setprecision(17)
            << "\n--- Grid parameters ---\n"
            << "Dimensions: "
            << this->volumeInfo.cellCountXYZ.x << " "
            << this->volumeInfo.cellCountXYZ.y << " "
            << this->volumeInfo.cellCountXYZ.z << "\n"
            << "Cell size: "
            << this->volumeInfo.cellSize << "\n"
            << "Box min: "
            << this->volumeInfo.boxMin.x << " "
            << this->volumeInfo.boxMin.y << " "
            << this->volumeInfo.boxMin.z << "\n"
            << "Box max: "
            << this->volumeInfo.boxMax.x << " "
            << this->volumeInfo.boxMax.y << " "
            << this->volumeInfo.boxMax.z << "\n"
            << "Main volume center: "
            << this->volumeInfo.mainVolumeCenter.x << " "
            << this->volumeInfo.mainVolumeCenter.y << " "
            << this->volumeInfo.mainVolumeCenter.z << "\n"
            << "Zero cell center: "
            << this->volumeInfo.zeroCellCenter.x << " "
            << this->volumeInfo.zeroCellCenter.y << " "
            << this->volumeInfo.zeroCellCenter.z << "\n"
            << "-----------------------\n";
        std::list<std::size_t>::iterator itlayerindex=model3D.modelFacesLayerIndex.begin();
        std::size_t layerIndex=1;
        for(std::list<formatRPLY::t_face>::iterator itface=model3D.modelFaces.begin();
            itface!=model3D.modelFaces.end();
            itface++)
        {
            if(itlayerindex!=model3D.modelFacesLayerIndex.end())
            {
                layerIndex=*itlayerindex;
                itlayerindex++;
            }
            this->SecondStep_PushTri(vertices_vec[(*itface).indicesSommets.a],
                vertices_vec[(*itface).indicesSommets.b],
                vertices_vec[(*itface).indicesSommets.c],
                layerIndex);
        }
        this->ThirdStep_VolumesCreator();
		return true;
    }
    void TriangleScalarFieldCreator::SecondStep_PushTri(
        const dvec3& A,
        const dvec3& B,
        const dvec3& C,
        const SpatialDiscretization::weight_t& marker)
	{
		#ifdef _DEBUG
		bool insideABox(false);
		#endif

        // Reject degenerate or numerically unstable triangles.
        // The degeneracy threshold remains based on float precision.
        const dvec3 edgeAB = B - A;
        const dvec3 edgeAC = C - A;
        const dvec3 edgeBC = C - B;

        dvec3 normal;
        normal.cross(edgeAB, edgeAC);

        const double edgeAB2 = edgeAB.dot(edgeAB);
        const double edgeAC2 = edgeAC.dot(edgeAC);
        const double edgeBC2 = edgeBC.dot(edgeBC);

        const double maxEdgeSquared = MAX(
            MAX(edgeAB2, edgeAC2),
            edgeBC2
        );

        const double normalSquared = normal.dot(normal);

        const double floatEpsilon =
            static_cast<double>(std::numeric_limits<float>::epsilon());

        if (maxEdgeSquared == 0.0 ||
            normalSquared <=
                floatEpsilon * floatEpsilon
                * maxEdgeSquared * maxEdgeSquared)
        {
            return;
        }

        this->volumeInfo.maximal_marker_index=MAX(this->volumeInfo.maximal_marker_index,marker);
		using namespace SpatialDiscretization;
		ivec3 minRange,maxRange;
		GetRangeIntersectedBoundingCubeByTri(
            this->volumeInfo.cellCountXYZ,
            this->volumeInfo.zeroCellCenter,
            this->volumeInfo.cellSize,
            A, B, C, minRange, maxRange);

        // The triangle does not intersect any voxel in the grid.
        if (minRange.x > maxRange.x
            || minRange.y > maxRange.y
            || minRange.z > maxRange.z)
        {
            #ifdef _DEBUG
            inoutbox.b++;
            #endif
            return;
        }
        double boxhalfsize[3];

        boxhalfsize[0] = this->volumeInfo.cellHalfSize.x;
        boxhalfsize[1] = this->volumeInfo.cellHalfSize.y;
        boxhalfsize[2] = this->volumeInfo.cellHalfSize.z;

        double triverts[3][3];

        // Todo multi-thread sur X ou X,Y
        double boxcenter[3] = {0.0, 0.0, 0.0};
        dvec3 boxcenter_double;

        for (int cell_x = minRange.x;
             cell_x <= maxRange.x;
             ++cell_x)
        {
            for (int cell_y = minRange.y;
                cell_y <= maxRange.y;
                ++cell_y)
            {
                for (int cell_z = minRange.z;
                    cell_z <= maxRange.z;
                    ++cell_z)
                {
                    boxcenter_double = CellIdToCenterCoordinate(
                        ivec3(cell_x, cell_y, cell_z),
                        this->volumeInfo.cellSize,
                        this->volumeInfo.zeroCellCenter
                    );

                    // Compute vertex coordinates relative to the box center in double precision.
                    triverts[0][0] = A.x - boxcenter_double.x;
                    triverts[0][1] = A.y - boxcenter_double.y;
                    triverts[0][2] = A.z - boxcenter_double.z;

                    triverts[1][0] = B.x - boxcenter_double.x;
                    triverts[1][1] = B.y - boxcenter_double.y;
                    triverts[1][2] = B.z - boxcenter_double.z;

                    triverts[2][0] = C.x - boxcenter_double.x;
                    triverts[2][1] = C.y - boxcenter_double.y;
                    triverts[2][2] = C.z - boxcenter_double.z;

                    if (boxtri_test::triBoxOverlap(
                            boxcenter,
                            boxhalfsize,
                            triverts) == 1)
                            {
						//TODO optimize to set a range of Z to the corresponding values
						(*this->fieldData)[cell_x][cell_y].SetData(cell_z,this->domainInformation,weight_t(marker));

						#ifdef _DEBUG
						insideABox=true;
						//Check Z length
						cell_id_t cell_z_test=0;
						zcell* currentCell=&((*this->fieldData)[cell_x][cell_y]);
						while(currentCell)
						{
							cell_z_test+=currentCell->GetSize();
							currentCell->Next(&currentCell);
						}
						if (cell_z_test != volumeInfo.cellCountXYZ.z)
							throw "error z length";
						#endif
					}
				}
			}
		}
		#ifdef _DEBUG
		if(!insideABox)
			inoutbox.b++;
		else
			inoutbox.a++;
		#endif
	}
}
