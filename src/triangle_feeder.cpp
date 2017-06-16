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
void GetRangeIntersectedBoundingCubeByTri(const SpatialDiscretization::cell_id_t& boxCellCount,const vec3& boxCenter,const decimal& cellSize,const vec3& triA,const vec3& triB,const vec3& triC, ivec3& minRange,ivec3& maxRange)
{
	//Ancienne m�thode
	vec3 bmin(MIN(MIN(triA.x,triB.x),triC.x),MIN(MIN(triA.y,triB.y),triC.y),MIN(MIN(triA.z,triB.z),triC.z));
	vec3 bmax(MAX(MAX(triA.x,triB.x),triC.x),MAX(MAX(triA.y,triB.y),triC.y),MAX(MAX(triA.z,triB.z),triC.z));
	vec3 tmpvec=((bmin-boxCenter)/cellSize);
	ivec3 halfCellCount(boxCellCount/2,boxCellCount/2,boxCellCount/2);
	minRange=ivec3((long)floor(tmpvec.x),(long)floor(tmpvec.y),(long)floor(tmpvec.z))+halfCellCount;
	tmpvec=((bmax-boxCenter)/cellSize);
	maxRange=ivec3((long)ceil(tmpvec.x),(long)ceil(tmpvec.y),(long)ceil(tmpvec.z))+halfCellCount;
}

namespace ScalarFieldBuilders
{

	TriangleScalarFieldCreator::TriangleScalarFieldCreator(const decimal& _resolution)
	:ScalarFieldCreator(_resolution)
	{


	}

    bool TriangleScalarFieldCreator::LoadPlyModel(const std::string& fileInput)
    {
        formatRPLY::t_model model3D;
        if(!formatRPLY::CPly::ImportPly(model3D,fileInput) || model3D.modelVertices.size()==0)
            return false;

        vec3 minBoundingBox=(*model3D.modelVertices.begin());
        vec3 maxBoundingBox=(*model3D.modelVertices.begin());
        std::vector<vec3> vertices_vec;
        vertices_vec.reserve(model3D.modelVertices.size());
        for(std::list<vec3>::iterator itvert=model3D.modelVertices.begin();itvert!=model3D.modelVertices.end();itvert++)
        {
            MAXVEC(maxBoundingBox,(*itvert));
            MINVEC(minBoundingBox,(*itvert));
            vertices_vec.push_back(*itvert);
        }

        model3D.modelVertices.clear();
        this->FirstStep_Params(minBoundingBox,maxBoundingBox);
        std::list<std::size_t>::iterator itlayerindex=model3D.modelFacesLayerIndex.begin();
        std::size_t layerIndex=1;
        for(std::list<formatRPLY::t_face>::iterator itface=model3D.modelFaces.begin();itface!=model3D.modelFaces.end();itface++)
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
	void TriangleScalarFieldCreator::SecondStep_PushTri(const vec3& A,const vec3& B,const vec3& C,const SpatialDiscretization::weight_t& marker)
	{
		#ifdef _DEBUG
		bool insideABox(false);
		#endif
        this->volumeInfo.maximal_marker_index=MAX(this->volumeInfo.maximal_marker_index,marker);
		using namespace SpatialDiscretization;
		ivec3 minRange,maxRange;
		GetRangeIntersectedBoundingCubeByTri(this->volumeInfo.cellCount,this->volumeInfo.mainVolumeCenter,this->volumeInfo.cellSize,A,B,C,minRange,maxRange);
		decimal boxhalfsize[3],triverts[3][3];

		memcpy(boxhalfsize,&this->volumeInfo.cellHalfSize,sizeof(vec3));
		memcpy(triverts[0],&A,sizeof(vec3));
		memcpy(triverts[1],&B,sizeof(vec3));
		memcpy(triverts[2],&C,sizeof(vec3));

		//Todo multi-thread sur X ou X,Y
		vec3 boxcenter;
		for(cell_id_t cell_x=minRange.x;cell_x<=(cell_id_t)maxRange.x;cell_x++)
		{
			for(cell_id_t cell_y=minRange.y;cell_y<=(cell_id_t)maxRange.y;cell_y++)
			{
				for(cell_id_t cell_z=minRange.z;cell_z<=(cell_id_t)maxRange.z;cell_z++)
				{
					boxcenter=CellIdToCenterCoordinate(ivec3(cell_x,cell_y,cell_z),this->volumeInfo.cellSize,this->volumeInfo.zeroCellCenter);
					if(boxtri_test::triBoxOverlap(boxcenter,boxhalfsize,triverts)==1)
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
						if(cell_z_test!=volumeInfo.cellCount)
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
