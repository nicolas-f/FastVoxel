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
#ifdef _DEBUG
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
	#include <crtdbg.h>
#endif

#include "triangle_feeder.hpp"

#include <iostream>
#include <input_output/ply/rply_interface.hpp>
#include <vector>
#include <string>
#include <string.h>
#include <stdio.h>
/**
 * Cette méthode permet le tirage aléatoire d'un nombre décimal
 * Return a float random number
 * @return Decimal from 0 to 1
 */
inline decimal GetRandValue()
{
	return ((decimal)rand()) / ((decimal)RAND_MAX+1);
}
inline int Range(const int& from,const int& to)
{
	int res=(int)floor(GetRandValue()*(to-from)+from);
	if(res<to)
		return res;
	else
		return from;
}


#ifndef MIN
	#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif
#ifndef MAX
	#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif

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
void PrintUsage(int argc, char* argv[])
{
	std::cout<<"Usage: "<< argv[0] <<" [-precDECIMAL] [-depthINTEGER] [-v] [-volstatsFILENAME] [-ivINTEGER] -iFILENAME -oFILENAME"<<std::endl;
	std::cout<<" -prec : Absolute cell size."<<std::endl;
	std::cout<<" -depth : [5-10] Relative cell size, cell subdivision count will be 2^depth . Default 5."<<std::endl;
	std::cout<<" -v : Verbose mode. Give more details about remeshing."<<std::endl;
	std::cout<<" -i : PLY input filename."<<std::endl;
	std::cout<<" -t : Coordinate translation. For each line, translate x,y,z coordinates into the corresponding i,j,k and volume id."<<std::endl;
	std::cout<<" -o : Output filename.Do not write extension. The format can be read by ParaView"<<std::endl;
	std::cout<<" -iv : [0-Volume Count] Export specified volumes,materials into separate files. You can add this parameter multiple times.You can use the option -volstats to find volumes ID."<<std::endl;
	std::cout<<" -volstats : Discretise with prec or depth and write volumes statistics to specified file then stop execution."<<std::endl;
	std::cout<<std::endl;
}

/**
 * State of program arguments, feed by the main program
 */
struct VolumeSelectionInfo_t
{
	VolumeSelectionInfo_t()
		:minimalVol(0.f),selectionFilter(FILTER_BY_VOLUME)
	{}
	enum FILTER
	{
		FILTER_EXT_ONLY,
		FILTER_SKIP_EXT,
		FILTER_BY_VOLUME,
		FILTER_BY_ID
	} selectionFilter;
	decimal minimalVol;
	std::list<SpatialDiscretization::cell_id_t> extractedVolumes;

};
#define sscanf_s sscanf

/**
 * FastVoxel in console mode
 */
int MainRemesh(int argc, char* argv[])
{
	using namespace ScalarFieldBuilders;

	VolumeSelectionInfo_t volumeSelectionInfo;
	decimal precision(0.f);
	unsigned long d(0);

	std::string fileInput;
	std::string fileOutput;
	std::string volStatsOutput;
	std::string translationFileInput;
	unsigned int depth(0);     //Domain subdivision in 2^n in each dimension
	unsigned int iv_buffer(0); //Extracted volume temporary variable
	bool verbose(false);
	//Scan user arguments
	try
	{
		while (--argc)
		{
			if (sscanf_s(argv[argc], "-prec%g", &precision) == 1);
			else if (sscanf_s(argv[argc], "-depth%lu", &d) == 1);
			else if (sscanf_s(argv[argc], "-iv%u", &iv_buffer) == 1)
				volumeSelectionInfo.extractedVolumes.push_back(iv_buffer);
			else if (sscanf_s(argv[argc], "-minvol%g", &volumeSelectionInfo.minimalVol) == 1);
			else if (strncmp(argv[argc], "-volstats", 9) == 0)
			  volStatsOutput = std::string(argv[argc] + 9);
			else if (strncmp(argv[argc], "-v", 2) == 0) verbose=true;
			else if (strncmp(argv[argc], "-i", 2) == 0)
			  fileInput = std::string(argv[argc] + 2);
			else if (strncmp(argv[argc], "-t", 2) == 0)
			  translationFileInput = std::string(argv[argc] + 2);
			else if (strncmp(argv[argc], "-o", 2) == 0)
			  fileOutput = std::string(argv[argc] + 2);
			else {
			  PrintUsage(argc,argv);
			  return -2;
			}
		}
	}catch(...)
	{
	  PrintUsage(argc,argv);
	  return -2;
	}
	//Incomplete arguments
	if(fileInput.empty() || (fileOutput.empty() && (volStatsOutput.empty() && translationFileInput.empty())))
	{
	  PrintUsage(argc,argv);
	  return -2;
	}
	vec3 minBoundingBox,maxBoundingBox;


	//Init the bounding box of the model
	std::cout<<"Open "<< fileInput<<std::endl;
	formatRPLY::t_model model3D;
	if(!formatRPLY::CPly::ImportPly(model3D,fileInput) || model3D.modelVertices.size()==0)
		return -2;

	minBoundingBox=(*model3D.modelVertices.begin());
	maxBoundingBox=(*model3D.modelVertices.begin());
	std::vector<vec3> vertices_vec;
	vertices_vec.reserve(model3D.modelVertices.size());
	for(std::list<vec3>::iterator itvert=model3D.modelVertices.begin();itvert!=model3D.modelVertices.end();itvert++)
	{
		MAXVEC(maxBoundingBox,(*itvert));
		MINVEC(minBoundingBox,(*itvert));
		vertices_vec.push_back(*itvert);
	}
	model3D.modelVertices.clear();


    //init the cell size, corresponding to specified arguments and bounding box
	if(precision==0 && d==0)
		d=5;
	if(d!=0)
	{
		vec3 axesPrecision=(maxBoundingBox-minBoundingBox)/vec3((decimal)pow((decimal)2,(int)d),(decimal)pow((decimal)2,(int)d),(decimal)pow((decimal)2,(int)d));
		precision=MAX(axesPrecision.x,MAX(axesPrecision.z,axesPrecision.y));
	}
	std::cout<<"Discretization of space with a "<<precision<<" m cell width"<<std::endl;

	ScalarFieldBuilders::TriangleScalarFieldCreator FromTriangleRemesh(precision);
	FromTriangleRemesh.FirstStep_Params(minBoundingBox,maxBoundingBox);

	std::size_t domainSize(FromTriangleRemesh.GetDomainSize());
	if(verbose)
		std::cout<<"Matrix size "<<domainSize<<"x"<<domainSize<<"x"<<domainSize<<" = "<<pow((decimal)domainSize,3)<<" cells max("<<pow((decimal)domainSize,2)<<"min)"<<std::endl;


	/////////////////////////////////////////////////////////////
	//Voxelisation of surfaces

	unsigned int idtri(0);
	std::size_t triCount(model3D.modelFaces.size());
	int lastprogression(0),progression(0);
	std::cout<<"Feeding matrix "<<std::endl;
	std::list<std::size_t>::iterator itlayerindex=model3D.modelFacesLayerIndex.begin();
	for(std::list<formatRPLY::t_face>::iterator itface=model3D.modelFaces.begin();itface!=model3D.modelFaces.end();itface++)
	{
	    std::size_t layerIndex=1;
	    //Iteration
	    if(itlayerindex!=model3D.modelFacesLayerIndex.end())
		{
            layerIndex=*itlayerindex;
			itlayerindex++;
		}
		//Add tri in voxel
		FromTriangleRemesh.SecondStep_PushTri(vertices_vec[(*itface).indicesSommets.a],
			vertices_vec[(*itface).indicesSommets.b],
			vertices_vec[(*itface).indicesSommets.c],
			layerIndex);
		if(verbose)
		{
			idtri++;
			progression=int(((float)idtri/(float)triCount)*100);
			if(progression!=lastprogression)
			{
				std::cout<<"Feeding matrix "<<progression<<"% face"<<idtri<<"/"<<triCount<<std::endl;
				lastprogression=progression;
			}
		}
	}
	if(verbose)
	{
		std::cout<<"Finish feeding matrix with intersection nodes"<<std::endl;
		std::cout<<"There are "<<FromTriangleRemesh.count()<<" cells in the matrix"<<std::endl;
	}


	/////////////////////////////////////////////////////////////
	//Init empty cells with unique volumetric information, starting with the higher material integer + 1


	FromTriangleRemesh.ThirdStep_VolumesCreator();
	if(verbose)
	{
		std::cout<<"Volumes defined :"<<std::endl;
		for(int i=0;i<FromTriangleRemesh.GetVolumeCount();i++)
		{
			std::cout<<"Volume id:"<<i<<" value:"<<FromTriangleRemesh.GetVolumeValue(i)<<" m^3"<<std::endl;
		}
	}else{
		std::cout<<"Volumes defined.."<<std::endl;
	}
	//Validation
	FromTriangleRemesh.CheckDiscretisation();
	/////////////////////////////////////////////////////////////
	//Save result in a file

	//Coordinates export
	if(!translationFileInput.empty())
	{
		FromTriangleRemesh.ExportIJKData(translationFileInput,translationFileInput+".vol");
	}
	if(!(volStatsOutput.empty()))
	{
		FromTriangleRemesh.ExportVolsStats(volStatsOutput);
		return 0;
	}

	if(!fileOutput.empty())
	{
		if(volumeSelectionInfo.extractedVolumes.empty())
		{
			//No selection done
			//Export the main volumes

			FromTriangleRemesh.ExportVTK(fileOutput);
		}else{
			std::size_t volumeCountToProcess(volumeSelectionInfo.extractedVolumes.size());

			std::cout<<"There are "<<volumeCountToProcess<<" volumes to export.."<<std::endl;
			for(std::list<SpatialDiscretization::cell_id_t>::iterator itvol=volumeSelectionInfo.extractedVolumes.begin();itvol!=volumeSelectionInfo.extractedVolumes.end();itvol++)
			{
				SpatialDiscretization::weight_t i(*itvol);
				std::cout<<"Export Volume "<<i<<std::endl;
				FromTriangleRemesh.ExportVTK(fileOutput,(SpatialDiscretization::weight_t)i);
			}
		}
	}
	return 0;
}
int main(int argc, char* argv[])
{
	int ret(MainRemesh(argc,argv));

    #ifdef _DEBUG
		_CrtDumpMemoryLeaks(); //Affiche les fuites mémoires
	#endif

	return ret;
}
