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
#include "scalar_field_creator.hpp"
#include <iostream>
#include <fstream>
#include <list>
#include <utility>
#include <cstring>
#include <stdexcept>
#include <input_output/progressionInfo.h>
namespace ScalarFieldBuilders
{
    inline unsigned int cell_addr(int i,int j,int k,int nj,int nk)
    {
        return k+j*nk+i*nk*nj;
    }
	inline SpatialDiscretization::cell_id_t At(const SpatialDiscretization::cell_id_t& X,const SpatialDiscretization::cell_id_t& Y, const SpatialDiscretization::cell_id_t& Size)
	{
		return X+Y*Size;
	}

	#ifndef MIN
		#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
	#endif
	#ifndef MAX
		#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
	#endif

    dvec3 CellIdToCenterCoordinate( const ivec3& cell_id, const double_t& cellSize, const dvec3& zeroCellCenter)
	{
        return  zeroCellCenter+dvec3(cellSize*cell_id.x,cellSize*cell_id.y,cellSize*cell_id.z);
	}

    ScalarFieldCreator::ScalarFieldCreator(const double_t& _resolution)
		:resolution(_resolution)
	{



	}
	ScalarFieldCreator::~ScalarFieldCreator()
	{
	}
    void ScalarFieldCreator::ComputeMatrixParams(const dvec3& boxMin,const dvec3& boxMax, const double_t& minResolution, mainVolumeConstruction_t& computedVolumeInfo)
	{
    const dvec3 boxsize = boxMax - boxMin;

    if (minResolution <= 0.)
        throw std::invalid_argument("minResolution must be positive");

    const ivec3 cellCount(
        static_cast<int>(ceil(boxsize.x / minResolution)),
        static_cast<int>(ceil(boxsize.y / minResolution)),
        static_cast<int>(ceil(boxsize.z / minResolution)));

    computedVolumeInfo.boxMin = boxMin;
    computedVolumeInfo.boxMax = boxMax;
    computedVolumeInfo.cellCountXYZ = cellCount;
    computedVolumeInfo.cellCount =
        MAX(MAX(cellCount.x, cellCount.y), cellCount.z);

    computedVolumeInfo.mainVolumeCenter = (boxMax + boxMin) / 2.0;

    // Keep the requested voxel size; do not shrink or enlarge it.
    computedVolumeInfo.cellSize = minResolution;
    computedVolumeInfo.cellHalfSize =
        dvec3(minResolution, minResolution, minResolution) / 2.0;

    // Center the complete voxel grid on the requested bounding box.
    // zeroCellCenter is the center of voxel (0, 0, 0).
    computedVolumeInfo.zeroCellCenter = computedVolumeInfo.mainVolumeCenter
										- dvec3((cellCount.x / 2) * minResolution,
												(cellCount.y / 2) * minResolution,
												(cellCount.z / 2) * minResolution);

    computedVolumeInfo.maximal_marker_index = 0;
}
    void ScalarFieldCreator::FirstStep_Params(const dvec3& boxMin,const dvec3& boxMax)
	{
		//Compute the bouding boxes
        dvec3 cellCubeSize(resolution,resolution,resolution);
		cellCubeSize*=2;
		ComputeMatrixParams(boxMin-cellCubeSize,boxMax+cellCubeSize,resolution,this->volumeInfo);
		//Allocate matrix

		domainInformation.domainSize=this->volumeInfo.cellCount;
		domainInformation.domainSizeX=this->volumeInfo.cellCountXYZ.x;
		domainInformation.domainSizeY=this->volumeInfo.cellCountXYZ.y;
		domainInformation.domainSizeZ=this->volumeInfo.cellCountXYZ.z;
		domainInformation.weight=0;
		fieldData=PTR<SpatialDiscretization::weight_matrix>(new SpatialDiscretization::weight_matrix(domainInformation));
	}
	unsigned int ScalarFieldCreator::count()
	{
		unsigned int cpt(0);
		(*fieldData).Count(cpt);
		return cpt;
	}

	std::size_t ScalarFieldCreator::GetDomainSize()
	{
		return this->volumeInfo.cellCount;
	}

	ivec3 ScalarFieldCreator::GetDomainShape()
	{
		return this->volumeInfo.cellCountXYZ;
	}
    void ScalarFieldCreator::GetMinMax(dvec3& minBox,dvec3& maxBox)
	{
		minBox=this->volumeInfo.boxMin;
		maxBox=this->volumeInfo.boxMax;
	}

    double_t ScalarFieldCreator::GetVolumeValue(const SpatialDiscretization::weight_t& volId)
	{
		if(volId>=0 && (std::size_t)volId<this->volumeInfo.volumeValue.size())
		{
            return this->volumeInfo.volumeValue[volId];
		}else{
			return -1.f;
		}
	}
	SpatialDiscretization::weight_t ScalarFieldCreator::GetFirstVolumeIndex()
	{
	    return this->volumeInfo.maximal_marker_index+1;
	}
	SpatialDiscretization::weight_t ScalarFieldCreator::GetVolumeCount()
	{
		return this->volumeInfo.volumeCount;
    }

	bool ScalarFieldCreator::IsContainsVol( const ivec2& xyCell, SpatialDiscretization::weight_t& volId)
	{
		using namespace SpatialDiscretization;
		zcell* currentCell=&((*this->fieldData)[xyCell.x][xyCell.y]);
		while(currentCell->Next(&currentCell))
		{
			if(currentCell->GetData()==volId)
				return true;
		}
		return false;
	}
	void ScalarFieldCreator::GetMinMaxOnZ( const ivec2& xyCell, SpatialDiscretization::weight_t& minVolId,SpatialDiscretization::weight_t& maxVolId)
	{
		using namespace SpatialDiscretization;
		zcell* currentCell=&((*this->fieldData)[xyCell.x][xyCell.y]);
		minVolId=currentCell->GetData();
		maxVolId=minVolId;
		while(currentCell->Next(&currentCell))
		{
			minVolId=MIN(minVolId,currentCell->GetData());
			maxVolId=MAX(maxVolId,currentCell->GetData());
		}
	}
    typedef std::pair<SpatialDiscretization::weight_t, double_t> listValue_t;
	int sortFunc(const listValue_t& left,const listValue_t& right)
	{
		return left.second>right.second;
	}
	SpatialDiscretization::weight_t ScalarFieldCreator::GetLargestVolumeId()
	{
		std::list<listValue_t> volumeList;
		for(int i=0;i<GetVolumeCount();i++)
		{
			volumeList.push_back(listValue_t(i,GetVolumeValue(i)));
		}
		volumeList.sort(&sortFunc);
		return volumeList.front().first;
	}
	void ScalarFieldCreator::ExportVolsStats(const std::string& fileName, const std::string& volsLabelsFileName)
	{
		std::ofstream statVolsFile;
		statVolsFile.open(fileName.c_str(),std::ios_base::out);
		if(!statVolsFile)
			return;

		statVolsFile<<"Volume id"<<";value(m^3)"<<std::endl;
		//On va trier par ordre d�croissant de volume
		std::list<listValue_t> volumeList;
		for(int i=0;i<GetVolumeCount();i++)
		{
			volumeList.push_back(listValue_t(i,GetVolumeValue(i)));
		}
		volumeList.sort(&sortFunc);
		for(std::list<listValue_t>::iterator itvol=volumeList.begin();itvol!=volumeList.end();itvol++)
		{
			statVolsFile<<itvol->first<<";"<<itvol->second<<std::endl;
		}
		statVolsFile.close();

	}
	std::size_t GetCellFromAdress(const std::size_t& x,const std::size_t& y,const std::size_t& z,const std::size_t& sizex,const std::size_t& sizey,const std::size_t& sizez)
	{
		//return x+y*sizex+z*sizex*sizey;
		return x*sizey*sizez+y*sizez+z;
	}
    ivec3 GetAdressFromCellId(const std::size_t& cellid,const std::size_t& sizex,const std::size_t& sizey,const std::size_t& sizez)
	{
		ivec3 adress;
		adress[2]=cellid / (sizey*sizex);
		std::size_t tcellid = cellid % (sizey*sizex);
		adress[1]=tcellid / sizex;
		adress[0]=tcellid % sizex;
		return adress;
	}

	void ScalarFieldCreator::MakeXYZ(const std::string& filename,const SpatialDiscretization::weight_t& idVol)
	{
		std::ofstream xyzFile;
		xyzFile.open(filename.c_str(),std::ios_base::app);
		if(!xyzFile)
			return;

		using namespace SpatialDiscretization;
		cell_id_t cell_z=0;
		for(cell_id_t cell_x=0;cell_x<volumeInfo.cellCountXYZ.x;cell_x++)
		{
			for(cell_id_t cell_y=0;cell_y<volumeInfo.cellCountXYZ.y;cell_y++)
			{
				cell_z=0;
				zcell* currentCell=&((*this->fieldData)[cell_x][cell_y]);
				while(currentCell)
				{
					if((*currentCell).GetData()==idVol)
					{
						for(cell_id_t cell_z_offset=0;cell_z_offset<currentCell->GetSize();cell_z_offset++)
						{
                            dvec3 cellCenter=CellIdToCenterCoordinate(ivec3(cell_x,cell_y,cell_z+cell_z_offset),this->volumeInfo.cellSize,this->volumeInfo.zeroCellCenter);
							xyzFile<<cellCenter.x<<" "<<cellCenter.y<<" "<<cellCenter.z<<std::endl;
						}
					}
					cell_z+=currentCell->GetSize();
					currentCell->Next(&currentCell);
				}
			}
		}
		xyzFile.close();

	}
	void ScalarFieldCreator::ExportIJKData(const std::string& Infilename,const std::string& Outfilename)
	{
		std::ifstream coordinatesFile;
		std::ofstream ijkFile;
		ijkFile.open(Outfilename.c_str(),std::ios_base::out);
		coordinatesFile.open(Infilename.c_str(),std::ios_base::in);
		if(!coordinatesFile.is_open())
		{
			std::cerr<<"Input coordinates translation file not found."<<std::endl;
			return;
		}
		while (! coordinatesFile.eof() )
		{
		  std::string line;
		  getline (coordinatesFile,line);
		  //id X Y Z
		  int objectIndex;
          dvec3 coordinate;
          sscanf(line.c_str(),"%i %lf %lf %lf",&objectIndex,&(coordinate.x),&(coordinate.y),&(coordinate.z));
		  ivec3 ijk=this->GetCellIdByCoord(coordinate);
		  ijkFile<<objectIndex<<" "<<ijk.a<<" "<<ijk.b<<" "<<ijk.c<<" "<<this->GetMatrixValue(ijk)<<std::endl;
		}
		coordinatesFile.close();
		ijkFile.close();
	}

	bool ScalarFieldCreator::CheckDiscretisation()

	{
		using namespace SpatialDiscretization;
		cell_id_t cell_z=0;
		for(cell_id_t cell_x=0;cell_x<volumeInfo.cellCountXYZ.x;cell_x++)
		{
			for(cell_id_t cell_y=0;cell_y<volumeInfo.cellCountXYZ.y;cell_y++)
			{
				cell_z=0;
				zcell* currentCell=&((*this->fieldData)[cell_x][cell_y]);
				while(currentCell)
				{
					SpatialDiscretization::weight_t cell_type=(*currentCell).GetData();
					cell_z+=currentCell->GetSize();
					currentCell->Next(&currentCell);
				}
				if(cell_z!=volumeInfo.cellCountXYZ.z)
					return false;
			}
		}
		return true;
	}
    void ScalarFieldCreator::GetCellValueBoundaries(ivec3& min,ivec3& max,const SpatialDiscretization::weight_t& volid)
    {
		using namespace SpatialDiscretization;
		cell_id_t min_x=volumeInfo.cellCountXYZ.x;
        cell_id_t min_y=volumeInfo.cellCountXYZ.y;
        cell_id_t min_z=volumeInfo.cellCountXYZ.z;
        cell_id_t max_x=0;
        cell_id_t max_y=0;
        cell_id_t max_z=0;
        std::size_t cell_z;
        for(cell_id_t cell_x=0;cell_x<volumeInfo.cellCountXYZ.x;cell_x++)
		{
		    for(cell_id_t cell_y=0;cell_y<volumeInfo.cellCountXYZ.y;cell_y++)
			{
				cell_z=0;
				zcell* currentCell=&((*this->fieldData)[cell_x][cell_y]);
				while(currentCell)
				{
					const SpatialDiscretization::weight_t& cell_type(currentCell->GetData());
                    if (cell_type==volid)
                    {
					    min_x=MIN(min_x,cell_x);
                        min_y=MIN(min_y,cell_y);
					    min_z=MIN(min_z,cell_z);
					    max_x=MAX(max_x,cell_x+1);
                        max_y=MAX(max_y,cell_y+1);
					    max_z=MAX(max_z,cell_z+currentCell->GetSize());
                    }
					cell_z+=currentCell->GetSize();
					currentCell->Next(&currentCell);
				}
			}
		}
		min.set(min_x,min_y,min_z);
		max.set(max_x,max_y,max_z);
    }

	void ScalarFieldCreator::ExportVTK(const std::string& filename,const SpatialDiscretization::weight_t& idVol)
	{
		progressionInfo exportProgressionInformation(3);

		std::ofstream xyzFile;
		xyzFile.open(filename.c_str(),std::ios_base::out);
		if(!xyzFile.is_open())
			return;

		using namespace SpatialDiscretization;
		cell_id_t cell_z=0;

		// RECHERCHE DES EXTREMAS
		cell_id_t min_x=volumeInfo.cellCountXYZ.x;
        cell_id_t min_y=volumeInfo.cellCountXYZ.y;
        cell_id_t min_z=volumeInfo.cellCountXYZ.z;
        cell_id_t max_x=0;
        cell_id_t max_y=0;
        cell_id_t max_z=0;
		bool somethingToExport=false;
		std::cout<<"Establishing volume bounding box"<<std::endl;
        for(cell_id_t cell_x=0;cell_x<volumeInfo.cellCountXYZ.x;cell_x++)
		{
		    for(cell_id_t cell_y=0;cell_y<volumeInfo.cellCountXYZ.y;cell_y++)
			{
				cell_z=0;
				zcell* currentCell=&((*this->fieldData)[cell_x][cell_y]);
				while(currentCell)
				{
					SpatialDiscretization::weight_t cell_type=(*currentCell).GetData();
					for(cell_id_t cell_z_offset=0;cell_z_offset<currentCell->GetSize();cell_z_offset++)
					{
					  if ((idVol==-1 && cell_type<=this->volumeInfo.maximal_marker_index && cell_type!=-1) || (idVol!=-1 && cell_type==idVol))
					  {
					    somethingToExport=true;
					    min_x=MIN(min_x,cell_x);
                        min_y=MIN(min_y,cell_y);
					    min_z=MIN(min_z,cell_z);
					    max_x=MAX(max_x,cell_x);
                        max_y=MAX(max_y,cell_y);
					    max_z=MAX(max_z,cell_z+currentCell->GetSize()-1);
					  }
					}
					cell_z+=currentCell->GetSize();
					currentCell->Next(&currentCell);
				}
			}
		}
		exportProgressionInformation.GetMainOperation()->Next();
		if(somethingToExport)
		{
			if(idVol!=-1)
			{
				if(min_x>0)
					min_x--;
				if(min_y>0)
					min_y--;
				if(min_z>0)
					min_z--;
				if(max_x<volumeInfo.cellCountXYZ.x-1)
					max_x++;
				if(max_y<volumeInfo.cellCountXYZ.y-1)
					max_y++;
				if(max_z<volumeInfo.cellCountXYZ.z-1)
					max_z++;
			}

			std::size_t sizex(max_x-min_x+1),sizey(max_y-min_y+1),sizez(max_z-min_z+1);
    		exportProgressionInformation.GetMainOperation()->Next();

			std::cout<<"Write file."<<std::endl;
			//ECRITURE DANS LE FICHIER
			xyzFile<<"# vtk DataFile Version 3.0"<<std::endl;
			xyzFile<<"Exemple STRUCTURED_POINTS"<<std::endl;
			xyzFile<<"ASCII"<<std::endl;
			xyzFile<<"DATASET STRUCTURED_POINTS"<<std::endl;
			xyzFile<<"DIMENSIONS "<<(max_x-min_x+1)<<" "<<(max_y-min_y+1)<<" "<<(max_z-min_z+1)<<std::endl;
            dvec3 origin=this->GetCenterCellCoordinates(ivec3(min_x, min_y, min_z));
			xyzFile<<"ORIGIN "<<origin.x<<" "<<origin.y<<" "<<origin.z<<std::endl;
			xyzFile<<"SPACING "<<volumeInfo.cellSize<<" "<<volumeInfo.cellSize<<" "<<volumeInfo.cellSize<<std::endl;
			xyzFile<<"POINT_DATA "<<(max_x-min_x+1)*(max_y-min_y+1)*(max_z-min_z+1)<<std::endl;
			xyzFile<<"SCALARS MATERIAL float"<<std::endl;
			xyzFile<<"LOOKUP_TABLE default"<<std::endl;

			progressOperation curProgress(exportProgressionInformation.GetMainOperation(),sizex);

			for(cell_id_t cell_z=min_z;cell_z<=max_z;cell_z++)
			{
				curProgress.Next();
				exportProgressionInformation.OutputCurrentProgression();
				for(cell_id_t cell_y=min_y;cell_y<=max_y;cell_y++)
				{
					for(cell_id_t cell_x=min_x;cell_x<=max_x;cell_x++)
					{
						SpatialDiscretization::weight_t cell_type=this->GetMatrixValue(ivec3(cell_x,cell_y,cell_z));
                        xyzFile<< cell_type <<std::endl;
					}
				}
			}

		}else{
			std::cerr<<"Nothing to export with theses parameters !"<<std::endl;
		}

		xyzFile.close();

	}

	void ScalarFieldCreator::InitExteriorVolumeId()
	{
		using namespace SpatialDiscretization;
		//Pour chaque s�rie Z
		// La premi�re et derni�re s�rie appartiendra au volume 2
		cell_id_t cell_z=0;
		for(cell_id_t cell_x=0;cell_x<volumeInfo.cellCountXYZ.x;cell_x++)
		{
			for(cell_id_t cell_y=0;cell_y<volumeInfo.cellCountXYZ.y;cell_y++)
			{
				zcell* currentCell=&((*this->fieldData)[cell_x][cell_y]);

				if(currentCell->GetData()==emptyValue)
					currentCell->SetData(weight_t(volumeInfo.maximal_marker_index+1));
				do
				{
					if(!currentCell->IsNextCell() && currentCell->GetData()==emptyValue)
						currentCell->SetData(weight_t(volumeInfo.maximal_marker_index+1));
				}while(currentCell->Next(&currentCell));

			}
		}

	}
	bool ScalarFieldCreator::CellToCellVolumePropagation(const ivec2& destinationPropa,const ivec2& sourcePropa,const SpatialDiscretization::weight_t& volumeId)
	{
		using namespace SpatialDiscretization;
		bool modification(false);
		zcell* sourceCell=&((*this->fieldData)[sourcePropa.x][sourcePropa.y]);
		zcell* destinationCell=&((*this->fieldData)[destinationPropa.x][destinationPropa.y]);
		cell_id_t sourceZ(0),destinationZ(0);
		do
		{
			if(sourceCell->GetData()==volumeId)
			{
				//On a trouv� une s�rie de Z correspondant � la valeur � �tendre
				//On navigue jusqu'� la position de la source
				while(destinationZ+destinationCell->GetSize()<sourceZ) //while(destinationZ<sourceZ)
				{
					destinationZ+=destinationCell->GetSize();
					if(!destinationCell->Next(&destinationCell))
						return modification;
				}
				//On affecte � partir de sourceZ jusqu'a sourceZ+sourceCell->GetSize() les noeuds dont la valeur est 0
				//while(destinationZ>= sourceZ || destinationZ+destinationCell->GetSize()<=sourceZ+sourceCell->GetSize())
				while(true)
				{
					if(destinationCell->GetData()==emptyValue)
					{
						destinationCell->SetData(volumeId);
						modification=true;
					}
					//On avance si notre pointer est encore dans la cellule source
					if(sourceZ+sourceCell->GetSize()>destinationZ+destinationCell->GetSize())	{

						destinationZ+=destinationCell->GetSize();
						if(!destinationCell->Next(&destinationCell))
							return modification;
					}else{
						break;
					}
				}
			}
			sourceZ+=sourceCell->GetSize();
		}while(sourceCell->Next(&sourceCell));
		return modification;
	}

    void ScalarFieldCreator::ComputeVolumesValue(std::vector<double_t>& volumeValue)
	{
		using namespace SpatialDiscretization;
        volumeValue=std::vector<double_t>(this->volumeInfo.volumeCount,0.);
        double_t cellVolume=pow(this->volumeInfo.cellSize, 3.);
		for(cell_id_t cell_x=0;cell_x<volumeInfo.cellCountXYZ.x;cell_x++)
		{
			for(cell_id_t cell_y=0;cell_y<volumeInfo.cellCountXYZ.y;cell_y++)
			{
				zcell* currentCell=&((*this->fieldData)[cell_x][cell_y]);
				do
				{
					SpatialDiscretization::weight_t cellValue((*currentCell).GetData()-this->volumeInfo.maximal_marker_index-1);
					if(cellValue>=0)
						volumeValue[cellValue]+=cellVolume*currentCell->GetSize();
				}while(currentCell->Next(&currentCell));
			}
		}
	}
	ivec3 ScalarFieldCreator::GetFirstCellByWeight(const SpatialDiscretization::weight_t& weight,SpatialDiscretization::zcell** foundCell )
	{
		using namespace SpatialDiscretization;
		cell_id_t cell_z;
		for(cell_id_t cell_x=0;cell_x<volumeInfo.cellCountXYZ.x;cell_x++)
		{
			for(cell_id_t cell_y=0;cell_y<volumeInfo.cellCountXYZ.y;cell_y++)
			{
				cell_z=0;
				zcell* currentCell=&((*this->fieldData)[cell_x][cell_y]);
				while(currentCell)
				{
					if((*currentCell).GetData()==weight)
					{
						*foundCell=currentCell;
						return ivec3(cell_x,cell_y,cell_z);
					}
					cell_z+=currentCell->GetSize();
					currentCell->Next(&currentCell);
				}
			}
		}
		*foundCell=NULL;
		return ivec3();
	}

    ivec3 ScalarFieldCreator::GetCellIdByCoord(const dvec3& position)
	{
        dvec3 tmpvec=((position-this->volumeInfo.mainVolumeCenter)/this->volumeInfo.cellSize);
		return ivec3(static_cast<int>(floor(tmpvec.x + this->volumeInfo.cellCountXYZ.x / 2.0)),static_cast<int>(floor(tmpvec.y + this->volumeInfo.cellCountXYZ.y / 2.0)),static_cast<int>(floor(tmpvec.z + this->volumeInfo.cellCountXYZ.z / 2.0)));
	}

    void ScalarFieldCreator::ExtandVolume(
        const SpatialDiscretization::weight_t& volumeId)
    {
        using namespace SpatialDiscretization;

        const long nx = this->volumeInfo.cellCountXYZ.x;
        const long ny = this->volumeInfo.cellCountXYZ.y;

        if (nx <= 0 || ny <= 0)
            return;

        const long neighLink[4][2] = {
            {0, 1}, {1, 0}, {0, -1}, {-1, 0}
        };

        const std::size_t sizeOfMatrixXY =
            static_cast<std::size_t>(nx) *
            static_cast<std::size_t>(ny);

        PTR_ARR<bool> cellsToCheck(new bool[sizeOfMatrixXY]);
        PTR_ARR<bool> nextCellsToCheck(new bool[sizeOfMatrixXY]);

        memset(cellsToCheck.get(), 1, sizeof(bool) * sizeOfMatrixXY);
        memset(nextCellsToCheck.get(), 0, sizeof(bool) * sizeOfMatrixXY);

        bool moreCellsToCheck = true;

        do
        {
            moreCellsToCheck = false;

            for (long cell_x = 0; cell_x < nx; ++cell_x)
            {
                for (long cell_y = 0; cell_y < ny; ++cell_y)
                {
                    const cell_id_t cellXY =
                        At(cell_x, cell_y, nx);

                    if (!cellsToCheck[cellXY])
                        continue;

                    const ivec2 currentid(
                        static_cast<int>(cell_x),
                        static_cast<int>(cell_y)
                    );

                    bool destCellUpdated = false;

                    for (unsigned short neigh = 0; neigh < 4; ++neigh)
                    {
                        const long fromX =
                            cell_x + neighLink[neigh][0];
                        const long fromY =
                            cell_y + neighLink[neigh][1];

                        if (fromX < 0 || fromX >= nx ||
                            fromY < 0 || fromY >= ny)
                            continue;

                        const ivec2 fromXY(
                            static_cast<int>(fromX),
                            static_cast<int>(fromY)
                        );

                        if (CellToCellVolumePropagation(
                                currentid, fromXY, volumeId))
                            destCellUpdated = true;
                    }

                    if (destCellUpdated)
                    {
                        moreCellsToCheck = true;

                        for (unsigned short neigh = 0; neigh < 4; ++neigh)
                        {
                            const long fromX =
                                cell_x + neighLink[neigh][0];
                            const long fromY =
                                cell_y + neighLink[neigh][1];

                            if (fromX < 0 || fromX >= nx ||
                                fromY < 0 || fromY >= ny)
                                continue;

                            const cell_id_t neighXY =
                                At(fromX, fromY, nx);

                            nextCellsToCheck[neighXY] = true;
                        }
                    }
                }
            }

            memcpy(
                cellsToCheck.get(),
                nextCellsToCheck.get(),
                sizeof(bool) * sizeOfMatrixXY
            );
            memset(
                nextCellsToCheck.get(),
                0,
                sizeof(bool) * sizeOfMatrixXY
            );

        } while (moreCellsToCheck);
    }

    dvec3 ScalarFieldCreator::GetCenterCellCoordinates( const ivec3& cell_id) const
	{
		return CellIdToCenterCoordinate(cell_id,this->volumeInfo.cellSize, this->volumeInfo.zeroCellCenter);
	}


    SpatialDiscretization::weight_t ScalarFieldCreator::GetMatrixValue(
        const ivec3 &index)
    {
        if (index.x >= 0 &&
            index.x < this->volumeInfo.cellCountXYZ.x &&
            index.y >= 0 &&
            index.y < this->volumeInfo.cellCountXYZ.y &&
            index.z >= 0 &&
            index.z < this->volumeInfo.cellCountXYZ.z)
        {
            return (*this->fieldData)[index.x][index.y][index.z];
        }

        throw std::out_of_range(
            "Requested index (x: " + std::to_string(index.x) +
            ", y: " + std::to_string(index.y) +
            ", z: " + std::to_string(index.z) +
            ") is out of bounds. Valid dimensions are: " +
            std::to_string(this->volumeInfo.cellCountXYZ.x) + " x " +
            std::to_string(this->volumeInfo.cellCountXYZ.y) + " x " +
            std::to_string(this->volumeInfo.cellCountXYZ.z)
        );
    }

	void ScalarFieldCreator::CopyMatrix(SpatialDiscretization::weight_t* data,int ni,int nj,int nk,const ivec3& extractPos)
	{
		const int sourceXBegin = extractPos.a < 0 ? 0 : extractPos.a;
		const int sourceYBegin = extractPos.b < 0 ? 0 : extractPos.b;
		const int sourceZBegin = extractPos.c < 0 ? 0 : extractPos.c;

		const int sourceXEnd = extractPos.a + ni < volumeInfo.cellCountXYZ.x ? extractPos.a + ni : volumeInfo.cellCountXYZ.x;
		const int sourceYEnd = extractPos.b + nj < volumeInfo.cellCountXYZ.y ? extractPos.b + nj : volumeInfo.cellCountXYZ.y;
		const int sourceZEnd = extractPos.c + nk < volumeInfo.cellCountXYZ.z ? extractPos.c + nk : volumeInfo.cellCountXYZ.z;

		if (sourceXBegin >= sourceXEnd || sourceYBegin >= sourceYEnd || sourceZBegin >= sourceZEnd)
		{
			return;
		}

		for (int sourceX = sourceXBegin;
			sourceX < sourceXEnd;
			++sourceX)
		{
			const int destX = sourceX - extractPos.a;

			for (int sourceY = sourceYBegin;
				sourceY < sourceYEnd;
				++sourceY)
			{
				const int destY = sourceY - extractPos.b;

				for (int sourceZ = sourceZBegin;
					sourceZ < sourceZEnd;
					++sourceZ)
				{
					const int destZ = sourceZ - extractPos.c;
					data[cell_addr(destX, destY, destZ, nj, nk)] = (*this->fieldData)[sourceX][sourceY][sourceZ];
				}
			}
		}
	}

	void ScalarFieldCreator::CopyMatrixFiltered(SpatialDiscretization::weight_t* data,int ni,int nj,int nk,const ivec3& extractPos,const SpatialDiscretization::weight_t* data_filter,int nindex)
	{
		const int sourceXBegin = extractPos.a < 0 ? 0 : extractPos.a;
		const int sourceYBegin = extractPos.b < 0 ? 0 : extractPos.b;
		const int sourceZBegin = extractPos.c < 0 ? 0 : extractPos.c;

		const int sourceXEnd = extractPos.a + ni < volumeInfo.cellCountXYZ.x ? extractPos.a + ni : volumeInfo.cellCountXYZ.x;
		const int sourceYEnd = extractPos.b + nj < volumeInfo.cellCountXYZ.y ? extractPos.b + nj : volumeInfo.cellCountXYZ.y;
		const int sourceZEnd = extractPos.c + nk < volumeInfo.cellCountXYZ.z ? extractPos.c + nk : volumeInfo.cellCountXYZ.z;

		if (sourceXBegin >= sourceXEnd || sourceYBegin >= sourceYEnd || sourceZBegin >= sourceZEnd)
		{
			return;
		}

		for (int sourceX = sourceXBegin;
			sourceX < sourceXEnd;
			++sourceX)
		{
			const int destX = sourceX - extractPos.a;
			for (int sourceY = sourceYBegin;
				sourceY < sourceYEnd;
				++sourceY)
			{
				const int destY = sourceY - extractPos.b;
				for (int sourceZ = sourceZBegin;
					sourceZ < sourceZEnd;
					++sourceZ)
				{
					const int destZ = sourceZ - extractPos.c;
					const SpatialDiscretization::weight_t cellData = (*this->fieldData)[sourceX][sourceY][sourceZ];
					if (cellData >= 0 && cellData < nindex)
					{
						data[cell_addr(destX, destY, destZ, nj, nk)] = data_filter[cellData];
					}
				}
			}
		}
	}

	void ScalarFieldCreator::ThirdStep_VolumesCreator()
	{
		using namespace SpatialDiscretization;
		//Initialisation du volume exterieur
		InitExteriorVolumeId();
		ExtandVolume(SpatialDiscretization::weight_t(this->volumeInfo.maximal_marker_index+1));
		zcell* foundCell(NULL);
		GetFirstCellByWeight(weight_t(SpatialDiscretization::emptyValue),&foundCell); //Find the first empty cell
		weight_t volId(this->volumeInfo.maximal_marker_index+2);
		while(foundCell!=NULL)
		{

			#ifdef _DEBUG
			std::cout<<"Propagation of volume id="<<volId<<std::endl;
			#endif
			//Initialisation du volume volId
			foundCell->SetData(volId);
			ExtandVolume(volId);
			//Passage au prochain volume

			volId++;
			GetFirstCellByWeight(weight_t(SpatialDiscretization::emptyValue),&foundCell);
		}
		volumeInfo.volumeCount=volId-this->volumeInfo.maximal_marker_index-1;
		ComputeVolumesValue(this->volumeInfo.volumeValue);
	}
}
