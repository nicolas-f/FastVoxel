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
#include <stdio.h>
#include <string.h>
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
        dvec3 boxsize=boxMax-boxMin;
		ivec3 cellCount=ivec3((long)ceil(boxsize.x/minResolution),(long)ceil(boxsize.y/minResolution),(long)ceil(boxsize.z/minResolution));

		computedVolumeInfo.boxMin=boxMin;
		computedVolumeInfo.boxMax=boxMax;
		computedVolumeInfo.cellCount=MAX(MAX(cellCount.x,cellCount.y),cellCount.z);
		computedVolumeInfo.mainVolumeCenter=(boxMax+boxMin)/2.f;
		computedVolumeInfo.cellSize=MAX(MAX(boxsize.x,boxsize.y),boxsize.z)/computedVolumeInfo.cellCount;
        computedVolumeInfo.cellHalfSize=dvec3(computedVolumeInfo.cellSize,computedVolumeInfo.cellSize,computedVolumeInfo.cellSize)/2;
        computedVolumeInfo.zeroCellCenter=computedVolumeInfo.mainVolumeCenter-dvec3(computedVolumeInfo.cellSize,computedVolumeInfo.cellSize,computedVolumeInfo.cellSize)*(decimal)(computedVolumeInfo.cellCount/2);
		computedVolumeInfo.maximal_marker_index=0;
	}
    void ScalarFieldCreator::FirstStep_Params(const dvec3& boxMin,const dvec3& boxMax)
	{
		//Compute the bouding boxes
        dvec3 cellCubeSize(resolution,resolution,resolution);
		cellCubeSize*=2;
		ComputeMatrixParams(boxMin-cellCubeSize,boxMax+cellCubeSize,resolution,this->volumeInfo);
		//Allocate matrix

		domainInformation.domainSize=this->volumeInfo.cellCount;
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
		//On va trier par ordre décroissant de volume
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
		for(cell_id_t cell_x=0;cell_x<volumeInfo.cellCount;cell_x++)
		{
			for(cell_id_t cell_y=0;cell_y<volumeInfo.cellCount;cell_y++)
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
 		for(cell_id_t cell_x=0;cell_x<volumeInfo.cellCount;cell_x++)
		{
		    for(cell_id_t cell_y=0;cell_y<volumeInfo.cellCount;cell_y++)
			{
				cell_z=0;
				zcell* currentCell=&((*this->fieldData)[cell_x][cell_y]);
				while(currentCell)
				{
					SpatialDiscretization::weight_t cell_type=(*currentCell).GetData();
					cell_z+=currentCell->GetSize();
					currentCell->Next(&currentCell);
				}
				if(cell_z!=volumeInfo.cellCount)
					return false;
			}
		}
		return true;
	}
    void ScalarFieldCreator::GetCellValueBoundaries(ivec3& min,ivec3& max,const SpatialDiscretization::weight_t& volid)
    {
		using namespace SpatialDiscretization;
		cell_id_t min_x=volumeInfo.cellCount;
        cell_id_t min_y=volumeInfo.cellCount;
        cell_id_t min_z=volumeInfo.cellCount;
        cell_id_t max_x=0;
        cell_id_t max_y=0;
        cell_id_t max_z=0;
        std::size_t cell_z;
 		for(cell_id_t cell_x=0;cell_x<volumeInfo.cellCount;cell_x++)
		{
		    for(cell_id_t cell_y=0;cell_y<volumeInfo.cellCount;cell_y++)
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
		cell_id_t min_x=volumeInfo.cellCount;
        cell_id_t min_y=volumeInfo.cellCount;
        cell_id_t min_z=volumeInfo.cellCount;
        cell_id_t max_x=0;
        cell_id_t max_y=0;
        cell_id_t max_z=0;
		bool somethingToExport=false;
		std::cout<<"Establishing volume bounding box"<<std::endl;
 		for(cell_id_t cell_x=0;cell_x<volumeInfo.cellCount;cell_x++)
		{
		    for(cell_id_t cell_y=0;cell_y<volumeInfo.cellCount;cell_y++)
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
				if(max_x<volumeInfo.cellCount-1)
					max_x++;
				if(max_y<volumeInfo.cellCount-1)
					max_y++;
				if(max_z<volumeInfo.cellCount-1)
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
		//Pour chaque série Z
		// La première et dernière série appartiendra au volume 2
		cell_id_t cell_z=0;
		for(cell_id_t cell_x=0;cell_x<volumeInfo.cellCount;cell_x++)
		{
			for(cell_id_t cell_y=0;cell_y<volumeInfo.cellCount;cell_y++)
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
				//On a trouvé une série de Z correspondant à la valeur à étendre
				//On navigue jusqu'à la position de la source
				while(destinationZ+destinationCell->GetSize()<sourceZ) //while(destinationZ<sourceZ)
				{
					destinationZ+=destinationCell->GetSize();
					if(!destinationCell->Next(&destinationCell))
						return modification;
				}
				//On affecte à partir de sourceZ jusqu'a sourceZ+sourceCell->GetSize() les noeuds dont la valeur est 0
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
		for(cell_id_t cell_x=0;cell_x<volumeInfo.cellCount;cell_x++)
		{
			for(cell_id_t cell_y=0;cell_y<volumeInfo.cellCount;cell_y++)
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
		for(cell_id_t cell_x=0;cell_x<volumeInfo.cellCount;cell_x++)
		{
			for(cell_id_t cell_y=0;cell_y<volumeInfo.cellCount;cell_y++)
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
		ivec3 halfCellCount(this->volumeInfo.cellCount/2,this->volumeInfo.cellCount/2,this->volumeInfo.cellCount/2);
		return ivec3((long)floor(tmpvec.x),(long)floor(tmpvec.y),(long)floor(tmpvec.z))+halfCellCount;
	}
	void ScalarFieldCreator::ExtandVolume(const SpatialDiscretization::weight_t& volumeId)
	{
		long neighLink[4][2]={{0,1},
							  {1,0},
		                      {0,-1},
			                  {-1,0} };

		using namespace SpatialDiscretization;
		std::size_t sizeOfMatrixXY(this->volumeInfo.cellCount*this->volumeInfo.cellCount);
		PTR_ARR<bool> cellsToCheck(new bool[sizeOfMatrixXY]); //Matrice X,Y de booléen indiquant les cellules à vérifier dans le cycle courant
		PTR_ARR<bool> NextCellsToCheck(new bool[sizeOfMatrixXY]); //Matrice X,Y de booléen indiquant les cellules à vérifier lors du prochain cycle.
		memset(cellsToCheck.get(),1,sizeof(bool)*sizeOfMatrixXY);
		memset(NextCellsToCheck.get(),0,sizeof(bool)*sizeOfMatrixXY);
		bool moreCellsToCheck(true); //Si sur un cycle de test aucune cellule ne s'est vu modifié alors le volume volumeId est complétement défini
		do
		{
			moreCellsToCheck=false;
			cell_id_t cellXY;
			for(cell_id_t cell_x=0;cell_x<volumeInfo.cellCount;cell_x++)
			{
				for(cell_id_t cell_y=0;cell_y<volumeInfo.cellCount;cell_y++)
				{

					cellXY=At(cell_x,cell_y,this->volumeInfo.cellCount);
					if(cellsToCheck[cellXY]==true)
					{

						ivec2 currentid(cell_x,cell_y);
						/////////////////////////////
						// En fonction des [2-4] voisins nous allons propager les identifiants de la série.
						bool destCellUpdated(false);
						for(unsigned short neigh=0;neigh<4;neigh++)
						{
							ivec2 FromXY=ivec2(neighLink[neigh])+currentid;
							if((unsigned int)FromXY.x>=0 && (unsigned int)FromXY.x<this->volumeInfo.cellCount && (unsigned int)FromXY.y>=0 && (unsigned int)FromXY.y<this->volumeInfo.cellCount)
							{
								cell_id_t neighXY=At(FromXY.x,FromXY.y,this->volumeInfo.cellCount);
								if(CellToCellVolumePropagation(currentid,FromXY,volumeId))
									destCellUpdated=true;
							}
						}

						if(destCellUpdated) //Cellule de destination modifié, toute les cellule voisines doivent être marqué comme à subir une propagation
						{
							moreCellsToCheck=true;
							for(unsigned short neigh=0;neigh<4;neigh++)
							{
								ivec2 FromXY=ivec2(neighLink[neigh])+currentid;
								if((unsigned int)FromXY.x>=0 && (unsigned int)FromXY.x<this->volumeInfo.cellCount && (unsigned int)FromXY.y>=0 && (unsigned int)FromXY.y<this->volumeInfo.cellCount)
								{
									cell_id_t neighXY=At(FromXY.x,FromXY.y,this->volumeInfo.cellCount);
									NextCellsToCheck[neighXY]=true;
								}
							}
						}
					}
				}
			}
			memcpy(cellsToCheck.get(),NextCellsToCheck.get(),sizeof(bool)*sizeOfMatrixXY);
			memset(NextCellsToCheck.get(),0,sizeof(bool)*sizeOfMatrixXY);
		}while(moreCellsToCheck);
	}

    dvec3 ScalarFieldCreator::GetCenterCellCoordinates( const ivec3& cell_id) const
	{
		return CellIdToCenterCoordinate(cell_id,this->volumeInfo.cellSize, this->volumeInfo.zeroCellCenter);
	}

	SpatialDiscretization::weight_t ScalarFieldCreator::GetMatrixValue(const ivec3& index)
	{
		return (*this->fieldData)[index.x][index.y][index.z];
	}
    void ScalarFieldCreator::CopyMatrix(SpatialDiscretization::weight_t* data,int ni,int nj,int nk,const ivec3& extractPos)
    {
        ivec3 extractPosEnd(MIN(ni+extractPos.a,volumeInfo.cellCount),MIN(nj+extractPos.b,volumeInfo.cellCount),MIN(nk+extractPos.c,volumeInfo.cellCount));
		using namespace SpatialDiscretization;
		cell_id_t cell_z=0;
		int i,j,k;
		for(cell_id_t cell_x=extractPos.a;cell_x<extractPosEnd.a;cell_x++)
		{
		    i=cell_x-extractPos.a;
			for(cell_id_t cell_y=extractPos.b;cell_y<extractPosEnd.b;cell_y++)
			{
			    j=cell_y-extractPos.b;
				cell_z=0;
				zcell* currentCell=&((*this->fieldData)[cell_x][cell_y]);
				while(currentCell)
				{
				    //Todo optimize with memset
				    const SpatialDiscretization::weight_t &cell_data(currentCell->GetData());
				    std::size_t cellsize=currentCell->GetSize();
				    if(cell_z+currentCell->GetSize()>extractPos.c)
				    {
                        for(cell_id_t cell_z_offset=0;cell_z_offset<cellsize;cell_z_offset++)
                        {
                            if(cell_z+cell_z_offset<extractPosEnd.c)
                            {
                                k=cell_z+cell_z_offset-extractPos.c;
                                if(k>=0)
                                    data[cell_addr(i,j,k,nj,nk)]=cell_data;
                            }else{
                                break;
                                break;
                            }
                        }
				    }
					cell_z+=currentCell->GetSize();
					currentCell->Next(&currentCell);
				}
			}
		}
    }
    void ScalarFieldCreator::CopyMatrixFiltered(SpatialDiscretization::weight_t* data,int ni,int nj,int nk,const ivec3& extractPos,const SpatialDiscretization::weight_t* data_filter,int nindex )
    {
        ivec3 extractPosEnd(MIN(ni+extractPos.a,volumeInfo.cellCount),MIN(nj+extractPos.b,volumeInfo.cellCount),MIN(nk+extractPos.c,volumeInfo.cellCount));
		using namespace SpatialDiscretization;
		cell_id_t cell_z=0;
		int i,j,k;
		for(cell_id_t cell_x=extractPos.a;cell_x<extractPosEnd.a;cell_x++)
		{
		    i=cell_x-extractPos.a;
			for(cell_id_t cell_y=extractPos.b;cell_y<extractPosEnd.b;cell_y++)
			{
			    j=cell_y-extractPos.b;
				cell_z=0;
				zcell* currentCell=&((*this->fieldData)[cell_x][cell_y]);
				while(currentCell)
				{
				    //Todo optimize with memset
				    const SpatialDiscretization::weight_t &cell_data(currentCell->GetData());
				    std::size_t cellsize=currentCell->GetSize();
				    if(cell_z+currentCell->GetSize()>extractPos.c)
				    {
                        for(cell_id_t cell_z_offset=0;cell_z_offset<cellsize;cell_z_offset++)
                        {
                            if(cell_z+cell_z_offset<extractPosEnd.c)
                            {
                                k=cell_z+cell_z_offset-extractPos.c;
                                if(k>=0 && cell_data<nindex)
                                    data[cell_addr(i,j,k,nj,nk)]=data_filter[cell_data];
                            }else{
                                break;
                                break;
                            }
                        }
				    }
					cell_z+=currentCell->GetSize();
					currentCell->Next(&currentCell);
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
