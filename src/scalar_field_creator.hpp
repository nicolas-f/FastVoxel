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
 * @author Nicolas Fortin , Judica�l Picaut judicael.picaut (home) ifsttar.fr
 * Official repository is https://github.com/nicolas-f/FastVoxel
 */

#include "spatial_discretization.hpp"
#include <vector>
#include <string>

#ifndef __SCALARFIELDBUILDERS__
#define __SCALARFIELDBUILDERS__

namespace ScalarFieldBuilders
{
	/**
	 * Retourne les coordonn�es du centre du cube correspondant � l'indice en param�tre
	 */
    dvec3 CellIdToCenterCoordinate( const ivec3& cell_id, const double_t& cellSize, const dvec3& zeroCellCenter);
	/**
	 * Cette classe permet de g�n�rer un espace discr�tis� en plusieurs volumes. C'est la premi�re �tape de la reconstruction du mod�le.
	 */
	class ScalarFieldCreator
	{
	protected:
		PTR<SpatialDiscretization::weight_matrix> fieldData; //Donn�es de la matrice X,Y,Z


		struct mainVolumeConstruction_t
		{
		    mainVolumeConstruction_t()
		    : cellSize(0.),cellCount(0),volumeCount(0),maximal_marker_index(0)
            {
            }
            double_t cellSize;
			unsigned int cellCount;
            dvec3 mainVolumeCenter;
            dvec3 cellHalfSize;
            dvec3 zeroCellCenter;
			SpatialDiscretization::weight_t volumeCount;
            std::vector<double_t> volumeValue;
            dvec3 boxMin;
            dvec3 boxMax;
			SpatialDiscretization::weight_t maximal_marker_index;
		} volumeInfo;
		SpatialDiscretization::domainInformation_t domainInformation;
        double_t resolution;
        static void ComputeMatrixParams(const dvec3& boxMin,const dvec3& boxMax, const double_t& minResolution, mainVolumeConstruction_t& computedVolumeInfo);
		/**
		 * Initialise les donn�es pour le volume ext�rieur
		 */
		void InitExteriorVolumeId();
		/**
		 * Propage les valeurs des identifiants de volume d'une cellule � l'autre
		 * Si la cellule cible est modifi� cette m�thode retourne vrai
		 */
		bool CellToCellVolumePropagation(const ivec2& destinationPropa,const ivec2& sourcePropa,const SpatialDiscretization::weight_t& volumeId);
		/**
		 * Propage un indice de volume dans toute la matrice
		 */
		void ExtandVolume(const SpatialDiscretization::weight_t& volumeId);
		/**
		 * Retourne la position de la premi�re cellule avec la valeur en param�tre
		 */
		ivec3 GetFirstCellByWeight(const SpatialDiscretization::weight_t& weight,SpatialDiscretization::zcell** foundCell);

		/**
		 * Calcul pour chaque volume sa valeur en m^3
		 * @param[out] volumeValue Un tableau de dimension �gale au nombre de volume dans le domaine. Dont la valeur est en m^3.
		 */
        void ComputeVolumesValue(std::vector<double_t>& volumeValue);


	public:
		/**
		 * Constructeur
		 * @param _resolution Dimension d'une cellule qui composera la matrice. Plus la r�solution est �lev�e plus le model g�n�r� sera proche du mod�le en entr�e et plus de triangles seront g�n�r�s.
		 */
        ScalarFieldCreator(const double_t& _resolution);
		/**
		 * Initialisation de la matrice selon la r�solution et la boite englobante pass� en param�tre.
		 * @param boxMin Coordonn�es minimale des objets qui alimenteront la matrice
		 * @param boxMax Coordonn�es maximale des objets qui alimenteront la matrice
		 */
        void FirstStep_Params(const dvec3& boxMin,const dvec3& boxMax);
		virtual ~ScalarFieldCreator();

		/**
		 * Une fois toutes les primitives renseign�es. Cette m�thode doit �tre appel�e afin de d�tecter les volumes d�limit� par les limites.
		 */
		void ThirdStep_VolumesCreator();

		/**
		 * Retourne la valeur de la matrice selon les indices des cellules
		 * @param index Entier positif d�signant le n� de cellule.
		 * @see GetDomainSize()
		 */
		SpatialDiscretization::weight_t GetMatrixValue(const ivec3& index);


        dvec3 GetCenterCellCoordinates( const ivec3& cell_id) const;


        /**
         * Get the first cell value associated to matrix, this index start from the highest material index of materials
         */
        SpatialDiscretization::weight_t GetFirstVolumeIndex();

		/**
		 * Pour toutes les valeurs de Z pour un x,y donn�. Retourne la valeur min,max.
		 */
		void GetMinMaxOnZ( const ivec2& xyCell, SpatialDiscretization::weight_t& minVolId,SpatialDiscretization::weight_t& maxVolId);
		bool IsContainsVol( const ivec2& xyCell, SpatialDiscretization::weight_t& volId);

        /**
         * Get the volume count
         */
		SpatialDiscretization::weight_t GetVolumeCount();
		/**
		 * Get the volume (m3) corresponding to volume ID [0-GetVolumeCount()[
		 */
        double_t GetVolumeValue(const SpatialDiscretization::weight_t& volId);

        void GetMinMax(dvec3& minBox,dvec3& maxBox);

		/**
		 * Exporte les indices et volumes des domaines
		 * @param fileName Nom et chemin du fichier de sortie
		 * @param volsLabelsFileName Nom et chemin du fichier d'entr�e des noms utilisateur de volume (fichier texte avec "x y z NomDuVolume" � chaque ligne)
		 */
		void ExportVolsStats(const std::string& fileName, const std::string& volsLabelsFileName=std::string());
		///////////////////////////////
		// Debug Functions
		unsigned int count();
		void MakeXYZ(const std::string& filename,const SpatialDiscretization::weight_t& idVol);
		void ExportIJKData(const std::string& Infilename,const std::string& Outfilename);
        void MakeXYZ(const std::string& filename,const double_t& minVol);
		void ExportVTK(const std::string& filename,const SpatialDiscretization::weight_t& idVol=-1);
		std::size_t GetDomainSize();
		bool CheckDiscretisation();
		SpatialDiscretization::weight_t GetLargestVolumeId();

        /**
         * Copy data from spatial discretisation to array in argument
         * @param data Previously allocated data
         * @param ni x maximum size
         * @param nj y maximum size
         * @param nk z maximum size
         * @param extractPos zero position extraction
         */
        void CopyMatrix(SpatialDiscretization::weight_t* data,int ni,int nj,int nk,const ivec3& extractPos);
        void CopyMatrixFiltered(SpatialDiscretization::weight_t* data,int ni,int nj,int nk,const ivec3& extractPos,const SpatialDiscretization::weight_t* data_filter,int nindex );
        /**
		 * Retourne l'indice de la cellule contenant le point pass� en param�tre
		 */
        ivec3 GetCellIdByCoord(const dvec3& position);

		void GetCellValueBoundaries(ivec3& min,ivec3& max,const SpatialDiscretization::weight_t& volid);
	};

}

#endif
