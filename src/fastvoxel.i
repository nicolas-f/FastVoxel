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

%module fastvoxel
%{
#define SWIG_FILE_WITH_INIT
#include "triangle_feeder.hpp"
%}
%include "std_string.i"
%include "typemaps.i"
%include "numpy.i"
%init %{
    import_array();
%}




namespace core_mathlib
{

    %rename(dvec3) base_vec3<double_t>;

    class dvec3 {
    public:
        dvec3(void);
        dvec3(const double& _x,const double& _y,const double& _z);
        dvec3(const dvec3 &_v);
        int operator==(const dvec3 &_v);
        int operator!=(const dvec3 &_v);
        const dvec3 operator*(double _f)const;
        const dvec3 operator/(double _f) const;
        dvec3 operator/(const dvec3 &_v) const;
        double length(void) const;
        const dvec3 operator+(const dvec3 &_v) const;
        const dvec3 operator-() const;
        const dvec3 operator-(const dvec3 &_v) const;
        dvec3 &operator*=(double _f);
        dvec3 &operator/=(double _f);
        dvec3 &operator+=(const dvec3 &_v);
        dvec3 &operator-=(const dvec3 &_v);
        float operator*(const dvec3 &_v) const;
        bool barelyEqual(const dvec3 &_v) const;
        void set(double _x,double _y,double _z);
        void reset(void);
        double normalize(void);
        void cross(const dvec3 &v1,const dvec3 &v2);
        void cross(const dvec3 &v2);
        double cosinus(const dvec3 &ac);
        double dot(const dvec3 &v) const;
        bool compare(const dvec3 &_v,double epsi=EPSILON);
        double angle(const dvec3 &v)  const;
        dvec3 closestPointOnLine(const dvec3 &vA, const dvec3 &vB);
        dvec3 closestPointOnSegment(const dvec3 &vA, const dvec3 &vB);
        double projectionOnLine(const dvec3 &vA, const dvec3 &vB);
        dvec3 lerp(dvec3 &u, dvec3 &v, double factor);
        double distance(const dvec3& a_vector) const;
        dvec3 Rotation(const dvec3 &n,const double &angle) const;

        %extend {
            char *__str__() {
                static char temp[256];
                sprintf(temp,"[ %g, %g, %g ]", $self->x,$self->y,$self->z);
                return &temp[0];
            }
            char *__repr__() {
                static char temp[256];
                sprintf(temp,"vec3(%g,%g,%g)", $self->x,$self->y,$self->z);
                return &temp[0];
            }
            int __len__() {
                return 3;
            }
            double __getitem__(const int& ind)
            {
                return $self->v[ind];
            }
        }
    };
    class ivec3
    {
        public:
            ivec3(void);
            ivec3(long _a,long _b,long _c);
            ivec3(const ivec3 &iv);

            int operator==(const ivec3 &iv);
            int operator!=(const ivec3 &iv);

            const ivec3 operator*(long _i) const;
            const ivec3 operator/(long _i) const;
            const ivec3 operator+(const ivec3 &iv) const;
            const ivec3 operator-() const;
            const ivec3 operator-(const ivec3 &iv) const;

            ivec3 &operator*=(long _i);
            ivec3 &operator/=(long _i);
            ivec3 &operator+=(const ivec3 &iv);
            ivec3 &operator-=(const ivec3 &iv);

            long operator*(const ivec3 &iv) const;
            void set(long _a,long _b,long _c);
            void set(int tab[3]);
            void reset(void);
            %extend {
            char *__str__() {
                static char temp[256];
                sprintf(temp,"[ %li, %li, %li ]", $self->a,$self->b,$self->c);
                return &temp[0];
            }
            char *__repr__() {
                static char temp[256];
                sprintf(temp,"vec3(%li,%li,%li)", $self->a,$self->b,$self->c);
                return &temp[0];
            }
            int __len__() {
                return 3;
            }
            long __getitem__(const int& ind)
            {
                return $self->i[ind];
            }
        }
   };
};
namespace ScalarFieldBuilders
{
    using namespace core_mathlib;

    class ScalarFieldCreator
	{
        public:
            ScalarFieldCreator(const double& resolution);
            %rename(first_step_params) FirstStep_Params;
            void FirstStep_Params(const dvec3& boxMin,const dvec3& boxMax);
            %rename(third_step_volumescreator) ThirdStep_VolumesCreator;
            void ThirdStep_VolumesCreator();
            %rename(get_volume_value) GetVolumeValue;
            double GetVolumeValue(short volId);
            %rename(get_volume_count) GetVolumeCount;
            short GetVolumeCount();
            %rename(get_center_cell_coordinates) GetCenterCellCoordinates;
            dvec3 GetCenterCellCoordinates( ivec3 cell_id) const;
            %rename(get_cell_id_by_coord) GetCellIdByCoord;
            ivec3 GetCellIdByCoord(const dvec3& position);
            %rename(get_matrix_value) GetMatrixValue;
            short GetMatrixValue(const ivec3& index);
            %rename(get_domain_size) GetDomainSize;
            unsigned int GetDomainSize();
            %rename(copy_matrix) CopyMatrix;
            void CopyMatrix(short* INPLACE_ARRAY3,int DIM1,int DIM2,int DIM3,const ivec3& extractPos);
            %rename(copy_matrix_filtered) CopyMatrixFiltered;
            void CopyMatrixFiltered(short* INPLACE_ARRAY3,int DIM1,int DIM2,int DIM3,const ivec3& extractPos,short* IN_ARRAY1,int DIM1 );
            %rename(get_cell_value_boundaries) GetCellValueBoundaries;
            void GetCellValueBoundaries(ivec3& min,ivec3& max,const short& volid);
            %rename(get_first_volume_index) GetFirstVolumeIndex;
            short GetFirstVolumeIndex();
    };
    class TriangleScalarFieldCreator : public ScalarFieldCreator
    {
        public:
            TriangleScalarFieldCreator(const double& _resolution);
            %rename(second_step_pushtri) SecondStep_PushTri;
            void SecondStep_PushTri(const dvec3& A,const dvec3& B,const dvec3& C,const short& marker=1);
            %rename(load_ply_model) LoadPlyModel;
            bool LoadPlyModel(const std::string& fileInput);
    };
};
