#include <math.h>
#include <stdio.h>
#include <string.h>

/* X-tests */
#define X 0
/* Y-tests */
#define Y 1
/* Z-tests */
#define Z 2

#define CROSS(dest, v1, v2) \
    dest[0] = v1[1] * v2[2] - v1[2] * v2[1]; \
    dest[1] = v1[2] * v2[0] - v1[0] * v2[2]; \
    dest[2] = v1[0] * v2[1] - v1[1] * v2[0];

#define DOT(v1, v2) \
    (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])

#define SUB(dest, v1, v2) \
    dest[0] = v1[0] - v2[0]; \
    dest[1] = v1[1] - v2[1]; \
    dest[2] = v1[2] - v2[2];

#define FINDMINMAX(x0, x1, x2, min, max) \
    min = max = x0; \
    if (x1 < min) min = x1; \
    if (x1 > max) max = x1; \
    if (x2 < min) min = x2; \
    if (x2 > max) max = x2;

namespace boxtri_test
{
	int planeBoxOverlap(double normal[3], double d, double maxbox[3])
	{
		int q;
		double vmin[3], vmax[3];

		for (q = X; q <= Z; q++)
		{
			if (normal[q] > 0.0)
			{
				vmin[q] = -maxbox[q];
				vmax[q] = maxbox[q];
			}
			else
			{
				vmin[q] = maxbox[q];
				vmax[q] = -maxbox[q];
			}
		}

		if (DOT(normal, vmin) + d > 0.0)
			return 0;

		if (DOT(normal, vmax) + d >= 0.0)
			return 1;

		return 0;
	}

	/*======================== X-tests ========================*/

	#define AXISTEST_X01(a, b, fa, fb) \
		p0 = a * v0[Y] - b * v0[Z]; \
		p2 = a * v2[Y] - b * v2[Z]; \
		if (p0 < p2) { min = p0; max = p2; } \
		else { min = p2; max = p0; } \
		rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z]; \
		if (min > rad || max < -rad) return 0;

	#define AXISTEST_X2(a, b, fa, fb) \
		p0 = a * v0[Y] - b * v0[Z]; \
		p1 = a * v1[Y] - b * v1[Z]; \
		if (p0 < p1) { min = p0; max = p1; } \
		else { min = p1; max = p0; } \
		rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z]; \
		if (min > rad || max < -rad) return 0;

	/*======================== Y-tests ========================*/

	#define AXISTEST_Y02(a, b, fa, fb) \
		p0 = -a * v0[X] + b * v0[Z]; \
		p2 = -a * v2[X] + b * v2[Z]; \
		if (p0 < p2) { min = p0; max = p2; } \
		else { min = p2; max = p0; } \
		rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z]; \
		if (min > rad || max < -rad) return 0;

	#define AXISTEST_Y1(a, b, fa, fb) \
		p0 = -a * v0[X] + b * v0[Z]; \
		p1 = -a * v1[X] + b * v1[Z]; \
		if (p0 < p1) { min = p0; max = p1; } \
		else { min = p1; max = p0; } \
		rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z]; \
		if (min > rad || max < -rad) return 0;

	/*======================== Z-tests ========================*/

	#define AXISTEST_Z12(a, b, fa, fb) \
		p1 = a * v1[X] - b * v1[Y]; \
		p2 = a * v2[X] - b * v2[Y]; \
		if (p2 < p1) { min = p2; max = p1; } \
		else { min = p1; max = p2; } \
		rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y]; \
		if (min > rad || max < -rad) return 0;

	#define AXISTEST_Z0(a, b, fa, fb) \
		p0 = a * v0[X] - b * v0[Y]; \
		p1 = a * v1[X] - b * v1[Y]; \
		if (p0 < p1) { min = p0; max = p1; } \
		else { min = p1; max = p0; } \
		rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y]; \
		if (min > rad || max < -rad) return 0;

	/*======================== Main function ========================*/

	/*
	* Function: int triBoxOverlap(double boxcenter[3],
	*                             double boxhalfsize[3],
	*                             double triverts[3][3])
	*
	* Description:
	*   Check whether a triangle overlaps an axis-aligned bounding box.
	*
	* Return:
	*   0 = no overlap
	*   1 = overlap
	*/

	int triBoxOverlap(
		double boxcenter[3],
		double boxhalfsize[3],
		double triverts[3][3])
	{
		/* Local variables */
		double v0[3], v1[3], v2[3];
		double axis[3];
		double min, max, p0, p1, p2, rad;
		double fex, fey, fez;
		double normal[3], e0[3], e1[3], e2[3];
		double d;

		/* Move the triangle so that the boxcenter is in the origin */
		SUB(v0, triverts[0], boxcenter);
		SUB(v1, triverts[1], boxcenter);
		SUB(v2, triverts[2], boxcenter);

		/* Compute triangle edges */
		SUB(e0, v1, v0);
		SUB(e1, v2, v1);
		SUB(e2, v0, v2);

		/* Compute absolute values of edge components */
		fex = fabs(e0[X]);
		fey = fabs(e0[Y]);
		fez = fabs(e0[Z]);

		/* Test the 9 axes given by the cross products of the edges and the box axes */

		/* First edge */
		AXISTEST_X01(e0[Z], e0[Y], fez, fey);
		AXISTEST_Y02(e0[Z], e0[X], fez, fex);
		AXISTEST_Z12(e0[Y], e0[X], fey, fex);

		/* Second edge */
		fex = fabs(e1[X]);
		fey = fabs(e1[Y]);
		fez = fabs(e1[Z]);

		AXISTEST_X01(e1[Z], e1[Y], fez, fey);
		AXISTEST_Y02(e1[Z], e1[X], fez, fex);
		AXISTEST_Z0(e1[Y], e1[X], fey, fex);

		/* Third edge */
		fex = fabs(e2[X]);
		fey = fabs(e2[Y]);
		fez = fabs(e2[Z]);

		AXISTEST_X2(e2[Z], e2[Y], fez, fey);
		AXISTEST_Y1(e2[Z], e2[X], fez, fex);
		AXISTEST_Z12(e2[Y], e2[X], fey, fex);

		/* Test overlap in the X direction */
		FINDMINMAX(v0[X], v1[X], v2[X], min, max);
		if (min > boxhalfsize[X] || max < -boxhalfsize[X])
			return 0;

		/* Test overlap in the Y direction */
		FINDMINMAX(v0[Y], v1[Y], v2[Y], min, max);
		if (min > boxhalfsize[Y] || max < -boxhalfsize[Y])
			return 0;

		/* Test overlap in the Z direction */
		FINDMINMAX(v0[Z], v1[Z], v2[Z], min, max);
		if (min > boxhalfsize[Z] || max < -boxhalfsize[Z])
			return 0;

		/* Test if the box intersects the plane of the triangle */

		CROSS(normal, e0, e1);
		d = -DOT(normal, v0);

		if (!planeBoxOverlap(normal, d, boxhalfsize))
			return 0;

		return 1;
	}
} // namespace boxtri_test
