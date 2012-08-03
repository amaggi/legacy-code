/*
 *******************************************************************************
 *
 *   Subroutine contr
 *
 *   Author - Danny Harvey
 *
 *******************************************************************************
 */

#include <stdio.h>

void
contrd_ ( nx, x, xtype, ny, y, ytype,
          nxmax, z, ztype, zstart, dz, nz, nptsm,
          xout, yout, 
          len_xtype, len_ytype, len_ztype)

int *     nx;
float *       x;
char *           xtype;
int *                   ny;
float *                     y;
char *                         ytype;
int *     nxmax;
float *          z;
char *              ztype;
float *                    zstart;
float *                            dz;
int *                                  nz;
int *                                      nptsm;
float *   xout;
float *         yout;

int       len_xtype, len_ytype, len_ztype;

/*
 *   Subroutine contrd will compute the X-Y coordinates of constant
 *   Z value contours for some Z = Z(X,Y).  Z(X,Y) must be sampled
 *   along some rectangular grid of X-Y values although the grid
 *   need not have equal spacing everywhere.
 *
 *   inputs  - nx     = number of different X values in the grid.
 *             x(nx)  = X sampling values.  These must be ordered
 *                      but can be either increasing or decreasing.
 *             xtype  = CHARACTER*3 flag which indicates whether to
 *                      return linear ('LIN') or logarithmic ('LOG')
 *                      contour values of X in xout
 *             ny     = number of different Y values in the grid.
 *             y(ny)  = Y sampling values.  These must be ordered
 *                      but can be either increasing or decreasing.
 *             ytype  = CHARACTER*3 flag which indicates whether to
 *                      return linear ('LIN') or logarithmic ('LOG')
 *                      contour values of Y in yout
 *             z(nx*ny) = sampled Z values for each X-Y sampling point.
 *                        This array is indexed by X values first.
 *             ztype  = CHARACTER*3 flag which indicates whether to
 *                      make the contour spacing linear ('LIN') or
 *                      logarithmic ('LOG')
 *             zstart = first contour Z value.
 *             dz     = contour Z value increment.
 *             nz     = number of contours to be computed.
 *             nptsm  = maximum number of values in the xout and yout
 *                      arrays.
 *             xout(nptsm)
 *	              = Array used for storing X-contour plot points.
 *             yout(nptsm)
 *	              = Array used for storing Y-contour plot points.
 *
 */

{
	static short int *izxy;
	static float *zs, *zxmin, *zxmax, *zymin, *zymax;
	int n;

	if ((*nx) < 2 || (*ny) < 2) return;
	n = ((*nx) >= (*ny)) ? (*nx) : (*ny);
	n = (n >= (*nz)) ? n : (*nz);
	zs = (float *) malloc (n*sizeof(float));
	if (zs == NULL) {
		fprintf (stderr, "contrd: Malloc error.\n");
		return;
	}
	zxmin = (float *) malloc (n*sizeof(float));
	if (zxmin == NULL) {
		fprintf (stderr, "contrd: Malloc error.\n");
		free (zs);
		return;
	}
	zxmax = (float *) malloc (n*sizeof(float));
	if (zxmax == NULL) {
		fprintf (stderr, "contrd: Malloc error.\n");
		free (zs);
		free (zxmin);
		return;
	}
	zymin = (float *) malloc (n*sizeof(float));
	if (zymin == NULL) {
		fprintf (stderr, "contrd: Malloc error.\n");
		free (zs);
		free (zxmin);
		free (zxmax);
		return;
	}
	zymax = (float *) malloc (n*sizeof(float));
	if (zymax == NULL) {
		fprintf (stderr, "contrd: Malloc error.\n");
		free (zs);
		free (zxmin);
		free (zxmax);
		free (zymin);
		return;
	}
	n = (*nx)*(*ny)*4;
	izxy = (short int *) malloc (n*sizeof(short int));
	if (izxy == NULL) {
		fprintf (stderr, "contrd: Malloc error.\n");
		free (zs);
		free (zxmin);
		free (zxmax);
		free (zymin);
		free (zymax);
		return;
	}
	contr_ ( nx, x, xtype, ny, y, ytype,
          	nxmax, z, ztype, zstart, dz, nz, nptsm,
          	xout, yout, 
     		izxy, zs, zxmin, zxmax, zymin, zymax,
          	len_xtype, len_ytype, len_ztype);
	free (zs);
	free (zxmin);
	free (zxmax);
	free (zymin);
	free (zymax);
	free (izxy);
}
