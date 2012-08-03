/*
 ******************************************************************************
 *
 *	Function hdfplt_
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

#include <stdio.h>
#include <X11/Xlib.h>
#ifdef STELLAR
#include <xfdi.h>
#endif

Display *hd_display;
Pixmap hd_pixmap;
GC hd_gc;
XRectangle hd_rect;
#ifdef STELLAR
XFDIGC hd_gc3;
float hd_matrix[16];
#endif

void
hdfplt_(n,x,y,xr,yr,wr,hr,xmin,xmax,ymin,ymax)

int *n;
float *x;
float *y;
float *xr,*yr,*wr,*hr;
float *xmin,*xmax,*ymin,*ymax;

/*
 *	hdfplt_ will cause a fast polyline plot.
 */

{
#ifdef STELLAR_OFF
	int nn;

	nn = *n;
	hdstrk_();
	hd_matrix[0] = *wr / ( *xmax - *xmin );
	hd_matrix[5] = *hr / ( *ymax - *ymin );
	hd_matrix[12] = *xr - ( *xmin * hd_matrix[0] );
	hd_matrix[13] = *yr - ( *ymin * hd_matrix[5] );
	XFDIClear (hd_display, hd_gc3, 0, 0, &hd_rect);
	XFDISetMatrix (hd_display, hd_gc3, XFDIObjectMatrix, hd_matrix, 
								XFDIMatrix2d);
	xfl2_ (&hd_display, &hd_gc3, x, y, &nn);
#else
#endif
}
