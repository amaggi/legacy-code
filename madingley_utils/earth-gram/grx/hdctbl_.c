/*
 ******************************************************************************
 *
 *	Function hdctbl_
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

#include <stdio.h>
#include <X11/Xlib.h>

Display *hd_display;
Colormap hd_cmap;
Pixmap hd_pixmap;
GC hd_gc;
Visual *hd_visual;
unsigned long int hd_back;

void
hdback_(nc, hue, light, sat, ctgflg, pixel, ierr)

int *nc;
float *hue;
float *light;
float *sat;
int *ctgflg;
int *pixel;
int *ierr;

/*
 *	hdctbl_ will define a color table.
 */

{
	unsigned long int hdhls_to_pixel();

	hd_back = hdhls_to_pixel (*hue, *light, *sat);
	XSetBackground (hd_display, hd_gc, hd_back);
}
