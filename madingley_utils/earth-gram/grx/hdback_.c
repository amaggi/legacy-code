/*
 ******************************************************************************
 *
 *	Function hdback_
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
Colormap hd_cmap;
Pixmap hd_pixmap;
GC hd_gc;
Visual *hd_visual;
unsigned long int hd_back;
#ifdef STELLAR
XFDIGC hd_gc3;
#endif

float hd_back_light;
float hd_back_sat;
float hd_back_red;
float hd_back_green;
float hd_back_blue;

void
hdback_(hue, light, sat)

float *hue;
float *light;
float *sat;

/*
 *	hdback_ will set the background color.
 */

{
	unsigned long int hdhls_to_pixel();

	hd_back = hdhls_to_pixel (*hue, *light, *sat);
	hdhls_to_rgb (*hue, *light, *sat, 
				&hd_back_red, &hd_back_green, &hd_back_blue);
	hd_back_light = *light;
	hd_back_sat = *sat;
	XSetBackground (hd_display, hd_gc, hd_back);
}
