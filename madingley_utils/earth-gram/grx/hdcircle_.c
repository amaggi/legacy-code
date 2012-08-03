/*
 ******************************************************************************
 *
 *	Function hdcircle
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

#include <X11/Xlib.h>

Display *hd_display;
Pixmap hd_pixmap;
Window hd_window;
GC hd_gc;
int hd_batch;

void
hdcircle_(x, y, r, iclose, iclip)

int *x;
int *y;
int *r;
int *iclose;
int *iclip;

/*
 *	hdcircle will cause a circle to be plotted
 */

{
	int w;

	hdstrk_();
	w = 2 * (*r);
	if (!(*iclip)) hdsetclp_ ();
	XDrawArc (hd_display, hd_pixmap, hd_gc, *x-*r, *y-*r, w, w, 0, 360*64);
	if (!hd_batch) {
		XDrawArc (hd_display, hd_window, hd_gc, *x-*r, *y-*r, 
							w, w, 0, 360*64);
		XFlush (hd_display);
	}
	if (!(*iclip)) hdclrclp_ ();
}

/*
 ******************************************************************************
 *
 *	Function hdcirclel_
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

int hd_fd;
int hd_lon;

void
hdcirclel_(x, y, r, iclose, iclip)

float *x;
float *y;
float *r;
int *iclose;
int *iclip;

/*
 *	hdcircle will cause a circle to be plotted
 */

{
	float xcl,ycb,fl;
	float xl, yl, rl, xscale;
	static char line[128];

	if (hd_fd < 1) return;
	if (!hd_lon) return;
	if (!(*iclip)) hdsetclpl_ ();
	getfl_ (&xcl,&ycb,&fl);
	getxscale_ (&xscale);
	getxmp_ (x, &xl);
	getymp_ (y, &yl);
	xl += xcl;
	yl += ycb;
	xl = xl*fl + 150.0;
	yl = yl*fl + 150.0;
	rl = (*r)*xscale*fl;
	sprintf (line, "newpath %.1f %.1f %.1f 0 360 arc closepath stroke\n", 
								xl, yl, rl);
	write (hd_fd, line, strlen(line));
	if (!(*iclip)) hdclrclpl_ ();
}
