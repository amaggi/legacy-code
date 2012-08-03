/*
 ******************************************************************************
 *
 *	Function hdsetclp
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

void
hdsetclp_()

/*
 *	hdsetclp_() will set X-window clipping according to the current
 *	viewport defined by setdim-setscl.
 */

{
	XRectangle rect;
	int clipx, clipy;
	float xmin, xmax, ymin, ymax;
	float xll, yll, xur, yur;

	getclp_ (&xmin, &xmax, &ymin, &ymax);
	getfrm_ (&xll, &xur, &yll, &yur);
	clipx = xmin + xll + 0.5;
	clipy = yur - (ymax + yll) + 0.5;
	rect.x = 0; rect.y = 0;
	rect.width = xmax - xmin + 1.5; rect.height = ymax - ymin + 1.5;
	XSetClipRectangles (hd_display, hd_gc, clipx, clipy,
							&rect, 1, Unsorted);
}

/*
 ******************************************************************************
 *
 *	Function hdclrclp
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

void
hdclrclp_()

/*
 *	hdclrclp_() will clear X-window clipping.
 */

{
	XSetClipMask (hd_display, hd_gc, None);
}

/*
 ******************************************************************************
 *
 *	Function hdsetclpl
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

int hd_fd;
int hd_lon;

void
hdsetclpl_()

/*
 *	hdsetclpl_() will set PostScript clipping according to the current
 *	viewport defined by setdim-setscl.
 */

{
	float xcl,ycb,fl;
	float x, y, w, h;
	static char line[128];

	if (hd_fd < 1) return;
	if (!hd_lon) return;
	getfl_ (&xcl,&ycb,&fl);
	getclp_ (&x, &w, &y, &h);
	x += xcl;
	y += ycb;
	x = x*fl + 150.0;
	y = y*fl + 150.0;
	w += xcl;
	h += ycb;
	w = w*fl + 150.0;
	h = h*fl + 150.0;
	w = w - x + 1.0;
	h = h - y + 1.0;
	sprintf (line, "gsave %.1f %.1f %.1f %.1f rectcl\n", x, y, w, h);
	write (hd_fd, line, strlen(line));
}

/*
 ******************************************************************************
 *
 *	Function hdclrclpl
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

void
hdclrclpl_()

/*
 *	hdclrclpl_() will clear PostScript clipping.
 */

{
	if (hd_fd < 1) return;
	if (!hd_lon) return;
	write (hd_fd, "grestore\n", strlen("grestore\n"));
}
