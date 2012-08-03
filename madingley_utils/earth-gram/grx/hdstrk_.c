/*
 ******************************************************************************
 *
 *	Function hdstrk_
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
int hd_npts;
XPoint hd_points[1001];
int hd_batch;

void
hdstrk_()

/*
 *	hdstrk_ will cause the plot to be stroked 
 */

{
	if (hd_npts > 1) {
		XDrawLines (hd_display, hd_pixmap, hd_gc, hd_points, hd_npts,
							CoordModeOrigin);
		if (!hd_batch)
			XDrawLines (hd_display, hd_window, hd_gc, hd_points, 
						hd_npts, CoordModeOrigin);
		hd_points[0].x = hd_points[hd_npts-1].x;
		hd_points[0].y = hd_points[hd_npts-1].y;
		hd_npts = 1;
	}
	if (!hd_batch) XFlush (hd_display);
}

long hd_fd;
int hd_lon;
int hd_nptsl;

void
hdstrkl_()

{
	if (hd_fd > 0 && hd_lon) {
		if (hd_nptsl > 0) {
			write (hd_fd, "stroke\n", 7);
			write (hd_fd, "newpath\n", 8);
			hd_nptsl = 0;
		}
	}
}
