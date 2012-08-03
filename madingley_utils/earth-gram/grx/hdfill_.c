/*
 ******************************************************************************
 *
 *	Function hdfill_
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
unsigned long int hd_fore;
unsigned long int hd_back;

void
hdfill_(iclr)

int *   iclr;

/*
 *	hdfill_ will cause a region to be filled
 */

{
	if (hd_npts > 1) {
		XSetForeground (hd_display, hd_gc, hd_back);
		XFillPolygon (hd_display, hd_pixmap, hd_gc, hd_points, hd_npts,
							Complex, CoordModeOrigin);
		if (!hd_batch)
			XFillPolygon (hd_display, hd_window, hd_gc, hd_points, 
						hd_npts, Complex, CoordModeOrigin);
		if (hd_npts >= 1000 || *iclr) {
			hd_points[0].x = hd_points[hd_npts-1].x;
			hd_points[0].y = hd_points[hd_npts-1].y;
			hd_npts = 1;
		}
		XSetForeground (hd_display, hd_gc, hd_fore);
	}
	if (!hd_batch) XFlush (hd_display);
}

int hd_pscolor;

int hd_fd;
int hd_lon;
int hd_nptsl;

float hd_fore_light;
float hd_fore_sat;
float hd_fore_red;
float hd_fore_green;
float hd_fore_blue;
float hd_back_light;
float hd_back_sat;
float hd_back_red;
float hd_back_green;
float hd_back_blue;

void
hdfilll_(iclr)

int *    iclr;

{
	static char line[80];

	if (hd_fd < 1) return;
	if (!hd_lon) return;
	if (hd_nptsl < 1) return;
/*	if (hd_pscolor != 1) return;*/
	sprintf (line, "gsave\n");
	write (hd_fd, line, strlen(line));
        if (hd_back_light == 0.0 || hd_back_light == 1.0 || hd_back_sat==0.0) {
                sprintf (line, "%.3f setgray\n", hd_back_light);
        } else {
                sprintf (line, "%.3f %.3f %.3f setrgbcolor\n",
                                hd_back_red, hd_back_green, hd_back_blue);
        }
	write (hd_fd, line, strlen(line));
	sprintf (line, "fill\n");
	write (hd_fd, line, strlen(line));
        if (hd_fore_light == 0.0 || hd_fore_light == 1.0 || hd_fore_sat==0.0) {
                sprintf (line, "%.3f setgray\n", hd_fore_light);
        } else {
                sprintf (line, "%.3f %.3f %.3f setrgbcolor\n",
                                hd_fore_red, hd_fore_green, hd_fore_blue);
        }
	write (hd_fd, line, strlen(line));
	sprintf (line, "grestore\n");
	write (hd_fd, line, strlen(line));
	if (*iclr) {
		sprintf (line, "newpath\n");
		write (hd_fd, line, strlen(line));
	}
}
