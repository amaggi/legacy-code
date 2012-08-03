/*
 ******************************************************************************
 *
 *	Function hdclrg_
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
unsigned long int hd_fore;
unsigned long int hd_back;
XRectangle hd_rect;
int hd_batch;

int hd_pscolor;

int hd_fd;
int hd_lon;

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
hdclrg_(x,y,w,h)

int *x;
int *y;
int *w;
int *h;

/*
 *	hdclrg_ will cause a region to be cleared
 */

{
	hd_rect.x = *x;
	hd_rect.y = *y;
	hd_rect.width = *w;
	hd_rect.height = *h;
	XSetForeground (hd_display, hd_gc, hd_back);
	XFillRectangle (hd_display, hd_pixmap, hd_gc, *x, *y, *w, *h);
	if (!hd_batch)
		XFillRectangle (hd_display, hd_window, hd_gc, *x, *y, *w, *h);
	XSetForeground (hd_display, hd_gc, hd_fore);
	if (!hd_batch)
		XFlush (hd_display);
}

void
hdclrgl_(x,y,w,h)

int *x;
int *y;
int *w;
int *h;

/*
 *	hdclrgl_ will cause a PostScript region to be cleared 
 */

{
	static char line[80];

	if (hd_fd < 1) return;
	if (!hd_lon) return;
	sprintf (line, "gsave\n");
	if (hd_pscolor != 1) return;
	write (hd_fd, line, strlen(line));
	sprintf (line, "newpath\n");
	write (hd_fd, line, strlen(line));
        if (hd_back_light == 0.0 || hd_back_light == 1.0 || hd_back_sat==0.0) {
                sprintf (line, "%.3f setgray\n", hd_back_light);
        } else {
                sprintf (line, "%.3f %.3f %.3f setrgbcolor\n",
                                hd_back_red, hd_back_green, hd_back_blue);
        }
	write (hd_fd, line, strlen(line));
	sprintf(line, "%d %d moveto\n", (*x)+150, (*y)+150);
	write (hd_fd, line, strlen(line));
	sprintf(line, "%d 0 rlineto\n", (*w));
	write (hd_fd, line, strlen(line));
	sprintf(line, "0 %d rlineto\n", (*h));
	write (hd_fd, line, strlen(line));
	sprintf(line, "%d neg 0 rlineto\n", (*w));
	write (hd_fd, line, strlen(line));
	sprintf (line, "closepath\n");
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
}

