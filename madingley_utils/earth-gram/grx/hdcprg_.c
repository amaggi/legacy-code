/*
 ******************************************************************************
 *
 *	Function hdcprg_
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
int hd_width;
int hd_height;

void
hdcprg_(x,y,w,h)

int *x;
int *y;
int *w;
int *h;

/*
 *	hdcprg_ will cause a region to be cleared
 */

{
	XCopyArea (hd_display, hd_pixmap, hd_window, hd_gc, *x, *y,
							*w, *h, *x, *y);
	XFlush (hd_display);
}

void
hdrefresh_ ()

/*
 *	hdrefresh_ will refresh the whole screen.
 */

{
	XCopyArea (hd_display, hd_pixmap, hd_window, hd_gc, 0, 0,
						hd_width, hd_height, 0, 0);
	XFlush (hd_display);
}
