/*
 ******************************************************************************
 *
 *	Function hdmptr_
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

#include <X11/Xlib.h>

Display *hd_display;
Window hd_window;

void
hdmptr_(ixmin, ixmax, iymin, iymax, ix, iy)

long *ixmin,*iymin;
long *ixmax,*iymax;
long *ix,*iy;

/*
 *	hdmptr_ will cause the pointer to be moved.
 */

{
	int x, y;

	hdstrk_();
	x = *ix;
	if (x < *ixmin) x = *ixmin;
	if (x > *ixmax) x = *ixmax;
	y = *iy;
	if (y < *iymin) y = *iymin;
	if (y > *iymax) y = *iymax;
	XWarpPointer (hd_display, None, hd_window, 0, 0, 0, 0, x, y);
}
