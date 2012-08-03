/*
 ******************************************************************************
 *
 *	Function hdclr_
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

#include <X11/Xlib.h>

Pixmap hd_pixmap;
int hd_width;
int hd_height;
int hd_fd;
int hd_lon;
int hd_lstart;

void
hdclr_()

/*
 *	hdclr_ will cause a new page to occur
 */

{
	int x, y, w, h;

	void hdclrl();

	hdstrk_();
	x = 0;
	y = 0;
	w = hd_width;
	h = hd_height;
	hdclrg_ (&x, &y, &w, &h);
}

int hd_ldraw;

void
hdclrl_(x, y, w, h)

int *x, *y, *w, *h;

{
	if (hd_fd < 0) return;
	if (!hd_lon) return;
	hdstrkl_();
	if (hd_ldraw && !hd_lstart) write (hd_fd, "showpage\n", 9);
	hd_ldraw = 0;
        hd_lstart = 0;
	hdclrgl_ (x, y, w, h);
}
