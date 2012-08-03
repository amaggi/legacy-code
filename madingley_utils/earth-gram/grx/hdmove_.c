/*
 ******************************************************************************
 *
 *	Function hdmove_
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

#include <X11/Xlib.h>

Pixmap hd_pixmap;
XPoint hd_points[1001];
int hd_npts;

void
hdmove_(fx,fy,fl)

float *fx,*fy,*fl;

/*
 *	hdmove_ will move to (ix,iy)
 */

{
	hdstrk_();
	hd_points[0].x = *fx + 0.5;
	hd_points[0].y = *fy + 0.5;
	hd_npts = 1;
}

int hd_fd;
int hd_lon;
int hd_oldx, hd_oldy;
int hd_von = 1;

void
hdmovel_(fx,fy,fl)

float *fx,*fy,*fl;


{
	int ix, iy;
	char line[20];

	if (hd_fd > 0 && hd_lon && hd_von) {
		hdstrkl_();
		ix = (*fx)*(*fl) + 0.5;
		iy = (*fy)*(*fl) + 0.5;
		sprintf(line, "%d %d m\n", ix+150, iy+150);
		write (hd_fd, line, strlen(line));
		hd_oldx = ix+150;
		hd_oldy = iy+150;
	}
}
