/*
 ******************************************************************************
 *
 *	Function hddraw_
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

#include <X11/Xlib.h>

Pixmap hd_pixmap;
int hd_npts;
XPoint hd_points[1001];

void
hddraw_(fx,fy,fl)

float *fx,*fy,*fl;

/*
 *	hddraw_ will draw to (ix,iy)
 */

{
	short x, y;

	if (hd_npts < 1) hdmove_ (fx,fy,fl);
	x = *fx + 0.5;
	y = *fy + 0.5;
	if (x != hd_points[hd_npts-1].x || y != hd_points[hd_npts-1].y) {
		hd_points[hd_npts].x = *fx + 0.5;
		hd_points[hd_npts].y = *fy + 0.5;
		hd_npts++;
	}
	if (hd_npts > 1000) {
		hdstrk_();
		hdmove_(fx,fy,fl);
	}
}

int hd_fd;
int hd_lon;
int hd_oldx, hd_oldy;
int hd_nptsl;
int hd_ldraw = 0;
int hd_von;

void
hddrawl_(fx,fy,fl)

float *fx,*fy,*fl;

{
	int ix, iy;
	int x, y, rx, ry;
	char line[20];

	if (hd_fd < 1) return;
	if (!hd_lon) return;
	if (!hd_von) return;
	ix = (*fx)*(*fl) + 0.5;
	iy = (*fy)*(*fl) + 0.5;
	x = ix+150;
	y = iy+150;
	rx = x - hd_oldx;
	ry = y - hd_oldy;
	switch (ry) {
	case 2:
		switch (rx) {
		case -2:
			strcpy (line, "A\n");
			break;
		case -1:
			strcpy (line, "B\n");
			break;
		case 0:
			strcpy (line, "C\n");
			break;
		case 1:
			strcpy (line, "D\n");
			break;
		case 2:
			strcpy (line, "E\n");
			break;
		default:
			sprintf(line, "%d %d r\n", rx, ry);
			break;
		}
		break;
	case 1:
		switch (rx) {
		case -2:
			strcpy (line, "F\n");
			break;
		case -1:
			strcpy (line, "G\n");
			break;
		case 0:
			strcpy (line, "H\n");
			break;
		case 1:
			strcpy (line, "I\n");
			break;
		case 2:
			strcpy (line, "J\n");
			break;
		default:
			sprintf(line, "%d %d r\n", rx, ry);
			break;
		}
		break;
	case 0:
		switch (rx) {
		case -2:
			strcpy (line, "K\n");
			break;
		case -1:
			strcpy (line, "L\n");
			break;
		case 0:
			return;
		case 1:
			strcpy (line, "M\n");
			break;
		case 2:
			strcpy (line, "N\n");
			break;
		default:
			sprintf(line, "%d x\n", rx);
			break;
		}
		break;
	case -1:
		switch (rx) {
		case -2:
			strcpy (line, "O\n");
			break;
		case -1:
			strcpy (line, "P\n");
			break;
		case 0:
			strcpy (line, "Q\n");
			break;
		case 1:
			strcpy (line, "R\n");
			break;
		case 2:
			strcpy (line, "S\n");
			break;
		default:
			sprintf(line, "%d %d r\n", rx, ry);
			break;
		}
		break;
	case -2:
		switch (rx) {
		case -2:
			strcpy (line, "T\n");
			break;
		case -1:
			strcpy (line, "U\n");
			break;
		case 0:
			strcpy (line, "V\n");
			break;
		case 1:
			strcpy (line, "W\n");
			break;
		case 2:
			strcpy (line, "X\n");
			break;
		default:
			sprintf(line, "%d %d r\n", rx, ry);
			break;
		}
		break;
	default:
		switch (rx) {
		case 0:
			sprintf(line, "%d y\n", ry);
			break;
		default:
			sprintf(line, "%d %d r\n", rx, ry);
			break;
		}
		break;
	}
	write (hd_fd, line, strlen(line));
	hd_ldraw = 1;
	hd_nptsl++;
	hd_oldx = x;
	hd_oldy = y;
	if (hd_nptsl > 1000) {
		hdstrkl_();
		hdmovel_ (fx,fy,fl);
	}
}
