/*
 ******************************************************************************
 *
 *	Function hdtext_
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

int hd_fd;
int hd_lon;
int hd_von;
int hd_lfont = 0;

void
hdtext_(x, y, angle, iref, string, iclip, lstr)

float *x;
float *y;
float *angle;
int *iref;
char *string;
int *iclip;
int lstr;

{
	static char line[256];
	float xl, yl, size, fl, xcl, xcr, ycb, yct;
	int icl;

	char *convert_string();
	float xmap_();
	float ymap_();

	if (hd_fd < 1) return;
	if (!hd_lon) return;
	sprintf (line, "%%%% text: x=%f y=%f angle=%f iref=%d\n",
						*x, *y, *angle, *iref);
	write (hd_fd, line, strlen(line));
	sprintf (line, "%%%%       string='%s' iclip=%d\n", 
					convert_string(string, lstr), *iclip);
	write (hd_fd, line, strlen(line));
	switch (hd_lfont) {
	default:
	case 0:
		hd_von = 1;
		break;
	case 112:
	case 115:
	case 116:
	case 130:
	case 131:
	case 132:
		hd_von = 0;
		gethit_ (&size);
		getfl_ (&xcl,&ycb,&fl);
		getxmp_ (x, &xl);
		getymp_ (y, &yl);
		xl += xcl;
		yl += ycb;
		xl = xl*fl + 150.0;
		yl = yl*fl + 150.0;
		size *= fl;
		if (*iclip) icl = 0; else icl = 1;
		getclp_ (&xcl, &xcr, &ycb, &yct);
		xcl = xcl*fl + 150.0;
		xcr = xcr*fl + 150.0;
		ycb = ycb*fl + 150.0;
		yct = yct*fl + 150.0;
		sprintf (line, "npfont%d\n", hd_lfont);
		write (hd_fd, line, strlen(line));
		sprintf (line, 
		"%.0f %.0f %.0f npchsc %.0f %d %.0f %.0f %.0f %.0f %d (%s) nptext\n", 
				xl, yl, size, *angle, *iref, xcl, xcr, ycb, yct,
				icl, convert_string(string, lstr));
		write (hd_fd, line, strlen(line));
		break;
	}
}

void
hdtxtf_ ()

{
	static char line[256];

	if (hd_fd < 1) return;
	if (!hd_lon) return;
	sprintf (line, "%%%% text: finished\n");
	write (hd_fd, line, strlen(line));
	hd_von = 1;
}

void
hdstlf_ (lfont)

int *lfont;

{
	hd_lfont = *lfont;
}
