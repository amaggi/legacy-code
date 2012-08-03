/*
 ******************************************************************************
 *
 *	Function hdopen_
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

#include <stdio.h>
#include <fcntl.h>
#include <X11/Xlib.h>
#include <X11/cursorfont.h>

#include "nxplotdlib_types.h"

#ifdef STELLAR
#include <xfdi.h>
#endif

#define	FONTNAME	"9x15"

Display *hd_display = NULL;
int hd_screen=0;
Pixmap hd_pixmap = 0;
Colormap hd_cmap = 0;
Window hd_iwindow = 0;
Window hd_window = 0;
Cursor hd_cursor = 0;
XFontStruct *hd_cfont = NULL;
int hd_width = 0;
int hd_height = 0;
int hd_depth = 0;
GC hd_gc = NULL;
Visual *hd_visual = NULL;
unsigned long hd_fore;
unsigned long hd_back;
int hd_batch = 0;
#ifdef STELLAR
XFDIGC hd_gc3;
float hd_matrix[16] = {1.0, 0.0, 0.0, 0.0,
		       0.0, 1.0, 0.0, 0.0,
		       0.0, 0.0, 1.0, 0.0,
		      -1.0, 0.0, 0.0, 1.0};
#endif

char hd_filename[256];
int hd_fd = 0;
int hd_lon = 0;
int hd_lstart = 0;
nxplotdLib hd_xpllib;

void 
hdopen_ (plot_fname, display_name, prog_name, xwin, ywin, wsize, hsize, fw, fh, 
							pid, lplt, ldis, lprg)

char *plot_fname;
char *display_name;
char *prog_name;
float *xwin;
float *ywin;
float *wsize;
float *hsize;
int lplt;
int ldis;
int lprg;
float *fw;
float *fh;
int *pid;

/*
 *	hdopen_ will connect to an X server, open an input only window
 *	for communication with the server, launch the plotting daemon,
 *	and initialize global variables relating to Xlib stuff.
 *
 *	Inputs  -	display_name	= Character string that defines
 *					  the X server to hook to (e.g.
 *					  "eldojr:0"). If this string is
 *					  blank, then the display is obtained
 *					  from the DISPLAY environment variable.
 *			prog_name	= Character string that defines the
 *					  program name. This is used in the
 *					  title bar and icon for the plotting
 *					  daemon window.
 *			wsize		= Floating number that defines the
 *					  width of the plotting window in total
 *					  screen heights. If wsize = 0, then
 *					  the full screen width is used.
 *			hsize		= Floating number that defines the
 *					  height of the plotting window in total
 *					  screen heights. If hsize = 0, then
 *					  the full screen height is used.
 */

{
	char *dname, *pname;
	static char args[256];
	Window root;
	int x, y, w, h, b, d;
	XColor fcol, bcol;

	char *convert_string();

/*
 *	Connect to the X server.
 */
	if (hd_display == NULL) {
		dname = convert_string (display_name, ldis);
		if (*dname == '\0') dname = NULL;
		hd_display = XOpenDisplay (dname);
		if (hd_display == NULL) {
			fprintf (stderr, 
				"hdopen: Cannot connect to X server %s.\n",
					XDisplayName (dname));
			exit (1);
		}
	} else {
		goto PS;
	}
	hd_screen = DefaultScreen (hd_display);
	hd_visual = DefaultVisual (hd_display, hd_screen);
	hd_cmap = DefaultColormap (hd_display, hd_screen);
	hd_depth = DefaultDepth (hd_display, hd_screen);
	w = DisplayWidth (hd_display, hd_screen);
	h = DisplayHeight (hd_display, hd_screen);
/*
 *	Create an input only window for this process.
 */
	hd_iwindow = XCreateWindow (hd_display, 
			RootWindow(hd_display, hd_screen),
			0, 0, 100, 100, 0, 0, InputOnly, CopyFromParent,
			NULL, NULL);
/*
 *	Launch the plotting daemon.
 */
	pname = convert_string (prog_name, lprg);
	if (*wsize > 0.0) {
		x = w * (*xwin);
		y = h * (*ywin);
		w = h * (*wsize);
	}
	if (*hsize > 0.0) {
		x = w * (*xwin);
		y = h * (*ywin);
		h *= (*hsize);
	}
	if (*wsize > 0.0 && *hsize > 0.0) {
		sprintf (args, "-geom %dx%d+%d+%d", w, h, x, y);
		nxplotd_launchdaemon (pname, hd_display, hd_iwindow, args,
								&hd_xpllib);
	} else {
		nxplotd_launchdaemon (pname, hd_display, hd_iwindow, NULL,
								&hd_xpllib);
	}
/*
 *	Fill in the other Xlib stuff.
 */
	*pid = hd_xpllib.daemon_pid;
	hd_pixmap = hd_xpllib.daemon_pixmap;
	hd_window = hd_xpllib.daemon_window;
	hd_gc = XCreateGC (hd_display, hd_pixmap, NULL, NULL);
	hd_cursor = XCreateFontCursor (hd_display, XC_crosshair);
	hd_cfont = XLoadQueryFont (hd_display, FONTNAME);
	if (hd_cfont == NULL) {
		fprintf (stderr, "hdopen: Cannot load cursor font %s.\n",
								FONTNAME);
		exit (1);
	}
	if (hd_depth > 2) {
		fcol.red = 0xffff;
		fcol.green = 0x0000;
		fcol.blue = 0x0000;
		bcol.red = 0xffff;
		bcol.green = 0x0000;
		bcol.blue = 0x0000;
		XRecolorCursor (hd_display, hd_cursor, &fcol, &bcol);
	} else {
		fcol.red = 0x0000;
		fcol.green = 0x0000;
		fcol.blue = 0x0000;
		bcol.red = 0x0000;
		bcol.green = 0x0000;
		bcol.blue = 0x0000;
		XRecolorCursor (hd_display, hd_cursor, &fcol, &bcol);
	}
	hd_fore = BlackPixel (hd_display, hd_screen);
	hd_back = WhitePixel (hd_display, hd_screen);
	XSetForeground (hd_display, hd_gc, hd_fore);
	XSetBackground (hd_display, hd_gc, hd_back);
	XGetGeometry (hd_display, hd_pixmap, &root, &x, &y, &w, &h, &b, &d);
	hd_width = w;
	hd_height = h;
	*fw = w;
	*fh = h;
#ifdef STELLAR_OFF
	hd_gc3 = XFDICreateGc3 (hd_display, hd_pixmap);
#endif
/*
 *	Open the postscript output file. 
 */
PS: 	pname = convert_string (plot_fname, lplt);
	if (!strcmp(pname, "none")) {
		if (hd_fd) close (hd_fd);
		hd_fd = 0;
		return;
	}
	hd_lon = 1;
	if (*pname == '\0') {
		strcpy (hd_filename, "niceplot.ps");
	} else {
		strcpy (hd_filename, pname);
	}
	if (hd_fd) close(hd_fd);
	hd_fd = open(hd_filename,(O_WRONLY|O_CREAT|O_TRUNC),0666);
	if (hd_fd<0) {
		fprintf(stderr,"hdopen: Error opening file %s.\n",hd_filename);
		exit (1);
	}
	lseek(hd_fd,0,0);
	hd_lstart = 1;
}

char *
convert_string (string, len)

char *string;
int len;

{
	static char out[512];
	int i, j, k;

	out[0] = '\0';
	if (len < 1) return (out);
	for (i=0; i<len; i++) if (string[i] != ' ') break;
	i=0;
	if (i == len) return (out);
	for (j=len; j>0; j--) if (string[j-1] != ' ') break;
	if (j == 0) return (out);
	for (k=0; i<j; i++) out[k++] = string[i];
	out[k] = '\0';
	return (out);
}
