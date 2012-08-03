/*
 ******************************************************************************
 *
 *	Function hdinit_
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

#include <stdio.h>
#include <time.h>
#include <sys/types.h>
#ifdef NeXT
#include <sys/time.h>
#endif

int hd_fd;
int hd_npts;
int hd_nptsl;

void
hdinit_(itran, fname, prname, fname_len, prname_len)

int *itran;
char *fname;
char *prname;
int fname_len;
int prname_len;

/*
 *	hdinit_ will initialize postscript for plotting
 */

{
	static char line[256];
	static char filn[256];
	static char prgn[256];
	static int ninit = 3;
	static char *init[] = {
		"initmatrix newpath initclip\n",
		"1 setlinewidth 0 setlinecap 0 setlinejoin\n",
		"[] 0 setdash 0 setgray 10 setmiterlimit\n"
	};
	static int npro = 80;
	static char *prolog[] = {
		"/npchsc 1.0 def\n",
		"/m {moveto} def\n",
		"/l {lineto} def\n",
		"/r {rlineto} def\n",
		"/x {0 rlineto} def\n",
		"/y {0 exch rlineto} def\n",
		"/A {-2 2 rlineto} def\n",
		"/B {-1 2 rlineto} def\n",
		"/C {0 2 rlineto} def\n",
		"/D {1 2 rlineto} def\n",
		"/E {2 2 rlineto} def\n",
		"/F {-2 1 rlineto} def\n",
		"/G {-1 1 rlineto} def\n",
		"/H {0 1 rlineto} def\n",
		"/I {1 1 rlineto} def\n",
		"/J {2 1 rlineto} def\n",
		"/K {-2 0 rlineto} def\n",
		"/L {-1 0 rlineto} def\n",
		"/M {1 0 rlineto} def\n",
		"/N {2 0 rlineto} def\n",
		"/O {-2 -1 rlineto} def\n",
		"/P {-1 -1 rlineto} def\n",
		"/Q {0 -1 rlineto} def\n",
		"/R {1 -1 rlineto} def\n",
		"/S {2 -1 rlineto} def\n",
		"/T {-2 -2 rlineto} def\n",
		"/U {-1 -2 rlineto} def\n",
		"/V {0 -2 rlineto} def\n",
		"/W {1 -2 rlineto} def\n",
		"/X {2 -2 rlineto} def\n",
		"/npfont112 {/font /Helvetica def} def\n",
		"/npfont115 {/font /Times-Roman def} def\n",
		"/npfont116 {/font /Times-Italic def} def\n",
		"/npfont130 {/font /Palatino-Roman def} def\n",
		"/npfont131 {/font /Palatino-Italic def} def\n",
		"/npfont132 {/font /Palatino-Bold def} def\n",
#ifdef NeXT
	"/nptext { %xp yp size szsc angle ref xl xr yb yt clip mytext nptext\n",
		"/mytext exch def\n",
		"/clip exch def\n",
#else
	"/nptext { %xp yp size szsc angle ref xl xr yb yt clp mytext nptext\n",
		"/mytext exch def\n",
		"/clp exch def\n",
#endif
		"/yt exch def\n",
		"/yb exch def\n",
		"/xr exch def\n",
		"/xl exch def\n",
		"/ref exch def\n",
		"/angle exch def\n",
		"/szsc exch def\n",
		"/size exch def\n",
		"/yp exch def\n",
		"/xp exch def\n",
		"ref 0 lt {/ref 0 def} if\n",
		"ref 8 gt {/ref 0 def} if\n",
		"size 1.5 mul szsc mul /size exch def\n",
		"ref cvi /ref exch def\n",
		"ref 0.5 add 3 div cvi /xjust exch def\n",
		"ref xjust 3 mul neg add /yjust exch def\n",
		"font findfont size scalefont setfont\n",
		"mytext stringwidth /wy exch def /wx exch def\n",
		"wx wx mul wy wy mul add sqrt /width exch def\n",
		"width 0.5 mul xjust mul /xj exch def\n",
		"size 0.333 mul yjust mul /yj exch def\n",
		"xj angle cos mul yj angle sin mul sub /xpp exch def\n",
		"xj angle sin mul yj angle cos mul add /ypp exch def\n",
		"xp xpp sub /xp exch def\n",
		"yp ypp sub /yp exch def\n",
		"angle rotate\n",
		"xp angle cos mul yp angle sin mul add /xpp exch def\n",
		"xp neg angle sin mul yp angle cos mul add /ypp exch def\n",
		"xpp ypp moveto\n",
		"mytext show\n",
		"angle neg rotate\n",
		"} def\n",
		"/rectcl {\n",
		"/ht exch def\n",
		"/wd exch def\n",
		"/yc exch def\n",
		"/xc exch def\n",
		"newpath xc yc moveto xc wd add yc lineto\n",
		"xc wd add yc ht add lineto xc yc ht add lineto\n",
		"closepath clip newpath\n",
		"} def\n"
	};
	int i;
	int llx, lly, urx, ury;

	char *convert_string ();

	if (hd_fd > 0) {
		sprintf (line, "%%!PS-Adobe-3.0\n");
		write (hd_fd, line, strlen(line));
		strcpy (prgn, convert_string(prname, prname_len));
		strcpy (filn, convert_string(fname, fname_len));
		if (!strcmp(filn, "")) strcpy (filn, "niceplot.ps");
		pscreate (hd_fd, prgn, filn);
		sprintf (line, "%%%%EndComments\n");
		write (hd_fd, line, strlen(line));
		for (i=0; i<ninit; i++)
			write (hd_fd, init[i], strlen(init[i]));
		if (*itran == 1) {
			sprintf (line, "603.5 0 translate 90 rotate\n");
			write (hd_fd, line, strlen(line));
#ifdef NeXT
			llx = 0.5*72.0;
			lly = 0.5*72.0;
			urx = 10.5*72.0;
			ury = 8.0*72.0;
		} else {
			llx = 0.5*72.0;
			lly = 0.5*72.0;
			ury = 10.5*72.0;
			urx = 8.0*72.0;
#else
			llx = 0;
			lly = 0;
			urx = 10.0*72.0;
			ury = 7.5*72.0;
		} else {
			llx = 0;
			lly = 0;
			ury = 10.0*72.0;
			urx = 7.5*72.0;
#endif
		}
		sprintf (line, "%%%%BeginDocument: %s\n", filn);
		write (hd_fd, line, strlen(line));
		epsfpro (hd_fd, llx, lly, urx, ury, prgn, filn);
		sprintf (line, "%%%%BeginProlog\n");
		write (hd_fd, line, strlen(line));
		sprintf (line, 
		"%% Start of niceplot.pro -- prolog for niceplot PostScript\n");
		write (hd_fd, line, strlen(line));
		for (i=0; i<npro; i++)
			write (hd_fd, prolog[i], strlen(prolog[i]));
		sprintf (line, "%%%%EndProlog\n");
		write (hd_fd, line, strlen(line));
		sprintf (line, ".236667 .236667 scale\n");
		write (hd_fd, line, strlen(line));
		sprintf (line, "0 setlinewidth\n");
		write (hd_fd, line, strlen(line));
		sprintf (line, "newpath\n");
		write (hd_fd, line, strlen(line));
	}
	hd_npts = 0;
	hd_nptsl = 0;
}

pscreate (fd, prname, fname)

int fd;
char *prname;
char *fname;

{
	static char line[256];
	static char host[32];
	char *user;
	char *login;
#ifdef NeXT
	time_t tloc;
#else
	int tloc;
#endif

	char *getlogin();

	if (gethostname (host, 32) < 0) strcpy (host, "NULL");
	login = getlogin();
#ifdef NeXT
	user = login;
#else
	user = cuserid(NULL);
#endif
	if (user && login)
		sprintf (line, "%%%%Creator: %s:%s (%s)\n", host, user, login);
	else if (user)
		sprintf (line, "%%%%Creator: %s:%s (NULL)\n", host, user);
	else if (login)
		sprintf (line, "%%%%Creator: %s:NULL (%s)\n", host, login);
	else sprintf (line, "%%%%Creator: %s:NULL (NULL)\n", host);
	write (fd, line, strlen(line));
	sprintf (line, "%%%%Title: %s (%s)\n", fname, prname);
	write (fd, line, strlen(line));
	time (&tloc);
	sprintf (line, "%%%%CreationDate: %s", (char *)ctime(&tloc));
	write (fd, line, strlen(line));
}

epsfpro (fd, llx, lly, urx, ury, prname, fname)

int fd;
int llx;
int lly;
int urx;
int ury;
char *prname;
char *fname;

{
	static char line[256];

	sprintf (line, "%%!PS-Adobe-3.0 EPSF-3.0\n");
	write (fd, line, strlen(line));
	sprintf (line, "%%%%BoundingBox: %d %d %d %d\n", llx, lly, urx, ury);
	write (fd, line, strlen(line));
	pscreate (fd, prname, fname);
	sprintf (line, "%%%%EndComments\n");
	write (fd, line, strlen(line));
}

void
hdcmnt_ (comment, comment_len)

char * comment;
int     comment_len;

{
	static char line[256];

	char *convert_string ();

	if (hd_fd > 0) {
		sprintf (line, "%%%s\n", convert_string(comment, comment_len));
		write (hd_fd, line, strlen(line));
	}
}
