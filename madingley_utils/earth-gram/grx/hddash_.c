/*
 ******************************************************************************
 *
 *	Function hddash_
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

#include <X11/Xlib.h>

Display *hd_display;
int hd_fd;
int hd_lon;
GC hd_gc;

void
hddash_(n, pat)

long *n, *pat;

/*
 *	hddash_ will set the line style
 */

{
	void hddashl();

	if (*n < 1) {
		XSetLineAttributes (hd_display, hd_gc, 0, LineSolid,
				CapButt, JoinMiter ) ;
	} else {
		XSetLineAttributes (hd_display, hd_gc, 0, LineOnOffDash,
				CapButt, JoinMiter ) ;
	}
	if (hd_fd > 0 && hd_lon) hddashl (n, pat);
}

void
hddashl(n, pat)

long *n, *pat;

{
	char line[80];
	int i,j;

	if (*n < 1) {
		write (hd_fd, "[] 0 setdash\n", 13);
	} else {
		line[0] = '[';
		line[1] = '\0';
		for (i=0; i<*n; i++) {
			j = strlen(line);
			sprintf(&line[j], "%d \0", pat[i]);
		}
		j = strlen(line);
		j--;
		line[j] = '\0';
		strcat (line, "] 0 setdash\n");
		write (hd_fd, line, strlen(line));
	}
}
