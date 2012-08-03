/*
 *******************************************************************************
 *
 *	nxplotd_getatoms
 *
 *	Author - Danny Harvey
 *
 *******************************************************************************
 */

#include <stdio.h>

#include "nxplotdlib_types.h"

int
nxplotd_getatoms (xplptr)

nxplotdLib *xplptr;

{
	/* get all of the window property atoms */
	xplptr->window_atom = XInternAtom (xplptr->display, 
						NXPLOTD_WINDOW_PROP, False);
	if (xplptr->window_atom == None) {
		fprintf (stderr, "%s: Unable to get %s atom.\n",
					NXPLOTD_WINDOW_PROP, xplptr->progname);
		exit (1);
	}
	xplptr->daemon_window_atom = XInternAtom (xplptr->display, 
					NXPLOTD_DAEMON_WINDOW_PROP, False);
	if (xplptr->daemon_window_atom == None) {
		fprintf (stderr, "%s: Unable to get %s atom.\n",
				NXPLOTD_DAEMON_WINDOW_PROP, xplptr->progname);
		exit (1);
	}
	xplptr->client_status_atom = XInternAtom (xplptr->display, 
					NXPLOTD_CLIENT_STATUS_PROP, False);
	if (xplptr->client_status_atom == None) {
		fprintf (stderr, "%s: Unable to get %s atom.\n",
				NXPLOTD_CLIENT_STATUS_PROP, xplptr->progname);
		exit (1);
	}
	xplptr->daemon_status_atom = XInternAtom (xplptr->display, 
					NXPLOTD_DAEMON_STATUS_PROP, False);
	if (xplptr->daemon_status_atom == None) {
		fprintf (stderr, "%s: Unable to get %s atom.\n",
				NXPLOTD_DAEMON_STATUS_PROP, xplptr->progname);
		exit (1);
	}
	xplptr->client_request_atom = XInternAtom (xplptr->display, 
					NXPLOTD_CLIENT_REQUEST_PROP, False);
	if (xplptr->client_request_atom == None) {
		fprintf (stderr, "%s: Unable to get %s atom.\n",
				NXPLOTD_CLIENT_REQUEST_PROP, xplptr->progname);
		exit (1);
	}
	xplptr->daemon_pixmap_atom = XInternAtom (xplptr->display, 
					NXPLOTD_DAEMON_PIXMAP_PROP, False);
	if (xplptr->daemon_pixmap_atom == None) {
		fprintf (stderr, "%s: Unable to get %s atom.\n",
				NXPLOTD_DAEMON_PIXMAP_PROP, xplptr->progname);
		exit (1);
	}
}

/*
 *******************************************************************************
 *
 *	nxplotd_getprop
 *
 *	Author - Danny Harvey
 *
 *******************************************************************************
 */

int
nxplotd_getprop (program, display, window, prop_string, atom, type, prop)

char *program;
Display *display;
Window window;
char *prop_string;
Atom atom;
Atom type;
unsigned char **prop;

{
	Status status;
	Atom actual_type;
	int actual_format;
	unsigned long nitems, bytes_after;

	status = XGetWindowProperty (display, window, atom,
				0, 1, False, type, &actual_type,
				&actual_format, &nitems, &bytes_after,
				prop);
	switch (status) {
	case Success:
		if (actual_type == None) {
			fprintf (stderr, "%s: %s does not exist.\n",
							program, prop_string);
			exit (1);
		}
		if (actual_type == type) {
			break;
		}
		fprintf (stderr, "%s: %s type mismatch.\n", 
							program, prop_string);
		exit (1);
	case BadWindow:
		fprintf (stderr, "%s: Bad window.\n", program);
		exit (1);
	default:
		fprintf (stderr, "%s: Unknown XGetWindowProperty error (%d).\n",
						program, status);
		exit (1);
	}
}
