/*
 *******************************************************************************
 *
 *	nxplotd_launchdaemon
 *
 *	Author - Danny Harvey
 *
 *******************************************************************************
 */

#include "nxplotdlib_types.h"
#include "my_system.h"

#include <stdio.h>

int
nxplotd_launchdaemon (program, display, window, args, xplptr)

char *program;
Display *display;
Window window;
char *args;
nxplotdLib *xplptr;

{
	static Window win_null = NULL;
	int status;
	char cmd[256];
	int coreflag;
	int termsig;
	XEvent report;
	int loop;
	Window *window_ptr;
	int *int_ptr;
	Pixmap *pixmap_ptr;
	Colormap *colormap_ptr;
	int i;

/*
 *	Initialize the nxplotLib structure by getting the window property atoms
 *	from the X11 server.
 */
	xplptr->progname = program;
	xplptr->display = display;
	nxplotd_getatoms (xplptr);
/*
 *	Initialize the properties for this window.
 */
	XChangeProperty (display, window, xplptr->window_atom, XA_WINDOW,
				32, PropModeReplace, &window, 1);
	XChangeProperty (display, window, xplptr->daemon_window_atom, XA_WINDOW,
				32, PropModeReplace, &win_null, 1);
	XChangeProperty (display, window, xplptr->daemon_pixmap_atom, XA_PIXMAP,
				32, PropModeReplace, &win_null, 1);
	status = STATUS_INITIAL;
	XChangeProperty (display, window, xplptr->client_status_atom,XA_INTEGER,
				32, PropModeReplace, &status, 1);
	XChangeProperty (display, window, xplptr->daemon_status_atom,XA_INTEGER,
				32, PropModeReplace, &status, 1);
	xplptr->client_window = window;
	xplptr->daemon_window = NULL;
	xplptr->daemon_pixmap = NULL;
/*
 *	Set the event mask to accept handshake from the daemon.
 */
	XSelectInput(display, window, PropertyChangeMask); 
	XSync (display, False);
/*
 *	Launch the daemon program.
 */
	if (args) {
		sprintf (cmd, "nxplotd -window %d -title %s %s &", 
							window, program, args);
	} else {
		sprintf (cmd, "nxplotd -window %d -title %s &", window, program);
	}
/*	printf ("%s[libgrx]: Launching plotting daemon...\n", program);*/
	xplptr->daemon_pid = my_system (cmd, 1, &coreflag, &termsig);
	if (xplptr->daemon_pid < 0) {
		switch (xplptr->daemon_pid) {
		case MY_SYSTEM_COMMAND_NOEXIST:
			fprintf (stderr, 
			"%s[nxplotd_launchdaemon]: nxplotd does not exist.\n",
								program);
			exit (1);
		default:
			fprintf (stderr, 
			"%s[nxplotd_launchdaemon]: Unable to launch nxplotd.\n",
								program);
			exit (1);
		}
	}
/*
 *	Look for the daemon window ID.
 */
	loop = 1;
/*	printf ("%s: looking for daemon window ID...\n", program);*/
	while (loop) {
		XNextEvent(display, &report);
		switch  (report.type) {
		case PropertyNotify:
			if (report.xproperty.atom ==xplptr->daemon_window_atom){
				nxplotd_getprop (program,
						display, window, 
						NXPLOTD_DAEMON_WINDOW_PROP, 
						report.xproperty.atom, 
						XA_WINDOW, &window_ptr);
				if (*window_ptr) loop = 0;
			}
			break;
		default:
			break;
		}
	}
/*	printf ("%s: daemon window ID = 0x%x\n", program, *window_ptr);*/
	xplptr->daemon_window = *window_ptr;
/*
 *	Change client status to ready for both the client and daemon windows.
 */
	status = STATUS_READY;
/*	printf ("%s: setting STATUS_READY (%d).\n", program, status);*/
	sleep (1);
	XChangeProperty (display, xplptr->daemon_window, 
			xplptr->client_status_atom, XA_INTEGER, 32, 
			PropModeReplace, &status, 1);
	XChangeProperty (display, xplptr->client_window, 
			xplptr->client_status_atom, XA_INTEGER, 32, 
			PropModeReplace, &status, 1);
/*
 *	Wait for daemon ready.
 */
	/* Wait for daemon ready */
	loop = 1;
/*	printf ("%s: waiting for daemon ready...\n", program);*/
	while (loop) {
		XNextEvent(display, &report);
		switch  (report.type) {
		case PropertyNotify:
			if (report.xproperty.atom ==xplptr->daemon_status_atom){
				nxplotd_getprop (program,
						display, window, 
						NXPLOTD_DAEMON_STATUS_PROP, 
						report.xproperty.atom, 
						XA_INTEGER, &int_ptr);
				if (*int_ptr == STATUS_READY) loop = 0;
			}
			break;
		default:
			break;
		}
	}
/*	printf ("%s[libgrx]: Plotting daemon ready.\n", program);*/
/*
 *	Get the daemon pixmap and graphics context.
 */
	nxplotd_getprop (program, display, xplptr->client_window, 
			NXPLOTD_DAEMON_PIXMAP_PROP, xplptr->daemon_pixmap_atom, 
			XA_PIXMAP, &pixmap_ptr);
	if (*pixmap_ptr == NULL) {
		fprintf (stderr, 
			"%s[nxplotd_launchdaemon]: Invalid daemon pixmap.\n", 
								program);
	}
	xplptr->daemon_pixmap = *pixmap_ptr;
/*	printf ("%s: daemon pixmap ID = 0x%x\n", program, *pixmap_ptr);*/
/*
 *	Normal exit.
 */
	return (1);
}
