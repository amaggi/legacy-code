/*
 *	This include file (nxplotdlib_types.h) contains all of the typedefs
 *	for the nxplotd library procedures.
 */

#include <X11/Xlib.h>
#include <X11/Xatom.h>

#define	NXPLOTD_WINDOW_PROP		"NXPLOTD_WINDOW_PROP"
#define	NXPLOTD_DAEMON_WINDOW_PROP	"NXPLOTD_DAEMON_WINDOW_PROP"
#define	NXPLOTD_CLIENT_STATUS_PROP	"NXPLOTD_CLIENT_STATUS_PROP"
#define	NXPLOTD_DAEMON_STATUS_PROP	"NXPLOTD_DAEMON_STATUS_PROP"
#define	NXPLOTD_CLIENT_REQUEST_PROP	"NXPLOTD_CLIENT_REQUEST_PROP"
#define	NXPLOTD_DAEMON_PIXMAP_PROP	"NXPLOTD_DAEMON_PIXMAP_PROP"

#define	STATUS_INITIAL	0
#define	STATUS_READY	1

#define	REQUEST_NULL	0
#define	REQUEST_COPY	1
#define	REQUEST_QUIT	2

typedef struct xplotd_lib_ {	/* nxplotdLib holds all of the state info
				   for the nxplotd library routines. */
	char *progname;			/* Program name. */
	Display *display;		/* Program display. */
	Atom window_atom;		/* Window ID atom for this window. */
	Atom daemon_window_atom;	/* Window ID atom for daemon window. */
	Atom client_status_atom;	/* Client status atom. */
	Atom daemon_status_atom;	/* Daemon status atom. */
	Atom client_request_atom;	/* Client request atom. */
	Atom daemon_pixmap_atom;	/* Pixmap ID atom for daemon pixmap. */
	int daemon_pid;			/* Daemon process id. */
	Window client_window;		/* Client window ID. */
	Window daemon_window;		/* Daemon window ID. */
	Pixmap daemon_pixmap;		/* Daemon Pixmap ID. */
} nxplotdLib;
