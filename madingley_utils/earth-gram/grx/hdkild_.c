/*
 ******************************************************************************
 *
 *	Function hdkild_
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

#include <X11/Xlib.h>

#include "nxplotdlib_types.h"

Display *hd_display;
Window hd_window;
nxplotdLib hd_xpllib;

void
hdkild_()

/*
 *	hdkild_ will kill the niceplot plotting daemon program.
 */

{
	int request;

	request = REQUEST_QUIT;
	XChangeProperty (hd_display, hd_xpllib.daemon_window, 
				hd_xpllib.client_request_atom,
				XA_INTEGER, 32, PropModeReplace, &request, 1);
	XFlush (hd_display);
}
