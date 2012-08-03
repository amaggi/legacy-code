/*
 ******************************************************************************
 *
 *	Function hdgptr_
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

#include <stdio.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>

Display *hd_display;
Window hd_window;
Pixmap hd_pixmap;
Cursor hd_cursor;
XFontStruct *hd_cfont;
GC hd_gc;
int hd_width;
int hd_height;
void (*hd_cursorcallback)() = NULL;
void (*hd_cursortrack)() = NULL;

void
hdgptr_(itype, iwidth, iheight, ixmin, ixmax, xscale, xmin, 
	iymin, iymax, yscale, ymin, ix, iy, c, lc)

long *itype,*iwidth,*iheight;
long *ixmin,*iymin;
long *ixmax,*iymax;
float *xscale, *xmin;
float *yscale, *ymin;
long *ix,*iy;
char *c;
int lc;

/*
 *	hdgptr_ will get the current pointer location after the user has
 *	positioned the pointer and typed a single keyboard character.
 */

{
	Window root, child;
	int root_x, root_y, win_x, win_y;
	int win1_x, win1_y, win1_w, win1_h;
	int win2_x, win2_y, win2_w, win2_h;
	unsigned int keys_buttons;
	XEvent event;
	int nbytes;
	char buf[32];
	KeySym keysym;
	XComposeStatus compose;
	int x,y;
	int loop = 1;
	int iwarp;
	static char string[256];
	int direction, ascent, descent;
	XCharStruct overall;
	int str_x, str_y, str_w, str_h;
	float fx, fy;
	int iw, ih;

	hdstrk_();
	XSync (hd_display, True);
	XSelectInput (hd_display, hd_window, 
					(KeyPressMask | PointerMotionMask));
	XDefineCursor (hd_display, hd_window, hd_cursor);
	iw = (*iwidth)/2;
	ih = (*iheight)/2;
	if (XQueryPointer (hd_display, hd_window, &root, &child,
		&root_x, &root_y, &win_x, &win_y, &keys_buttons) == True) {
		x = win_x;
		y = win_y;
		if (x-iw < *ixmin) {
			x = *ixmin+iw;
		}
		if (x+iw > *ixmax) {
			x = *ixmax-iw;
		}
		if (y-ih < *iymin) {
			y = *iymin+ih;
		}
		if (y+ih > *iymax) {
			y = *iymax-ih;
		}
		fx = (x - (*ixmin))*(*xscale) + (*xmin);
		fy = (y - (*iymin))*(*yscale) + (*ymin);
		if (hd_cursortrack) {
			(*hd_cursortrack) (&fx, &fy);
		}
		XUndefineCursor (hd_display, hd_window);
		XWarpPointer (hd_display, None,
					hd_window, 0,0,0,0, x, y);
		XDefineCursor (hd_display, hd_window, hd_cursor);
		switch (*itype) {
		default:
			XDrawLine (hd_display, hd_window, hd_gc, 
						*ixmin, y, *ixmax, y);
			XDrawLine (hd_display, hd_window, hd_gc, 
						x, *iymin, x, *iymax);
			win1_x = *ixmin;
			win1_y = y;
			win1_w = (*ixmax)-(*ixmin)+1;
			win1_h = 1;
			win2_x = x;
			win2_y = *iymin;
			win2_w = 1;
			win2_h = (*iymax)-(*iymin)+1;
			break;
		case 1:
			XDrawLine (hd_display, hd_window, hd_gc, 
						x-iw, *iymin, x-iw, *iymax);
			XDrawLine (hd_display, hd_window, hd_gc, 
						x+iw, *iymin, x+iw, *iymax);
			win1_x = x-iw;
			win1_y = *iymin;
			win1_w = 1;
			win1_h = (*iymax)-(*iymin)+1;
			win2_x = x+iw;
			win2_y = *iymin;
			win2_w = 1;
			win2_h = (*iymax)-(*iymin)+1;
			break;
		case 2:
			XDrawLine (hd_display, hd_window, hd_gc, 
						*ixmin, y-ih, *ixmax, y-ih);
			XDrawLine (hd_display, hd_window, hd_gc, 
						*ixmin, y+ih, *ixmax, y+ih);
			win1_x = *ixmin;
			win1_y = y-ih;
			win1_w = (*ixmax)-(*ixmin)+1;
			win1_h = 1;
			win2_x = *ixmin;
			win2_y = y+ih;
			win2_w = (*ixmax)-(*ixmin)+1;
			win2_h = 1;
			break;
		case 3:
			XDrawLine (hd_display, hd_window, hd_gc, 
						x-iw, y-ih, x-iw, y+ih);
			XDrawLine (hd_display, hd_window, hd_gc, 
						x+iw, y-ih, x+iw, y+ih);
			XDrawLine (hd_display, hd_window, hd_gc, 
						x-iw, y-ih, x+iw, y-ih);
			XDrawLine (hd_display, hd_window, hd_gc, 
						x-iw, y+ih, x+iw, y+ih);
			win1_x = x-iw;
			win1_y = y-ih;
			win1_w = 2*iw;
			win1_h = 2*ih;
			win2_w = 0;
			break;
		}
		win_x = x;
		win_y = y;
		sprintf (string, "%fx%f", fx, fy);
		XTextExtents (hd_cfont, string, strlen(string),
				&direction, &ascent, &descent, &overall);
		XSetFont (hd_display, hd_gc, hd_cfont->fid);
		str_x = 0;
		str_y = ascent;
		str_w = overall.width;
		str_h = ascent + descent;
		XDrawImageString (hd_display, hd_window, hd_gc,
				str_x, str_y, string, strlen(string));
	} else {
		XUndefineCursor (hd_display, hd_window);
		XSelectInput (hd_display, hd_window, PropertyChangeMask);
		return;
	}
	while (loop) {
		XNextEvent (hd_display, &event);
		if (event.xany.window != hd_window) continue;
		switch (event.type) {
		case KeyPress:
			if (event.xkey.x < 0) break;
			if (event.xkey.y < 0) break;
			nbytes = XLookupString (&event, buf, sizeof(buf),
							&keysym, &compose);
			if (nbytes == 1) {
				*c = buf[0];
				for (lc--, c++; lc > 0; lc--, c++) *c = ' ';
				*ix = event.xkey.x;
				*iy = event.xkey.y;
				loop = 0;
			} else {
				switch (keysym) {
				case XK_Left:
					x = event.xkey.x - 1;
					y = event.xkey.y;
					if (x-iw < *ixmin) x = *ixmin+iw;
					XUndefineCursor (hd_display, hd_window);
					XWarpPointer (hd_display, None,
						hd_window, 0,0,0,0, x, y);
					XDefineCursor (hd_display, hd_window, 
								hd_cursor);
					break;
				case XK_Right:
					x = event.xkey.x + 1;
					y = event.xkey.y;
					if (x+iw > *ixmax) x = *ixmax-iw;
					XUndefineCursor (hd_display, hd_window);
					XWarpPointer (hd_display, None,
						hd_window, 0,0,0,0, x, y);
					XDefineCursor (hd_display, hd_window, 
								hd_cursor);
					break;
				case XK_Up:
					x = event.xkey.x;
					y = event.xkey.y - 1;
					if (y-ih < *iymin) y = *iymin+ih;
					XUndefineCursor (hd_display, hd_window);
					XWarpPointer (hd_display, None,
						hd_window, 0,0,0,0, x, y);
					XDefineCursor (hd_display, hd_window, 
								hd_cursor);
					break;
				case XK_Down:
					x = event.xkey.x;
					y = event.xkey.y + 1;
					if (y+ih > *iymax) y = *iymax-ih;
					XUndefineCursor (hd_display, hd_window);
					XWarpPointer (hd_display, None,
						hd_window, 0,0,0,0, x, y);
					XDefineCursor (hd_display, hd_window, 
								hd_cursor);
					break;
				default:
					break;
				}
			}
			break;
		case MotionNotify:
			iwarp = 0;
			x = event.xmotion.x;
			y = event.xmotion.y;
			if (x-iw < *ixmin) {
				iwarp = 1;
				x = *ixmin+iw;
			}
			if (x+iw > *ixmax) {
				iwarp = 1;
				x = *ixmax-iw;
			}
			if (y-ih < *iymin) {
				iwarp = 1;
				y = *iymin+ih;
			}
			if (y+ih > *iymax) {
				iwarp = 1;
				y = *iymax-ih;
			}
			if (iwarp) {
				XUndefineCursor (hd_display, hd_window);
				XWarpPointer (hd_display, None,
						hd_window, 0,0,0,0, x, y);
				XDefineCursor (hd_display, hd_window, 
								hd_cursor);
				break;
			}
			XCopyArea (hd_display, hd_pixmap, hd_window, hd_gc,
				win1_x, win1_y, win1_w, win1_h, win1_x, win1_y);
			if (win2_w)
			XCopyArea (hd_display, hd_pixmap, hd_window, hd_gc,
				win2_x, win2_y, win2_w, win2_h, win2_x, win2_y);
			XCopyArea (hd_display, hd_pixmap, hd_window, hd_gc,
					0, 0, str_w, str_h, 0, 0);
			fx = (x - (*ixmin))*(*xscale) + (*xmin);
			fy = (y - (*iymin))*(*yscale) + (*ymin);
			if (hd_cursortrack) {
				(*hd_cursortrack) (&fx, &fy);
			}
			switch (*itype) {
			default:
				XDrawLine (hd_display, hd_window, hd_gc, 
						*ixmin, y, *ixmax, y);
				XDrawLine (hd_display, hd_window, hd_gc, 
						x, *iymin, x, *iymax);
				win1_x = *ixmin;
				win1_y = y;
				win1_w = (*ixmax)-(*ixmin)+1;
				win1_h = 1;
				win2_x = x;
				win2_y = *iymin;
				win2_w = 1;
				win2_h = (*iymax)-(*iymin)+1;
				break;
			case 1:
				XDrawLine (hd_display, hd_window, hd_gc, 
						x-iw, *iymin, x-iw, *iymax);
				XDrawLine (hd_display, hd_window, hd_gc, 
						x+iw, *iymin, x+iw, *iymax);
				win1_x = x-iw;
				win1_y = *iymin;
				win1_w = 1;
				win1_h = (*iymax)-(*iymin)+1;
				win2_x = x+iw;
				win2_y = *iymin;
				win2_w = 1;
				win2_h = (*iymax)-(*iymin)+1;
				break;
			case 2:
				XDrawLine (hd_display, hd_window, hd_gc, 
						*ixmin, y-ih, *ixmax, y-ih);
				XDrawLine (hd_display, hd_window, hd_gc, 
						*ixmin, y+ih, *ixmax, y+ih);
				win1_x = *ixmin;
				win1_y = y-ih;
				win1_w = (*ixmax)-(*ixmin)+1;
				win1_h = 1;
				win2_x = *ixmin;
				win2_y = y+ih;
				win2_w = (*ixmax)-(*ixmin)+1;
				win2_h = 1;
				break;
			case 3:
				XDrawLine (hd_display, hd_window, hd_gc, 
						x-iw, y-ih, x-iw, y+ih);
				XDrawLine (hd_display, hd_window, hd_gc, 
						x+iw, y-ih, x+iw, y+ih);
				XDrawLine (hd_display, hd_window, hd_gc, 
						x-iw, y-ih, x+iw, y-ih);
				XDrawLine (hd_display, hd_window, hd_gc, 
						x-iw, y+ih, x+iw, y+ih);
				win1_x = x-iw;
				win1_y = y-ih;
				win1_w = 2*iw;
				win1_h = 2*ih;
				win2_w = 0;
				break;
			}
			win_x = x;
			win_y = y;
			sprintf (string, "%fx%f", fx, fy);
			XTextExtents (hd_cfont, string, strlen(string),
				&direction, &ascent, &descent, &overall);
			XSetFont (hd_display, hd_gc, hd_cfont->fid);
			str_x = 0;
			str_y = ascent;
			str_w = overall.width;
			str_h = ascent + descent;
			XDrawImageString (hd_display, hd_window, hd_gc,
				str_x, str_y, string, strlen(string));
			break;
		default:
			break;
		}
	}
	XCopyArea (hd_display, hd_pixmap, hd_window, hd_gc,
				win1_x, win1_y, win1_w, win1_h, win1_x, win1_y);
	if (win2_w) XCopyArea (hd_display, hd_pixmap, hd_window, hd_gc,
				win2_x, win2_y, win2_w, win2_h, win2_x, win2_y);
	XCopyArea (hd_display, hd_pixmap, hd_window, hd_gc,
					0, 0, str_w, str_h, 0, 0);
	XUndefineCursor (hd_display, hd_window);
	XSelectInput (hd_display, hd_window, PropertyChangeMask);
	XSync (hd_display, False);
}

void
hdtptr_(callback)

void (*callback)();

/*
 *	hdtptr_ will set the cursor tracking callback proceedure.
 */

{
	hd_cursortrack = callback;
}

void
hdtptroff_()

/*
 *	hdtptroff_ will turn off cursor tracking callback.
 */

{
	hd_cursortrack = NULL;
}
