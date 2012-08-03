/*
 ******************************************************************************
 *
 *	Function hdbg and hdfg
 *
 *	Author - Ben Kohl
 *
 ******************************************************************************
 */

#include <stdio.h>
#include <X11/Xlib.h>

#ifdef STELLAR
#include <xfdi.h>
#endif

Display *hd_display;
Colormap hd_cmap;
Pixmap hd_pixmap;
GC hd_gc;
Visual *hd_visual;
unsigned long int hd_fore;
unsigned long int hd_back;
#ifdef STELLAR
XFDIGC hd_gc3;
#endif

char *convert_string();

int
hdbg_( colorname, h, l, s, length )

char *colorname;
float *h, *l, *s;
int length;

/*
 *	hdbg will set the background color.
 */

{
	XColor color;
	float r, g, b;
	char cname[256];

	strcpy( cname, convert_string(colorname, length) );

	if ( !XParseColor( hd_display, hd_cmap, cname, &color ) ) {
	  fprintf( stderr, "Error parsing named color:  %s\n", cname );
	  return(1);
	}

	if ( !XAllocColor (hd_display, hd_cmap, &color) ) {
	  fprintf( stderr, "Error allocating named color:  %s\n", cname );
	  return(1);
	}

	r = color.red / 65535.0;
	g = color.green / 65535.0;
	b = color.blue / 65535.0;
	hdrgb_to_hls( r, g, b, h, l, s );

	hd_back = color.pixel;

	XSetBackground (hd_display, hd_gc, hd_back);

	return(0);
}

int
hdfg_( colorname, h, l, s, length )

char *colorname;
float *h, *l, *s;
int length;

/*
 *	hdfg will set the foreground color.
 */

{
	XColor color;
	float r, g, b;
	char cname[256];

	strcpy( cname, convert_string(colorname, length) );

	if ( !XParseColor( hd_display, hd_cmap, cname, &color ) ) {
	  fprintf( stderr, "Error parsing named color:  %s\n", cname );
	  return(1);
	}

	if ( !XAllocColor (hd_display, hd_cmap, &color) ) {
	  fprintf( stderr, "Error allocating named color:  %s\n", cname );
	  return(1);
	}

	r = color.red / 65535.0;
	g = color.green / 65535.0;
	b = color.blue / 65535.0;
	hdrgb_to_hls( r, g, b, h, l, s );

	hd_fore = color.pixel;

	XSetForeground (hd_display, hd_gc, hd_fore);

	return(0);
}
