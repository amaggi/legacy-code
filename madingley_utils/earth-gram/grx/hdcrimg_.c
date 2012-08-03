/*
 *******************************************************************************
 *
 *	FORTRAN callable subroutine hdcrimg
 *
 *	Author - Danny Harvey
 *
 *******************************************************************************
 */

#include <stdio.h>
#include <X11/Xlib.h>

Display *hd_display;
int hd_screen;
Visual *hd_visual;
Pixmap hd_pixmap;
Colormap hd_cmap;
Window hd_window;
GC hd_gc;
int hd_batch;
int hd_depth;
Colormap hd_cmap;

typedef struct image_ {
	int width;
	int height;
	int depth;
	int ctable[256];
	int red[256];
	int green[256];
	int blue[256];
	char *image;
} Image;

int hd_rv = 0;
int hd_pscolor = 0;

int hd_fd;
int hd_lon;

void hdpscoloron_ ();
void hdpscoloroff_ ();

void
hdpscolorset_ (i)

int *i;

{
	switch (*i) {
	case 1:
		hdpscoloron_ ();
		break;
	case 2:
	case 3:
		hd_pscolor = (*i);
		break;
	default:
		hdpscoloroff_ ();
		break;
	}
}

void
hdpscoloron_ ()

{
	hd_pscolor = 1;
}

void
hdpscoloroff_ ()

{
	static char line[128];

	hd_pscolor = 0;
	if (hd_fd < 1) return;
	if (!hd_lon) return;
	sprintf (line, "0 setgray\n");
	write (hd_fd, line, strlen(line));
}

void
hdcrimg_ (width, height, depth, image)

int *width;
int *height;
int *depth;
Image **image;

/*
 *	hdcrimg will create an 8-bit deep image for use in later drawing
 *	and painting routines.
 *
 *	Inputs  -	width	= Width of image in pixels.
 *			height	= Height of image in pixels.
 *			depth	= Depth of image as an integer count. This
 *				  should be in the range of 1 to 255.
 */

{
	Image *im;
	int size;
	char *imm;
	int i, j;
	float h, l, s;
	float r, g, b;
	float x;

	unsigned long hdhls_to_pixel();

	im = NULL;
	*image = im;
	if (*width < 1 || *height < 1) return;
	size = (*width)*(*height);
	imm = (char *) malloc (size);
	if (imm == NULL) return;
	im = (Image *) malloc (sizeof(Image));
	if (im == NULL) {
		free (imm);
		return;
	}
	im->width = *width;
	im->height = *height;
	im->depth = *depth;
	if (im->depth < 1) im->depth = 1;
	if (im->depth > 256) im->depth = 256;
	im->image = imm;
	for (i=0; i<size; i++) imm[i] = 0;
	for (i=0; i<im->depth; i++) {
/* bck commented out
		if (hd_depth < 2) {
			if (i%2) {
			im->ctable[i] = WhitePixel (hd_display, hd_screen);
			} else {
			im->ctable[i] = BlackPixel (hd_display, hd_screen);
			}
		} else {
 */
			x = (im->depth-1);
			x = i / x;
			x *= 128.0;
			j = x + 0.5;
			x = j / 128.0;
			l = 0.3*(1.0 - x) + 0.5;
			s = 1.0*(x) + 0.0;
			h = 240.0 * (1.0 - x);
			s = 1.0;
			im->ctable[i] = hdhls_to_pixel (h, l, s);
			h = (im->depth-1);
			l = 0.1*(1.0 - i / h) + 0.5;
			s = 1.0;
			h = 240.0 * (1.0 - i / h);
			hdhls_to_rgb (h, l, s, &r, &g, &b);
			im->red[i] = r*255.9999;
			im->green[i] = g*255.9999;
			im->blue[i] = b*255.9999;
/* bck commented out
		}
 */
	}
	*image = im;
}

/*
 *******************************************************************************
 *
 *	FORTRAN callable subroutine hdfrimg
 *
 *	Author - Danny Harvey
 *
 *******************************************************************************
 */

void
hdfrimg_ (image)

Image **image;

/*
 *	hdfrimg will free up an image resource.
 */

{
	if (*image) {
		if ( (*image)->image ) free ( (*image)->image );
		free (*image);
	}
}

/*
 *******************************************************************************
 *
 *	FORTRAN callable subroutine hdinimg
 *
 *	Author - Danny Harvey
 *
 *******************************************************************************
 */

void
hdinimg_ (image, nx, nxmax, ny, x, y, z, xmin, xmax, ymin, ymax, zmin, dz)

Image **image;
int *nx;
int *ny;
int *nxmax;
float *x;
float *y;
float *z;
float *xmin;
float *xmax;
float *ymin;
float *ymax;
float *zmin;
float *dz;

/*
 *	hdinimg will interpolate a sampled x-y grid of points to an image
 *	pixel mesh.
 *
 *	Inputs -	image	= Image pointer.
 *			nx	= No. of x grid nodes.
 *			ny	= No. of y grid nodes.
 *			x(nx)	= X grid coordinates.
 *			y(ny)	= Y grid coordinates.
 *			z(nx*ny)	
 *				= Dependent value, z = z(x,y), which is
 *				  to be interpolated. 
 *			xmin	= X-value which maps to the left edge of the
 *				  pixel image.
 *			xmax	= X-value which maps to the right edge of the
 *				  pixel image.
 *			ymin	= Y-value which maps to the bottom edge of the
 *				  pixel image.
 *			ymax	= Y-value which maps to the top edge of the
 *				  pixel image.
 *			zmin	= Z-value which maps to pixel value 0.
 *			dz	= Z increment which maps to a pixel change of 1.
 */

{
	Image *im;
	int nxmx;
	int i, j, k;
	int y0, y1, x0, x1;
	float xscale, yscale, zscale;
	char *buf;
	float xint, yint, zint;
	int zpix;
	float fy0, fy1, fx0, fx1;
	float fz00, fz01, fz10, fz11;
	float xslope, yslope;
	float zyx0, zyx1;
	float xold, yold;

	if (!(*image)) return;
	im = *image;
	nxmx = *nxmax;
	xscale = (*xmax - *xmin) / (im->width - 1);
	yscale = (*ymax - *ymin) / (im->height - 1);
	zscale = 1.0 / (*dz);
	for (i=0,buf=im->image; i<im->height; i++) {
/*		if (!(i%100)) {
			printf ("Interpolating row %d...\n", i);
			fflush (stdout);
		}*/
		yint = ((im->height - 1 - i) * yscale) + (*ymin);
		if (y[*ny-1] > y[0]) {
			if (i == 0) {
				for (k=0; k<*ny; k++) {
					if (yint < y[k]) {
						y1 = k;
						y0 = k-1;
						break;
					}
					y1 = k+1;
					y0 = k;
				}
			} else {
				if (yint > yold) {
					for (k=y1; k<*ny; k++) {
						if (yint < y[k]) {
							y1 = k;
							y0 = k-1;
							break;
						}
						y1 = k+1;
						y0 = k;
					}
				} else {
					for (k=y1-1; k>=0; k--) {
						if (yint > y[k]) {
							y1 = k+1;
							y0 = k;
							break;
						}
						y1 = 0;
						y0 = 0;
					}
				}
			}
		} else {
			if (i == 0) {
				for (k=0; k<*ny; k++) {
					if (yint > y[k]) {
						y1 = k;
						y0 = k-1;
						break;
					}
					y1 = k+1;
					y0 = k;
				}
			} else {
				if (yint > yold) {
					for (k=y1; k<*ny; k++) {
						if (yint > y[k]) {
							y1 = k;
							y0 = k-1;
							break;
						}
						y1 = k+1;
						y0 = k;
					}
				} else {
					for (k=y1-1; k>=0; k--) {
						if (yint < y[k]) {
							y1 = k+1;
							y0 = k;
							break;
						}
						y1 = 0;
						y0 = 0;
					}
				}
			}
		}
		yold = yint;
		if (y0 < 0) y0 = 0;
		if (y1 > *ny-1) y1 = *ny-1;
		fy0 = y[y0];
		fy1 = y[y1];
		if (y0 == y1) {
			yslope = 0.0;
		} else {
			yslope = (yint - fy0) / (fy1 - fy0);
		}
		for (j=0; j<im->width; j++,buf++) {
			xint = (j*xscale) + (*xmin);
			if (x[*nx-1] > x[0]) {
				if (j == 0) {
					for (k=0; k<*nx; k++) {
						if (xint < x[k]) {
							x1 = k;
							x0 = k-1;
							break;
						}
						x1 = k+1;
						x0 = k;
					}
				} else {
					if (xint > xold) {
						for (k=x1; k<*nx; k++) {
							if (xint < x[k]) {
								x1 = k;
								x0 = k-1;
								break;
							}
							x1 = k+1;
							x0 = k;
						}
					} else {
						for (k=x1-1; k>=0; k--) {
							if (xint > x[k]) {
								x1 = k+1;
								x0 = k;
								break;
							}
							x1 = 0;
							x0 = 0;
						}
					}
				}
			} else {
				if (j == 0) {
					for (k=0; k<*nx; k++) {
						if (xint > x[k]) {
							x1 = k;
							x0 = k-1;
							break;
						}
						x1 = k+1;
						x0 = k;
					}
				} else {
					if (xint > xold) {
						for (k=x1; k<*nx; k++) {
							if (xint > x[k]) {
								x1 = k;
								x0 = k-1;
								break;
							}
							x1 = k+1;
							x0 = k;
						}
					} else {
						for (k=x1-1; k>=0; k--) {
							if (xint < x[k]) {
								x1 = k+1;
								x0 = k;
								break;
							}
							x1 = 0;
							x0 = 0;
						}
					}
				}
			}
			xold = xint;
			if (x0 < 0) x0 = 0;
			if (x1 > *nx-1) x1 = *nx-1;
			fx0 = x[x0];
			fx1 = x[x1];
			if (x0 == x1) {
				xslope = 0.0;
			} else {
				xslope = (xint - fx0) / (fx1 - fx0);
			}
			fz00 = z[x0+nxmx*y0];
			fz01 = z[x1+nxmx*y0];
			fz10 = z[x0+nxmx*y1];
			fz11 = z[x1+nxmx*y1];
			zyx0 = fz00 + (fz10 - fz00) * yslope;
			zyx1 = fz01 + (fz11 - fz01) * yslope;
			zint = zyx0 + (zyx1 - zyx0) * xslope;
			zpix = (zint - *zmin) * zscale + 1.0;
			if (zpix < 0) zpix = 0;
			if (zpix >= im->depth) zpix = im->depth-1;
			*buf = zpix;
		}
	}
/*	printf ("done.\n");
	fflush (stdout);*/
}

/*
 *******************************************************************************
 *
 *	FORTRAN callable subroutine hdpnimg
 *
 *	Author - Danny Harvey
 *
 *******************************************************************************
 */

void
hdpnimg_ (image, x, y)

Image **image;
int *x;
int *y;

{
	Image *im;
	XImage *ximage;
	int i, n;
	unsigned char *buf;

	if (!(*image)) return;
	im = *image;

	if ( hd_depth < 8 ) {

	  int wdth, hght;
	  int xoff, yoff;
	  int colors[4*17], bpl;
	  int *data, *pixdata;
	  int fg;

	  fg = BlackPixel (hd_display, hd_screen);

	  wdth = im->width / 4;
	  hght = im->height / 4;
	  xoff = ( im->width - wdth * 4 ) / 2;
	  yoff = ( im->height - hght * 4 ) / 2;

	  switch ( hd_depth ) {
	    case 1:
	      setmncol_( colors, &fg );
	      break;
	    case 2:
	      setgrcol_( colors );
	      break;
	  }

	  /*
	   * "bpl" is the number of bytes per line.
	   * "pixdata" is the data part of the XImage structure.  Holds all
	   * the pixel information for rendering an image.
	   */
	  bpl = setbpl( 4 * wdth, hd_depth );
	  pixdata = (int *) malloc( bpl * 4 * hght );
	  data = (int *) malloc( wdth * hght * sizeof(int));

	  resamp_( data, im->image, &wdth, &hght, &(im->width), &(im->height),
		  &xoff, &yoff, &(im->depth), &hd_depth );

	  switch (hd_depth ) {
	    case 1:
	      setmono_( data, pixdata, &bpl, &wdth, &hght, colors, &hd_rv );
	      break;
	    case 2:
	      setgrey_( data, pixdata, &bpl, &wdth, &hght, colors, &hd_rv );
	      break;
	  }

	  ximage = XCreateImage( hd_display, hd_visual, hd_depth, ZPixmap, 0,
		  pixdata, 4 * wdth, 4 * hght, 32, 0 );
	  XPutImage( hd_display, hd_pixmap, hd_gc, ximage, 0, 0,
		  *x + xoff, *y + yoff, 4 * wdth, 4 * hght );
	  if (!hd_batch)
		  XPutImage( hd_display, hd_window, hd_gc, ximage, 0, 0,
			  *x + xoff, *y + yoff, 4 * wdth, 4 * hght );
	  if (!hd_batch)
		XFlush (hd_display);
	  ximage->data = NULL;
	  XDestroyImage (ximage);
	  free (data);
	  free (pixdata);
	  return;
	}

/*
 *	Paint the image onto the pixmap and screen.
 */
	n = im->width*im->height;
	buf = (unsigned char *) malloc (n);
	for (i=0; i<n; i++) buf[i] = im->ctable[im->image[i]];
	ximage = XCreateImage (hd_display, hd_visual, 8, ZPixmap,
					0, buf, im->width, im->height,
					8, 0);
	XPutImage (hd_display, hd_pixmap, hd_gc, ximage, 0, 0, *x, *y,
					im->width, im->height);
	if (!hd_batch)
		XPutImage (hd_display, hd_window, hd_gc, ximage, 0, 0, *x, *y,
					im->width, im->height);
	if (!hd_batch)
		XFlush (hd_display);
	ximage->data = NULL;
	XDestroyImage (ximage);
	free (buf);
}

/*
 *******************************************************************************
 *
 *	FORTRAN callable subroutine hdpnimgl
 *
 *	Author - Danny Harvey
 *
 *******************************************************************************
 */

void
hdpnimgl_ (image)

Image **image;

{
	Image *im;
	float x, y, w, h, xcl, ycb, fl;
	unsigned char buf;
	int i, j, n, d;
	static char line[128];

	if (!(*image)) return;
	if (hd_fd < 1) return;
	if (!hd_lon) return;
	im = *image;
	getfl_ (&xcl,&ycb,&fl);
	getclp_ (&x, &w, &y, &h);
	x += xcl;
	y += ycb;
	x = x*fl + 150.0;
	y = y*fl + 150.0;
	w += xcl;
	h += ycb;
	w = w*fl + 150.0;
	h = h*fl + 150.0;
	w = w - x + 1.0;
	h = h - y + 1.0;
/*
 *	Paint the image into the PostScript file.
 */
	if (hd_pscolor) {
		sprintf (line, "gsave\n");
		write (hd_fd, line, strlen(line));
		sprintf (line, "/picstr %d string def\n", 3*im->width);
		write (hd_fd, line, strlen(line));
		sprintf (line, "%.1f %.1f translate %f %f scale\n", x, y, w, h);
		write (hd_fd, line, strlen(line));
		sprintf (line, "%d %d 8 [%d 0 0 -%d 0 %d]\n", im->width, 
				im->height, im->width, im->height, im->height);
		write (hd_fd, line, strlen(line));
		sprintf (line,"{currentfile picstr readhexstring pop} ");
		write (hd_fd, line, strlen(line));
		sprintf (line,"false 3 colorimage\n");
		write (hd_fd, line, strlen(line));
		n = im->height*im->width;
		d = im->depth - 1;
		if (d < 1) d = 1;
		for (i=0,j=0; i<n; i++) {
			sprintf (&line[j], "%2.2x", im->red[im->image[i]]);
			j += 2;
			sprintf (&line[j], "%2.2x", im->green[im->image[i]]);
			j += 2;
			sprintf (&line[j], "%2.2x", im->blue[im->image[i]]);
			j += 2;
			if (j > 66) {
				sprintf (&line[j], "\n");
				write (hd_fd, line, strlen(line));
				j = 0;
			}
		}
		if (j > 0) {
			sprintf (&line[j], "\n");
			write (hd_fd, line, strlen(line));
		}
		sprintf (line, "grestore\n");
		write (hd_fd, line, strlen(line));
	} else {
		sprintf (line, "gsave\n");
		write (hd_fd, line, strlen(line));
		sprintf (line, "/picstr %d string def\n", im->width);
		write (hd_fd, line, strlen(line));
		sprintf (line, "%.1f %.1f translate %f %f scale\n", x, y, w, h);
		write (hd_fd, line, strlen(line));
		sprintf (line, "%d %d 8 [%d 0 0 -%d 0 %d]\n", im->width, 
				im->height, im->width, im->height, im->height);
		write (hd_fd, line, strlen(line));
		sprintf (line,"{currentfile picstr readhexstring pop} image\n");
		write (hd_fd, line, strlen(line));
		n = im->height*im->width;
		d = im->depth - 1;
		if (d < 1) d = 1;
		for (i=0,j=0; i<n; i++) {
			double intensity;
			double sqrt();
	
			if ( hd_rv == 0 ) {
				intensity = (double)(d-im->image[i])/(double)d;
			} else {
				intensity = (double) im->image[i] / (double) d;
			}
			intensity = sqrt(intensity);
			buf = 255 * intensity + 0.5;
			sprintf (&line[j], "%2.2x", buf);
			j += 2;
			if (j > 72) {
				strcat (line, "\n");
				write (hd_fd, line, strlen(line));
				j = 0;
			}
		}
		if (j > 0) {
			strcat (line, "\n");
			write (hd_fd, line, strlen(line));
		}
		sprintf (line, "grestore\n");
		write (hd_fd, line, strlen(line));
	}
}

cset_( mychar, myint )
char *mychar;
int *myint;
{
	*mychar = *myint;
}

cadd_( mychar, myint )
char *mychar;
int *myint;
{
	*mychar = *mychar + *myint;
}

void
cewrite_( string, len )
char	*string;
int	len;
{
	fprintf( stderr, "%s\n", string );
}

int
setbpl( im_w, depth )
/*
 *	setbpl()
 *
 *	Calculate the number of bytes needed to fully give pixel values on one
 *	image scanline.
 */
	int	im_w;		/* width of the image in pixels */
	int	depth;		/* screen depth */
{
	int	bpl;		/* result of calculation = bytes_per_line */

	switch ( depth ) {

	  case 1:
	    bpl = im_w / 32;
	    if ( ( bpl * 32 ) != im_w ) bpl = bpl + 1;
	    bpl = bpl * 4;
	    break;

	  case 2:
	    bpl = im_w / 16;
	    if ( ( bpl * 16 ) != im_w ) bpl = bpl + 1;
	    bpl = bpl * 4;
	    break;

	}

	return( bpl );
}

void
hdreverse_video(rv)
int rv;
{
	hd_rv = rv;
}
