/*
 *******************************************************************************
 *
 *	FORTRAN callable subroutine hdptxwd
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

void
hdptxwd_ (xwd_file, imgx, imgy, n)

char *xwd_file;
int *imgx;
int *imgy;
int n;

/*
 *	hdptxwd will put the window image from the output of xwd into the
 *	pixmap.
 *
 *	Inputs  -	xwd_file	= xwd output file name.
 *			imgx		= Starting X-pixel in image grid.
 *			imgy		= Starting Y-pixel in image grid.
 */

{
	static XImage *image=NULL;
	int width;
	int height;
	int depth;
	static char file[512] = "";
	static char o_file[512] = "";

	XImage *read_xwd();
	char *convert_string();

	strcpy (file, convert_string(xwd_file, n));
	if (strcmp(file, o_file)) {
		if (image) {
			if (image->data) free (image->data);
			image->data = NULL;
			XDestroyImage (image);
		}
		image = read_xwd (hd_display, hd_visual, hd_depth, file, 
							&width, &height);
		if (image == NULL) return;
		strcpy (o_file, file);
	}
	XPutImage (hd_display, hd_pixmap, hd_gc, image, *imgx, *imgy, 0, 0, 
								width, height);
	if (!hd_batch) {
		XPutImage (hd_display, hd_window, hd_gc, image, *imgx,*imgy,0,0,
								width, height);
		XFlush (hd_display);
	}
}

#include <fcntl.h>
#include "XWDFile.h"

XImage *
read_xwd (display, visual, depth, file_name, width, height)

Display * display;
Visual *           visual;
int                        depth;
char *                            file_name;
int *                                        width;
int *                                               height;

{
	int fd;
	int i, j, k, n;
	static XWDFileHeader head;
	static XWDColor color;
	unsigned char *bytes;
	XImage *image;

	fd = open (file_name, O_RDONLY);
	if (fd < 0) {
		fprintf (stderr, "read_xwd: Unable to open file %s.\n",
							file_name);
		return (NULL);
	}
	n = read (fd, &head, sizeof(XWDFileHeader));
	if (n != sizeof(XWDFileHeader)) {
		fprintf (stderr, "read_xwd: Read error on file %s.\n",
							file_name);
		close (fd);
		return (NULL);
	}
	lseek (fd, head.header_size, 0);
	for (i=0; i<head.colormap_entries; i++) {
		n = read (fd, &color, sizeof(XWDColor));
		if (n != sizeof(XWDColor)) {
			fprintf (stderr, "read_xwd: Read error on file %s.\n",
								file_name);
			close (fd);
			return (NULL);
		}
	}
	bytes = (unsigned char *) 
		malloc (head.bytes_per_line*head.pixmap_height);
	if (bytes == NULL) {
		fprintf (stderr, "read_xwd: Malloc error.\n");
		close (fd);
		return (NULL);
	}
	for (i=0,k=0; i<head.pixmap_height; i++) {
		n = read (fd, &bytes[k], head.bytes_per_line);
		if (n != head.bytes_per_line) {
			fprintf (stderr, "read_xwd: Read error on file %s.\n",
							file_name);
			free (bytes);
			close (fd);
			return (NULL);
		}
		k += head.bytes_per_line;
	}
	close (fd);
	if (head.pixmap_depth == depth) {
		image = XCreateImage (display, visual, depth, 
						head.pixmap_format,
						0, bytes, 
						head.pixmap_width,
						head.pixmap_height,
						head.bitmap_pad,
						head.bytes_per_line);
	} else {
		if (head.pixmap_depth == 1) {
			image = XCreateImage (display, visual, 1,
						XYBitmap,
						0, bytes, 
						head.pixmap_width,
						head.pixmap_height,
						head.bitmap_pad,
						head.bytes_per_line);
		} else {
			image = NULL;
		}
	}
	if (image == NULL) {
		fprintf (stderr, "read_xwd: XCreateImage() error.\n");
		return (NULL);
	}
	*width = head.pixmap_width;
	*height = head.pixmap_height;
	return (image);
}
