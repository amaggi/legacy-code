/*
 ******************************************************************************
 *
 *	Function hdfore_
 *
 *	Author - Danny Harvey
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
#ifdef STELLAR
XFDIGC hd_gc3;
#endif

int hd_fd;
int hd_lon;
int hd_pscolor;

float hd_fore_light;
float hd_fore_sat;
float hd_fore_red;
float hd_fore_green;
float hd_fore_blue;

void
hdfore_(hue, light, sat)

float *hue;
float *light;
float *sat;

/*
 *	hdfore_ will set the foreground color.
 */

{
	static char line[128];

	unsigned long int hdhls_to_pixel();

	hd_fore = hdhls_to_pixel (*hue, *light, *sat);
	hdhls_to_rgb(*hue,*light,*sat, 
		&hd_fore_red, &hd_fore_green, &hd_fore_blue);
	hd_fore_light = *light;
	hd_fore_sat = *sat;
	XSetForeground (hd_display, hd_gc, hd_fore);
#ifdef STELLAR_OFF
	XFDISetLineColor (hd_display, hd_gc3, hd_fore);
#endif
	if (hd_fd < 1) return;
	if (!hd_lon) return;
	if (hd_pscolor == 0 || hd_pscolor == 2) return;
	if (hd_pscolor == 3) {
		if (hd_fore_light > 0.5) hd_fore_light = 0.0;
		if ((*hue) > 30.0 && (*hue) < 90.0) {
			hd_fore_light = 0.0;
		}
		hdhls_to_rgb(*hue,hd_fore_light,*sat, 
			&hd_fore_red, &hd_fore_green, &hd_fore_blue);
	}
	if (hd_fore_light == 0.0 || hd_fore_light == 1.0 || hd_fore_sat==0.0) {
		sprintf (line, "%.3f setgray\n", hd_fore_light);
	} else {
		sprintf (line, "%.3f %.3f %.3f setrgbcolor\n", 
				hd_fore_red, hd_fore_green, hd_fore_blue);
	}
	write (hd_fd, line, strlen(line));
}

/*
 ******************************************************************************
 *
 *	hdhls_to_rgb()
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

int
hdhls_to_rgb(h,l,s,r,g,b)

float h,l,s;
float *r,*g,*b;

/*
 *	hls_to_rgb will do a hue-lightness-saturation to red-green-blue
 *	transformation.
 *
 *	Inputs  -	h	= A floating hue value. This should be in
 *				  the range of 0. to 360. with both 0. and
 *				  360. being red. This describes an angle
 *				  around the color wheel. 60. is yellow,
 *				  120. is green, 180. is cyan, 240. is
 *				  blue and 300. is magenta.
 *			l	= A floating lightness value. This should be in
 *				  the range of 0. to 1. with 0. being black
 *				  and 1. being white.
 *			s	= A floating saturation value. This should be in
 *				  the range of 0. to 1. with 0. being a gray
 *				  shade and 1. being a pure color.
 *
 *	Outputs -	r	= A pointer to a floating red value. This should
 *				  be in the range of 0. to 1. with 0. being
 *				  no red and 1. being maximum red.
 *			g	= A pointer to a floating green value. This 
 *				  should be in the range of 0. to 1. with 
 *				  0. being no green and 1. being maximum green.
 *			b	= A pointer to a floating blue value. This 
 *				  should be in the range of 0. to 1. with 
 *				  0. being no blue and 1. being maximum blue.
 */

{
	float m1,m2;

	float hdvalue();

	if (l < 0.0) l = 0.0;
	if (l > 1.0) l = 1.0;
	if (s < 0.0) s = 0.0;
	if (s > 1.0) s = 1.0;
	if (s == 0.0) {
		*r = l;
		*g = l;
		*b = l;
	} else {
		if (l <= 0.5) {
			m2 = l*(1+s);
		} else {
			m2 = l + s - l*s;
		}
		m1 = 2*l - m2;
		*r = hdvalue(m1,m2,h+120.0);
		*g = hdvalue(m1,m2,h);
		*b = hdvalue(m1,m2,h-120.0);
	}
	return (0);
}

float
hdvalue(m1,m2,h)

float m1,m2,h;

{
	while (h > 360.0) h -= 360.0;
	while (h < 0.0) h += 360.0;
	if (h < 60.0) {
		return (m1+(m2-m1)*h/60.0);
	} else if (h < 180.0) {
		return (m2);
	} else if (h < 240.0) {
		return (m1+(m2-m1)*(240.0-h)/60.0);
	} else {
		return (m1);
	}
}

/*
 ******************************************************************************
 *
 *	hdrgb_to_hls()
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

int
hdrgb_to_hls(r,g,b,h,l,s)

float r,g,b;
float *h,*l,*s;

/*
 *	rgb_to_hls will do a red-green-blue to hue-lightness-saturation
 *	transformation.
 *
 *	Inputs -	r	= A floating red value. This should
 *				  be in the range of 0. to 1. with 0. being
 *				  no red and 1. being maximum red.
 *			g	= A floating green value. This 
 *				  should be in the range of 0. to 1. with 
 *				  0. being no green and 1. being maximum green.
 *			b	= A floating blue value. This 
 *				  should be in the range of 0. to 1. with 
 *				  0. being no blue and 1. being maximum blue.
 *
 *	Outputs  -	h	= A pointer to a floating hue value. This should
 *				  be in the range of 0. to 360. with both 0. and
 *				  360. being red. This describes an angle
 *				  around the color wheel. 60. is yellow,
 *				  120. is green, 180. is cyan, 240. is
 *				  blue and 300. is magenta.
 *			l	= A pointer to a floating lightness value. This 
 *				  should be in the range of 0. to 1. with 0. 
 *				  being black and 1. being white.
 *			s	= A pointer to a floating saturation value. This
 *				  should be in the range of 0. to 1. with 0. 
 *				  being a gray shade and 1. being a pure color.
 */

{
	float max, min, rc, gc, bc;

	if (r < 0.0) r = 0.0;
	if (r > 1.0) r = 1.0;
	if (g < 0.0) g = 0.0;
	if (g > 1.0) g = 1.0;
	if (b < 0.0) b = 0.0;
	if (b > 1.0) b = 1.0;
	min = r;
	if (g < min) min = g;
	if (b < min) min = b;
	max = r;
	if (g > max) max = g;
	if (b > max) max = b;
	*l = (max+min) * 0.5;
	if (max == min) {
		*s = 0.0;
		*h = 0.0;
	} else {
		if (*l <= 0.5) *s = (max - min)/(max + min);
			else   *s = (max - min)/(2.0 - max - min);
		rc = (max - r) / (max - min);
		gc = (max - g) / (max - min);
		bc = (max - b) / (max - min);
		if (r == max) {
			*h = bc - gc;
		} else if (g == max) {
			*h = 2.0 + rc - bc;
		} else {
			*h = 4.0 + gc - rc;
		}
		*h *= 60.0;
		if (*h < 0.0) *h += 360.0;
	}
}

unsigned long int
hdhls_to_pixel (h, l, s)

float h;
float l;
float s;

{
	float r, g, b;
	int ir, ig, ib;
	XColor colorcell;

	hdhls_to_rgb (h, l, s, &r, &g, &b);
	if (hd_visual->class == TrueColor) {
		colorcell.red = r*255.9999;
		colorcell.green = g*255.9999;
		colorcell.blue = b*255.9999;
		colorcell.pixel = ((colorcell.red<<16)&0xff0000)
			|((colorcell.green<<8)&0xff00)|(colorcell.blue&0xff);
	} else {
		ir = r*65535.9;
		ig = g*65535.9;
		ib = b*65535.9;
		if (ir < 0) ir = 0;
		if (ig < 0) ig = 0;
		if (ib < 0) ib = 0;
		if (ir > 65535) ir = 65535;
		if (ig > 65535) ig = 65535;
		if (ib > 65535) ib = 65535;
		colorcell.red = ir;
		colorcell.green = ig;
		colorcell.blue = ib;
		XAllocColor (hd_display, hd_cmap, &colorcell);
	}
	return (colorcell.pixel);
}

unsigned long int
hdhls_to_pixelrgb (h, l, s, ir, ig, ib)

float h;
float l;
float s;
int *ir;
int *ig;
int *ib;

{
	float r, g, b;
	int jr, jg, jb;
	XColor colorcell;

	hdhls_to_rgb (h, l, s, &r, &g, &b);
	if (hd_visual->class == TrueColor) {
		colorcell.red = r*255.9999;
		colorcell.green = g*255.9999;
		colorcell.blue = b*255.9999;
		colorcell.pixel = ((colorcell.red<<16)&0xff0000)
			|((colorcell.green<<8)&0xff00)|(colorcell.blue&0xff);
	} else {
		jr = r*65535.9;
		jg = g*65535.9;
		jb = b*65535.9;
		if (jr < 0) jr = 0;
		if (jg < 0) jg = 0;
		if (jb < 0) jb = 0;
		if (jr > 65535) jr = 65535;
		if (jg > 65535) jg = 65535;
		if (jb > 65535) jb = 65535;
		colorcell.red = jr;
		colorcell.green = jg;
		colorcell.blue = jb;
		XAllocColor (hd_display, hd_cmap, &colorcell);
	}
	*ir = colorcell.red;
	*ig = colorcell.green;
	*ib = colorcell.blue;
	return (colorcell.pixel);
}
