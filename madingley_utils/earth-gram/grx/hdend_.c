/*
 ******************************************************************************
 *
 *	Function hdend_
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

long hd_fd;
long hd_ldraw;

void
hdend_()

/*
 *	hdend_ will cause the plot to be stroked and the page to be displayed
 */

{
	void hdendl();

	if (hd_fd > 0) hdendl();
}

void
hdendl()

{
/*	write (hd_fd, "%%EOF\n", 6);*/
	write (hd_fd, "%%EndDocument\n", 14);
	write (hd_fd, "showpage\n", 9);
	write (hd_fd, "%%EOF\n", 6);
	close (hd_fd);
	hd_fd = 0;
	hd_ldraw = 0;
}
