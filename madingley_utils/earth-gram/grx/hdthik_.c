/*
 ******************************************************************************
 *
 *	Function hdthik_
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

int hd_fd;
int hd_lon;

int hd_lastth = -1;

void
hdthik_(n)

long *n;

/*
 *	hdthik_ will set the line thickness
 */

{
	char line[80];
	int i;

	if (hd_fd < 1) return;
	if (!hd_lon) return;
	if (*n < 1) {
		i = 0;
	} else {
		i = *n;
	}
	if (i == hd_lastth) return;
	if (i < 1) {
		write (hd_fd, "0 setlinewidth\n", 15);
	} else {
		sprintf (line, "%d setlinewidth\n", i);
		write (hd_fd, line, strlen(line));
	}
	hd_lastth = i;
}
