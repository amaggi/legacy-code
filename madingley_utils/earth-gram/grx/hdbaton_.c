/*
 ******************************************************************************
 *
 *	Function hdbaton_
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

int hd_batch;

void
hdbaton_()

/*
 *	hdbaton_ will enable plot batching which causes plotting to ONLY occur
 *	in the pixmap and not in the current window.
 */

{
	hd_batch = 1;
}
