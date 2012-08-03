/*
 ******************************************************************************
 *
 *	Function hdbatoff_
 *
 *	Author - Danny Harvey
 *
 ******************************************************************************
 */

int hd_batch;

void
hdbatoff_()

/*
 *	hdbatoff_ will disable plot batching which causes plotting to occur
 *	in the pixmap and also in the current window.
 */

{
	hd_batch = 0;
}
