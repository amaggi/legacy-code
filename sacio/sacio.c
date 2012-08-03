/*
 *	SAC IO functions for ANSI C
 *
 * $Id: sacio.c,v 1.2 2005/02/09 01:48:21 alessia Exp $
 *
 *
 *	To use the C functions:
 *  ----------------------------------
 *
 *		you need the files: sacio.c and sacio.h
 *
 *		If you intend to manipulate the headers, you must include sacio.h
 *		in the c-code that does the manipulating since the sac file header
 *		is stored in a structure defined in sac.h. It's much more convenient
 *		than the old getfhv()... in FORTRAN.
 *
 *		If you are setting a lot of header values, use ReadSACfile
 *		and WriteSACfile.
 *
 *		If you are doing bare bones headers, like wsac1 and wsac2 in FORTRAN,
 *			just use: WriteSAC1, WriteSAC2
 *
 *	To use the FORTRAN stubs:
 *  ----------------------------------
 *
 *              The header access functions are incomplete, although most 
 *                common header values are included in the [g,s]et[f,n,l]hv
 *                functions. I did not write a [g,s]etkhv routine at all. These  
 *                are what I needed right now to get some of my codes running on 
 *                machines where sac is not available.
 *
 *                To use these functions, make the library and then link your
 *                 FORTRAN code with this lib in place of the sac library.
 *
 *
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "sacio.h"


/*****************************************************************************/

void NewSacHeader(struct SacHeader *hd)
{
	*hd = sac_null;
	hd->nvhdr = 6;
	hd->ninf = 0;
	hd->nhst = 0;
	hd->lpspol = FALSE;
	hd->lcalda = TRUE;
	hd->unused27 = FALSE;
}

/*****************************************************************************/
/*
 * FUNCTION TO READ IN A SAC HEADER
 * 
 * read in a sac file header filename is the name of the SAC file header is the
 * SAC file header io_error = 0 upon success, -1 on failure
 */
/*****************************************************************************/

void ReadSACHeader(char *filename, struct SacHeader *header, int *io_error)
{
	FILE           *strm;
	int             nread;

	*io_error = kNO_ERROR;

	/* open the stream */
	strm = fopen(filename, "rb");
	if (strm == NULL) {
		*io_error = kIO_ERROR;
		return;
	}
	
	/* read in the header values */
	nread = fread(header, sizeof(struct SacHeader), 1, strm);
	
	fclose(strm);

	if (nread != 1)
		*io_error = kIO_ERROR;
}

/*****************************************************************************/
/*
 * FUNCTION TO READ IN AN EVENLY_SPACED SAC FILE
 * 
 *	filename is the name of the SAC file header
 *	x is the array for the time samples
 *	header is the sac header structure
 *	maxpts is the max number of points to input
 *	io_error = 0 upon success, -1 on failure
 * 
 */
/*****************************************************************************/

void ReadSACfile(char *filename,float *x,struct SacHeader *header,int maxpts, int *io_error)
{
	FILE           *strm;
	int             nread;

	*io_error = kNO_ERROR;

	/* open the stream */
	strm = fopen(filename, "rb");
	if (strm == NULL) {
		*io_error = kIO_ERROR;
		return;
	}
	/* read in the file header */

	nread = fread(header, sizeof(struct SacHeader), 1, strm);
	if (nread != 1) {
		*io_error = kIO_ERROR;
		fclose(strm);
		return;
	}
	
	/* check to make sure we don't overfill memory */
	if (header->npts > maxpts)
	{
		header->npts = maxpts;
	}
	
	/* read in the floating point amplitudes */

	nread = fread(x, header->npts * sizeof(float), 1, strm);

	fclose(strm);
	if (nread != 1){
		*io_error = kIO_ERROR;
        }
}

/*****************************************************************************/
/*
 * FUNCTION TO READ IN AN EVENLY_SPACED ASCII SAC FILE
 * 
 *	filename is the name of the ascii SAC file 
 *	x is the array for the time samples
 *	header is the sac header structure
 *	maxpts is the max number of points to input
 *	io_error = 0 upon success, -1 on failure
 * 
 */
/*****************************************************************************/

void ReadAsciiSACfile(char *filename,float *x,struct SacHeader *header,int maxpts, int *io_error)
{
	FILE           *strm;
	int             nfives, nleft, i;

	*io_error = kNO_ERROR;

	/* open the stream */
	strm = fopen(filename, "r"); /* ascii file */
	if (strm == NULL) {
		*io_error = kIO_ERROR;
		return;
	}

	/* read in the file header */
        fscanf(strm,"%g %g %g %g %g", &header->delta, &header->depmin, &header->depmax, &header->scale, &header->odelta);
        fscanf(strm,"%g %g %g %g %g", &header->b, &header->e, &header->o, &header->a, &header->internal1);
        fscanf(strm,"%g %g %g %g %g", &header->t0, &header->t1, &header->t2, &header->t3, &header->t4);
        fscanf(strm,"%g %g %g %g %g", &header->t5, &header->t6, &header->t7, &header->t8, &header->t9);
        fscanf(strm,"%g %g %g %g %g", &header->f, &header->resp0, &header->resp1, &header->resp2, &header->resp3);
        fscanf(strm,"%g %g %g %g %g", &header->resp4, &header->resp5, &header->resp6, &header->resp7, &header->resp8);
        fscanf(strm,"%g %g %g %g %g", &header->resp9, &header->stla, &header->stlo, &header->stel, &header->stdp);
        fscanf(strm,"%g %g %g %g %g", &header->evla, &header->evlo, &header->evel, &header->evdp, &header->unused1);
        fscanf(strm,"%g %g %g %g %g", &header->user0, &header->user1, &header->user2, &header->user3, &header->user4);
        fscanf(strm,"%g %g %g %g %g", &header->user5, &header->user6, &header->user7, &header->user8, &header->user9);
        fscanf(strm,"%g %g %g %g %g", &header->dist, &header->az, &header->baz, &header->gcarc, &header->sb);
        fscanf(strm,"%g %g %g %g %g", &header->sdelta, &header->depmen, &header->cmpaz, &header->cmpinc, &header->xminimum);
        fscanf(strm,"%g %g %g %g %g", &header->xmaximum, &header->yminimum, &header->ymaximum, &header->unused6, &header->unused7);
        fscanf(strm,"%g %g %g %g %g", &header->unused8, &header->unused9, &header->unused10, &header->unused11, &header->unused12);
        fscanf(strm,"%ld %ld %ld %ld %ld", &header->nzyear, &header->nzjday, &header->nzhour, &header->nzmin, &header->nzsec);
        fscanf(strm,"%ld %ld %ld %ld %ld", &header->nzmsec, &header->nvhdr, &header->ninf, &header->nhst, &header->npts);
        fscanf(strm,"%ld %ld %ld %ld %ld", &header->nspts, &header->nsn, &header->nxsize, &header->nysize, &header->unused15);
        fscanf(strm,"%ld %ld %ld %ld %ld", &header->iftype, &header->idep, &header->iztype, &header->unused16, &header->iinst);
        fscanf(strm,"%ld %ld %ld %ld %ld", &header->istreg, &header->ievreg, &header->ievtyp, &header->iqual, &header->isynth);
        fscanf(strm,"%ld %ld %ld %ld %ld", &header->unused17, &header->unused18, &header->unused19, &header->unused20, &header->unused21);
        fscanf(strm,"%ld %ld %ld %ld %ld", &header->unused22, &header->unused23, &header->unused24, &header->unused25, &header->unused26);
        fscanf(strm,"%ld %ld %ld %ld %ld\n", &header->leven, &header->lpspol, &header->lovrok, &header->lcalda, &header->unused27);
        fscanf(strm,"%8c %16c\n",header->kstnm, header->kevnm);
        fscanf(strm,"%8c%8c%8c\n",header->khole, header->ko, header->ka);
        fscanf(strm,"%8c%8c%8c\n",header->kt0, header->kt1, header->kt2);
        fscanf(strm,"%8c%8c%8c\n",header->kt3, header->kt4, header->kt5);
        fscanf(strm,"%8c%8c%8c\n",header->kt6, header->kt7, header->kt8);
        fscanf(strm,"%8c%8c%8c\n",header->kt9, header->kf, header->kuser0);
        fscanf(strm,"%8c%8c%8c\n",header->kuser1, header->kuser2, header->kcmpnm);
        fscanf(strm,"%8c%8c%8c\n",header->knetwk, header->kdatrd, header->kinst);
	
	/* check to make sure we don't overfill memory */
	if (header->npts > maxpts)
	{
		header->npts = maxpts;
	}


        /* find how many complete rows we have by integer division */
        nfives = header->npts / 5;
        nleft = header->npts % 5;

	/* read in the floating point amplitudes */
        for(i=0;i<nfives;i++){
            fscanf(strm,"%g %g %g %g %g", &x[i*5], &x[i*5+1], &x[i*5+2], &x[i*5+3], &x[i*5+4]);
        }

        switch(nleft){
          case 0:
            break;
          case 1:
            fscanf(strm,"%g", &x[(i+1)*5]);
            break;
          case 2:
            fscanf(strm,"%g %g", &x[(i+1)*5], &x[(i+1)*5+1]);
            break;
          case 3:
            fscanf(strm,"%g %g %g", &x[(i+1)*5], &x[(i+1)*5+1], &x[(i+1)*5+2]);
            break;
          case 4:
            fscanf(strm,"%g %g %g %g", &x[(i+1)*5], &x[(i+1)*5+1], &x[(i+1)*5+2], &x[(i+1)*5+3]);
            break;
        }

	fclose(strm);

        return;
}
/*****************************************************************************/
/*
 * FUNCTION TO READ IN AN EVENLY_SPACED SAC FILE
 * 
 * Assuming that you don't want header details
 *
 * 	filename is the name of the SAC file 
 * 	x is the array into which the data are read
 *	npts = number of points read in (<= maxpts)
 *	dt is the sample rate
 *	b is the time of the first sample
 *	io_error = 0 upon success, -1 on failure
 */
/*****************************************************************************/

void ReadSAC1(char *filename,float *x,int *npts,float *dt,float *b,int maxpts,int *io_error)
{
	struct SacHeader	hd;

	ReadSACfile(filename,x,&hd,maxpts,io_error);
	
	if(*io_error == kIO_ERROR)
	{
		return;
	}	

	*npts = hd.npts;
	*dt = hd.delta;
	*b =  hd.b;

}
/*****************************************************************************/
/*
 * FUNCTION TO READ AN UNEVENLY_SPACED SAC FILE
 * 
 *	filename is the name of the SAC file
 *	n is the number of data points read in
 *	x and y are the data read in
 *	maxpts is the maximum number of data to read in (dimensions of x and y)
 *	io_error = 0 upon success, -1 on failure
 * 
 */
/*****************************************************************************/
void ReadSAC2(char *filename, int *n, float *x, float *y, int maxpts, int *io_error)
{
	FILE *strm;
	struct SacHeader hd;
	int nread;
	
	/* OPEN THE stream for reading the file */

	*io_error = kNO_ERROR;

	/* open the stream */
	strm = fopen(filename, "r");
	if (strm == NULL) {
		*io_error = kIO_ERROR;
		return;
	}	
	
	/* Now read the sac header to the file */
	nread = fread(&hd, sizeof(struct SacHeader), 1, strm);
		if (nread != 1) {
		*io_error = kIO_ERROR;
		fclose(strm);
		return;
	}
	
	/* set up the number of points in the file, and apply the limits */
	*n = hd.npts;
	if(*n > maxpts) *n = maxpts;
	
	/* Now read the sac data points x, then y, from the file */
	nread = fread(x, sizeof(float), *n, strm);
		if (nread != 1) {
		*io_error = kIO_ERROR;
		fclose(strm);
		return;
	}
	nread = fread(y, sizeof(float), *n, strm);
		if (nread != 1) {
		*io_error = kIO_ERROR;
		fclose(strm);
		return;
	}
	
	fclose(strm);
}

/*****************************************************************************/
/*
 * FUNCTION TO READ AN UNEVENLY_SPACED 3D SAC FILE
 * 
 *	filename is the name of the SAC file
 *	n is the number of data points read in
 *      nx, ny are the number of points in the x and y directions
 *      xmin, xmax, ymin, ymax are the x and y extrema
 *	z are the data read in
 *	maxpts is the maximum number of data to read in (max dimension of z)
 *	io_error = 0 upon success, -1 on failure
 * 
 */
/*****************************************************************************/
void ReadSAC3(char *filename, int *n, int *nx, int *ny, float *xmin, float *xmax, float *ymin, float *ymax, float *z, int maxpts, int *io_error)
{
	FILE *strm;
	struct SacHeader hd;
	int nread;
	
	/* OPEN THE stream for reading the file */

	*io_error = kNO_ERROR;

	/* open the stream */
	strm = fopen(filename, "r");
	if (strm == NULL) {
		*io_error = kIO_ERROR;
		return;
	}	
	
	/* Now read the sac header to the file */
	nread = fread(&hd, sizeof(struct SacHeader), 1, strm);
		if (nread != 1) {
		*io_error = kIO_ERROR;
		fclose(strm);
		return;
	}
	
	/* set up the number of points in the file, and apply the limits */

	*n = hd.npts;
	*nx = hd.nxsize;
	*ny = hd.nysize;
	*xmin = hd.xminimum;
	*xmax = hd.xmaximum;
	*ymin = hd.yminimum;
	*ymax = hd.ymaximum;

	if(*n > maxpts) *n = maxpts;
	
	/* Now read the sac data points z, from the file */
	nread = fread(z, sizeof(float), *n, strm);
		if (nread != 1) {
		*io_error = kIO_ERROR;
		fclose(strm);
		return;
	}
	
	fclose(strm);
}

/*****************************************************************************/



/*****************************************************************************/
/*
 * FUNCTION TO WRITE AN EVENLY_SPACED SAC FILE
 * 
 * Assuming that the header is all set up
 * 
 * filename is the name of the SAC file header is the file header ar is the
 * array into which the data are read io_error = 0 upon success, -1 on
 * failure
 */
/*****************************************************************************/

void WriteSACfile(char *filename,float *x,float *y,struct SacHeader header,int *io_error) {
	FILE           *strm;
	int             nwritten, n;

	*io_error = kNO_ERROR;

	/* open the stream */
	strm = fopen(filename, "w");
	if (strm == NULL) {
		*io_error = kIO_ERROR;
		return;
	}
	/* write in the file header */

	nwritten = fwrite(&header, sizeof(struct SacHeader), 1, strm);
	if (nwritten != 1) {
		*io_error = kIO_ERROR;
		fclose(strm);
		return;
	}
        
        if((n=header.npts) == -12345){
          fprintf(stderr,"Npts not set.  Aborting writing sac file\n");
        }
	/* Now write the sac data points y, then x, to the file */
	nwritten = fwrite(y, sizeof(float), n, strm);
		if (nwritten != n) {
		*io_error = kIO_ERROR;
		fclose(strm);
		return;
	}
	/* only write the x out if leven is false */
        if(header.leven == FALSE){
	  nwritten = fwrite(x, sizeof(float), n, strm);
		if (nwritten != n) {
		*io_error = kIO_ERROR;
		fclose(strm);
		return;
  	  }
        }

	fclose(strm);

}

/*****************************************************************************/
/*
 * FUNCTION TO WRITE AN EVENLY_SPACED ASCII SAC FILE
 * 
 * Assuming that the header is all set up
 * 
 * filename is the name of the SAC file header is the file header ar is the
 * array into which the data are read io_error = 0 upon success, -1 on
 * failure
 */
/*****************************************************************************/

void WriteAsciiSACfile(char *filename,float *x,struct SacHeader header,int *io_error) {
	FILE           *strm;
	int             nfives, nleft, i;

	*io_error = kNO_ERROR;

	/* open the stream */
	strm = fopen(filename, "w");
	if (strm == NULL) {
		*io_error = kIO_ERROR;
		return;
	}

	/* read in the file header */
        fprintf(strm,"%15.7g%15.7g%15.7g%15.7g%15.7g\n", header.delta, header.depmin, header.depmax, header.scale, header.odelta);
        fprintf(strm,"%15.7g%15.7g%15.7g%15.7g%15.7g\n", header.b, header.e, header.o, header.a, header.internal1);
        fprintf(strm,"%15.7g%15.7g%15.7g%15.7g%15.7g\n", header.t0, header.t1, header.t2, header.t3, header.t4);
        fprintf(strm,"%15.7g%15.7g%15.7g%15.7g%15.7g\n", header.t5, header.t6, header.t7, header.t8, header.t9);
        fprintf(strm,"%15.7g%15.7g%15.7g%15.7g%15.7g\n", header.f, header.resp0, header.resp1, header.resp2, header.resp3);
        fprintf(strm,"%15.7g%15.7g%15.7g%15.7g%15.7g\n", header.resp4, header.resp5, header.resp6, header.resp7, header.resp8);
        fprintf(strm,"%15.7g%15.7g%15.7g%15.7g%15.7g\n", header.resp9, header.stla, header.stlo, header.stel, header.stdp);
        fprintf(strm,"%15.7g%15.7g%15.7g%15.7g%15.7g\n", header.evla, header.evlo, header.evel, header.evdp, header.unused1);
        fprintf(strm,"%15.7g%15.7g%15.7g%15.7g%15.7g\n", header.user0, header.user1, header.user2, header.user3, header.user4);
        fprintf(strm,"%15.7g%15.7g%15.7g%15.7g%15.7g\n", header.user5, header.user6, header.user7, header.user8, header.user9);
        fprintf(strm,"%15.7g%15.7g%15.7g%15.7g%15.7g\n", header.dist, header.az, header.baz, header.gcarc, header.sb);
        fprintf(strm,"%15.7g%15.7g%15.7g%15.7g%15.7g\n", header.sdelta, header.depmen, header.cmpaz, header.cmpinc, header.xminimum);
        fprintf(strm,"%15.7g%15.7g%15.7g%15.7g%15.7g\n", header.xmaximum, header.yminimum, header.ymaximum, header.unused6, header.unused7);
        fprintf(strm,"%15.7g%15.7g%15.7g%15.7g%15.7g\n", header.unused8, header.unused9, header.unused10, header.unused11, header.unused12);
        fprintf(strm,"%10ld%10ld%10ld%10ld%10ld\n", header.nzyear, header.nzjday, header.nzhour, header.nzmin, header.nzsec);
        fprintf(strm,"%10ld%10ld%10ld%10ld%10ld\n", header.nzmsec, header.nvhdr, header.ninf, header.nhst, header.npts);
        fprintf(strm,"%10ld%10ld%10ld%10ld%10ld\n", header.nspts, header.nsn, header.nxsize, header.nysize, header.unused15);
        fprintf(strm,"%10ld%10ld%10ld%10ld%10ld\n", header.iftype, header.idep, header.iztype, header.unused16, header.iinst);
        fprintf(strm,"%10ld%10ld%10ld%10ld%10ld\n", header.istreg, header.ievreg, header.ievtyp, header.iqual, header.isynth);
        fprintf(strm,"%10ld%10ld%10ld%10ld%10ld\n", header.unused17, header.unused18, header.unused19, header.unused20, header.unused21);
        fprintf(strm,"%10ld%10ld%10ld%10ld%10ld\n", header.unused22, header.unused23, header.unused24, header.unused25, header.unused26);
        fprintf(strm,"%10ld%10ld%10ld%10ld%10ld\n", header.leven, header.lpspol, header.lovrok, header.lcalda, header.unused27);
        fprintf(strm,"%8.8s %16.16s\n",header.kstnm, header.kevnm);
        fprintf(strm,"%8.8s%8.8s%8.8s\n",header.khole, header.ko, header.ka);
        fprintf(strm,"%8.8s%8.8s%8.8s\n",header.kt0, header.kt1, header.kt2);
        fprintf(strm,"%8.8s%8.8s%8.8s\n",header.kt3, header.kt4, header.kt5);
        fprintf(strm,"%8.8s%8.8s%8.8s\n",header.kt6, header.kt7, header.kt8);
        fprintf(strm,"%8.8s%8.8s%8.8s\n",header.kt9, header.kf, header.kuser0);
        fprintf(strm,"%8.8s%8.8s%8.8s\n",header.kuser1, header.kuser2, header.kcmpnm);
        fprintf(strm,"%8.8s%8.8s%8.8s\n",header.knetwk, header.kdatrd, header.kinst);
	
        /* find how many complete rows we have by integer division */
        nfives = header.npts / 5;
        nleft = header.npts % 5;

	/* read in the floating point amplitudes */
        for(i=0;i<nfives;i++){
            fprintf(strm,"%15.7g%15.7g%15.7g%15.7g%15.7g\n", x[i*5], x[i*5+1], x[i*5+2], x[i*5+3], x[i*5+4]);
        }

        switch(nleft){
          case 0:
            break;
          case 1:
            fprintf(strm,"%15.7g\n", x[(i+1)*5]);
            break;
          case 2:
            fprintf(strm,"%15.7g%15.7g\n", x[(i+1)*5], x[(i+1)*5+1]);
            break;
          case 3:
            fprintf(strm,"%15.7g%15.7g%15.7g\n", x[(i+1)*5], x[(i+1)*5+1], x[(i+1)*5+2]);
            break;
          case 4:
            fprintf(strm,"%15.7g%15.7g%15.7g%15.7g\n", x[(i+1)*5], x[(i+1)*5+1], x[(i+1)*5+2], x[(i+1)*5+3]);
            break;
        }

	fclose(strm);

 
        return;

}


/*****************************************************************************/
/*
 * FUNCTION TO WRITE AN EVENLY_SPACED SAC FILE
 * 
 * Assuming that the header is NOT set up ---
 * 
 * filename is the name of the SAC file, x is the
 * array from which the data are written,npts is the number of pts,
 * dt is the sample increment, b is the time fo the first sample (secs)
 * io_error = 0 upon success, -1 on failure
 * 
 */
/*****************************************************************************/

void WriteSAC1(char *filename,float *x,int npts,float dt,float b,int *io_error)
{
	struct SacHeader hd;

	/* set up a minimum SAC header for the output file */

	NewSacHeader(&hd);

	hd.npts = npts;
	hd.iftype = ITIME;
	hd.delta = dt;
	hd.b = b;
	hd.e = b + hd.npts * hd.delta;
	hd.iztype = IB;
	hd.leven = TRUE;

	WriteSACfile(filename,x,x,hd,io_error);

}

/*****************************************************************************/
/*
 * FUNCTION TO WRITE AN UNEVENLY_SPACED SAC FILE
 * 
 * Assuming that the header is NOT set up ---
 * 
 * filename is the name of the SAC file, x and y are the
 * arrays from which the data are written, io_error = 1 upon success, -1 on
 * failure
 */
/*****************************************************************************/
void WriteSAC2(char *filename, int n, float *x, float *y, int *io_error)
{
	FILE *strm;
	struct SacHeader hd;
	int nwritten;
	
	/* construct a minimum header */
	NewSacHeader(&hd);
	hd.npts = n;
	hd.iftype = IXY;
	hd.leven = FALSE;
	
	/* OPEN THE stream for writing the file */

	*io_error = kNO_ERROR;

	/* open the stream */
	strm = fopen(filename, "w");
	if (strm == NULL) {
		*io_error = kIO_ERROR;
		return;
	}	
	
	/* Now write the sac header to the file */
	nwritten = fwrite(&hd, sizeof(struct SacHeader), 1, strm);
		if (nwritten != 1) {
		*io_error = kIO_ERROR;
		fclose(strm);
		return;
	}
	/* Now write the sac data points y, then x, to the file */
	nwritten = fwrite(y, sizeof(float), n, strm);
		if (nwritten != n) {
		*io_error = kIO_ERROR;
		fclose(strm);
		return;
	}
	nwritten = fwrite(x, sizeof(float), n, strm);
		if (nwritten != n) {
		*io_error = kIO_ERROR;
		fclose(strm);
		return;
	}
	
	fclose(strm);
}

/*****************************************************************************/
/*
 * FORTRAN Stubs
 * 
 */
/*****************************************************************************/

/* this is a global value used for header access in fortran */

struct SacHeader ghdr;

void newhdr_()
{
	ghdr = sac_null;
	ghdr.nvhdr = 6;
}

void wsac1_(char *filename,float *x,int *npts,float *b,float *dt,int *io_error, long int flen)
{
	char token[FILENAME_MAX];

	trim(filename, token);
	WriteSAC1(token, x, *npts, *dt, *b, io_error);

}

void wsac2_(char *filename, float *y, int *npts, float *x, int *io_error, long int flen)
{

	char token[FILENAME_MAX];

	trim(filename, token);
	WriteSAC2(token, *npts, x, y, io_error);

}

void wsac0_(char *filename,float *x,float *y, int *io_error, long int flen)
{

	char token[FILENAME_MAX];

	trim(filename, token);

	/* write out the values in the global header, ghdr */
		
	/* ANYTHING you want in the header had better be 
	 * SET BEFORE calling this function.
	 * This includes the number of points, delta, etc.
	*/
	WriteSACfile(token, x, y, ghdr, io_error);

}
void rsac1_(char *filename,float *x,int *npts,float *b,float *dt,int *maxpts,int *io_error,long int flen)
{

	char token[FILENAME_MAX];

	trim(filename, token);
	ReadSACHeader(token, &ghdr, io_error);
		
	ReadSAC1(token,x,npts,dt,b,*maxpts,io_error);
	
}
void getfhv_(char *label, float *value, int *access_error, long int llen)
{	
	char thelabel[16], token[16];
	int i;
	
	*access_error = kNO_ERROR;

	/* convert to upper case */
	for(i=0;i<llen;i++)
	{
	    token[i] = toupper((int) label[i]);
	}
	token[i] = '\0';
	trim(token, thelabel);
//	printf("DEBUG: label: %s, token: %s, thelabel: %s.\n",label, token, thelabel);
	
	if(strcmp(thelabel,"B") == 0)
	{
             if(ghdr.b == sac_null.b){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.b;
	      return;
	}
	else if(strcmp(thelabel,"E") == 0)
	{
             if(ghdr.e == sac_null.e){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.e;
	      return;
	}
	else if(strcmp(thelabel,"DELTA") == 0)
	{
             if(ghdr.delta == sac_null.delta){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.delta;
	      return;
	}
	else if(strcmp(thelabel,"DEPMIN") == 0)
	{
             if(ghdr.depmin == sac_null.depmin){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.depmin;
	      return;
	}
	else if(strcmp(thelabel,"DEPMEN") == 0)
	{
             if(ghdr.depmen == sac_null.depmen){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.depmen;
	      return;
	}
	else if(strcmp(thelabel,"DEPMAX") == 0)
	{
             if(ghdr.depmax == sac_null.depmax){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.depmax;
	      return;
	}
	else if(strcmp(thelabel,"EVLA") == 0)
	{
             if(ghdr.evla == sac_null.evla){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.evla;
	      return;
	}
	else if(strcmp(thelabel,"EVLO") == 0)
	{
             if(ghdr.evlo == sac_null.evlo){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.evlo;
	      return;
	}
	else if(strcmp(thelabel,"EVDP") == 0)
	{
             if(ghdr.evdp == sac_null.evdp){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.evdp;
	      return;
	}
	else if(strcmp(thelabel,"STLA") == 0)
	{
             if(ghdr.stla == sac_null.stla){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.stla;
	      return;
	}
	else if(strcmp(thelabel,"STLO") == 0)
	{
             if(ghdr.stlo == sac_null.stlo){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.stlo;
	      return;
	}
	else if(strcmp(thelabel,"STEL") == 0)
	{
             if(ghdr.stel == sac_null.stel){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.stel;
	      return;
	}
	else if(strcmp(thelabel,"CMPAZ") == 0)
	{
             if(ghdr.cmpaz == sac_null.cmpaz){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.cmpaz;
	      return;
	}
	else if(strcmp(thelabel,"CMPINC") == 0)
	{
             if(ghdr.cmpinc == sac_null.cmpinc){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.cmpinc;
	      return;
	}
	else if(strcmp(thelabel,"DIST") == 0)
	{
             if(ghdr.dist == sac_null.dist){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.dist;
	      return;
	}
	else if(strcmp(thelabel,"GCARC") == 0)
	{
             if(ghdr.gcarc == sac_null.gcarc){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.gcarc;
	      return;
	}
	else if(strcmp(thelabel,"AZ") == 0)
	{
             if(ghdr.az == sac_null.az){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.az;
	      return;
	}
	else if(strcmp(thelabel,"BAZ") == 0)
	{
             if(ghdr.baz == sac_null.baz){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.baz;
	      return;
	}
	else if(strcmp(thelabel,"O") == 0)
	{
             if(ghdr.o == sac_null.o){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.o;
	      return;
	}
	else if(strcmp(thelabel,"A") == 0)
	{
             if(ghdr.a == sac_null.a){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.a;
	      return;
	}
	else if(strcmp(thelabel,"F") == 0)
	{
             if(ghdr.f == sac_null.f){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.f;
	      return;
	}
	else if(strcmp(thelabel,"T0") == 0)
	{
             if(ghdr.t0 == sac_null.t0){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.t0;
	      return;
	}
	else if(strcmp(thelabel,"T1") == 0)
	{
             if(ghdr.t1 == sac_null.t1){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.t1;
	      return;
	}
	else if(strcmp(thelabel,"T2") == 0)
	{
             if(ghdr.t2 == sac_null.t2){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.t2;
	      return;
	}
	else if(strcmp(thelabel,"T3") == 0)
	{
             if(ghdr.t3 == sac_null.t3){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.t3;
	      return;
	}
	else if(strcmp(thelabel,"T4") == 0)
	{
             if(ghdr.t4 == sac_null.t4){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.t4;
	      return;
	}
	else if(strcmp(thelabel,"T5") == 0)
	{
             if(ghdr.t5 == sac_null.t5){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.t5;
	      return;
	}
	else if(strcmp(thelabel,"T6") == 0)
	{
             if(ghdr.t6 == sac_null.t6){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.t6;
	      return;
	}
	else if(strcmp(thelabel,"T7") == 0)
	{
             if(ghdr.t7 == sac_null.t7){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.t7;
	      return;
	}
	else if(strcmp(thelabel,"T8") == 0)
	{
             if(ghdr.t8 == sac_null.t8){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.t8;
	      return;
	}
	else if(strcmp(thelabel,"T9") == 0)
	{
             if(ghdr.t9 == sac_null.t9){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.t9;
	      return;
	}
	else if(strcmp(thelabel,"USER0") == 0)
	{
             if(ghdr.user0 == sac_null.user0){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.user0;
	      return;
	}
	else if(strcmp(thelabel,"USER1") == 0)
	{
             if(ghdr.user1 == sac_null.user1){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.user1;
	      return;
	}
	else if(strcmp(thelabel,"USER2") == 0)
	{
             if(ghdr.user2 == sac_null.user2){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.user2;
	      return;
	}
	else if(strcmp(thelabel,"USER3") == 0)
	{
             if(ghdr.user3 == sac_null.user3){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.user3;
	      return;
	}
	else if(strcmp(thelabel,"USER4") == 0)
	{
             if(ghdr.user4 == sac_null.user4){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.user4;
	      return;
	}
	else if(strcmp(thelabel,"USER5") == 0)
	{
             if(ghdr.user5 == sac_null.user5){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.user5;
	      return;
	}
	else if(strcmp(thelabel,"USER6") == 0)
	{
             if(ghdr.user6 == sac_null.user6){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.user6;
	      return;
	}
	else if(strcmp(thelabel,"USER7") == 0)
	{
             if(ghdr.user7 == sac_null.user7){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.user7;
	      return;
	}
	else if(strcmp(thelabel,"USER8") == 0)
	{
             if(ghdr.user8 == sac_null.user8){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.user8;
	      return;
	}
	else if(strcmp(thelabel,"USER9") == 0)
	{
             if(ghdr.user9 == sac_null.user9){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.user9;
	      return;
	}
	else if(strcmp(thelabel,"XMAXIMUM") == 0)
	{
             if(ghdr.xmaximum == sac_null.xmaximum){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.xmaximum;
	      return;
	}
	else if(strcmp(thelabel,"XMINIMUM") == 0)
	{
             if(ghdr.xminimum == sac_null.xminimum){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.xminimum;
	      return;
	}
	else if(strcmp(thelabel,"YMAXIMUM") == 0)
	{
             if(ghdr.ymaximum == sac_null.ymaximum){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.ymaximum;
	      return;
	}
	else if(strcmp(thelabel,"YMINIMUM") == 0)
	{
             if(ghdr.yminimum == sac_null.yminimum){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.yminimum;
	      return;
	}
	else
	{
	       *value = -12345;
	       printf("Warning %s not matched in getfhv\n",thelabel);
               *access_error = kIO_ERROR;
	       return;
	}	

}
void getkhv_(char *label, char *value, int *access_error, long int llen, long int vlen)
{	
	char thelabel[16], token[16];
	int i;
	
	*access_error = kNO_ERROR;
	/* convert to upper case */
	for(i=0;i<llen;i++)
	{
	    token[i] = toupper((int) label[i]);
	}
	token[i] = '\0';
	trim(token, thelabel);
	
	if(strcmp(thelabel,"KSTNM") == 0)
	{
              if(strcmp(ghdr.kstnm,sac_null.kstnm)==0){
                *access_error = NOT_DEF;
                return;
              }
	      strncpy(value,ghdr.kstnm,8);
	      return;
	}
	else if(strcmp(thelabel,"KEVNM") == 0)
	{
              if(strcmp(ghdr.kevnm,sac_null.kevnm)==0){
                *access_error = NOT_DEF;
                return;
              }
	      strncpy(value,ghdr.kevnm,16);
	      return;
	}
	else if(strcmp(thelabel,"KHOLE") == 0)
	{
              if(strcmp(ghdr.khole,sac_null.khole)==0){
                *access_error = NOT_DEF;
                return;
              }
	      strncpy(value,ghdr.khole,8);
	      return;
	}
	else if(strcmp(thelabel,"KO") == 0)
	{
              if(strcmp(ghdr.ko,sac_null.ko)==0){
                *access_error = NOT_DEF;
                return;
              }
	      strncpy(value,ghdr.ko,8);
	      return;
	}
	else if(strcmp(thelabel,"KA") == 0)
	{
              if(strcmp(ghdr.ka,sac_null.ka)==0){
                *access_error = NOT_DEF;
                return;
              }
	      strncpy(value,ghdr.ka,8);
	      return;
	}
	else if(strcmp(thelabel,"KT0") == 0)
	{
              if(strcmp(ghdr.kt0,sac_null.kt0)==0){
                *access_error = NOT_DEF;
                return;
              }
	      strncpy(value,ghdr.kt0,8);
	      return;
	}
	else if(strcmp(thelabel,"KT1") == 0)
	{
              if(strcmp(ghdr.kt1,sac_null.kt1)==0){
                *access_error = NOT_DEF;
                return;
              }
	      strncpy(value,ghdr.kt1,8);
	      return;
	}
	else if(strcmp(thelabel,"KT2") == 0)
	{
              if(strcmp(ghdr.kt2,sac_null.kt2)==0){
                *access_error = NOT_DEF;
                return;
              }
	      strncpy(value,ghdr.kt2,8);
	      return;
	}
	else if(strcmp(thelabel,"KT3") == 0)
	{
              if(strcmp(ghdr.kt3,sac_null.kt3)==0){
                *access_error = NOT_DEF;
                return;
              }
	      strncpy(value,ghdr.kt3,8);
	      return;
	}
	else if(strcmp(thelabel,"KT4") == 0)
	{
              if(strcmp(ghdr.kt4,sac_null.kt4)==0){
                *access_error = NOT_DEF;
                return;
              }
	      strncpy(value,ghdr.kt4,8);
	      return;
	}
	else if(strcmp(thelabel,"KT5") == 0)
	{
              if(strcmp(ghdr.kt5,sac_null.kt5)==0){
                *access_error = NOT_DEF;
                return;
              }
	      strncpy(value,ghdr.kt5,8);
	      return;
	}
	else if(strcmp(thelabel,"KT6") == 0)
	{
              if(strcmp(ghdr.kt6,sac_null.kt6)==0){
                *access_error = NOT_DEF;
                return;
              }
	      strncpy(value,ghdr.kt6,8);
	      return;
	}
	else if(strcmp(thelabel,"KT7") == 0)
	{
              if(strcmp(ghdr.kt7,sac_null.kt7)==0){
                *access_error = NOT_DEF;
                return;
              }
	      strncpy(value,ghdr.kt7,8);
	      return;
	}
	else if(strcmp(thelabel,"KT8") == 0)
	{
              if(strcmp(ghdr.kt8,sac_null.kt8)==0){
                *access_error = NOT_DEF;
                return;
              }
	      strncpy(value,ghdr.kt8,8);
	      return;
	}
	else if(strcmp(thelabel,"KT9") == 0)
	{
              if(strcmp(ghdr.kt9,sac_null.kt9)==0){
                *access_error = NOT_DEF;
                return;
              }
	      strncpy(value,ghdr.kt9,8);
	      return;
	}
	else if(strcmp(thelabel,"KF") == 0)
	{
              if(strcmp(ghdr.kf,sac_null.kf)==0){
                *access_error = NOT_DEF;
                return;
              }
	      strncpy(value,ghdr.kf,8);
	      return;
	}
	else if(strcmp(thelabel,"KUSER0") == 0)
	{
              if(strcmp(ghdr.kuser0,sac_null.kuser0)==0){
                *access_error = NOT_DEF;
                return;
              }
	      strncpy(value,ghdr.kuser0,8);
	      return;
	}
	else if(strcmp(thelabel,"KUSER1") == 0)
	{
              if(strcmp(ghdr.kuser1,sac_null.kuser1)==0){
                *access_error = NOT_DEF;
                return;
              }
	      strncpy(value,ghdr.kuser1,8);
	      return;
	}
	else if(strcmp(thelabel,"KUSER2") == 0)
	{
              if(strcmp(ghdr.kuser2,sac_null.kuser2)==0){
                *access_error = NOT_DEF;
                return;
              }
	      strncpy(value,ghdr.kuser2,8);
	      return;
	}
	else if(strcmp(thelabel,"KCMPNM") == 0)
	{
              if(strcmp(ghdr.kcmpnm,sac_null.kcmpnm)==0){
                *access_error = NOT_DEF;
                return;
              }
	      strncpy(value,ghdr.kcmpnm,8);
	      return;
	}
	else if(strcmp(thelabel,"KNETWK") == 0)
	{
              if(strcmp(ghdr.knetwk,sac_null.knetwk)==0){
                *access_error = NOT_DEF;
                return;
              }
	      strncpy(value,ghdr.knetwk,8);
	      return;
	}
	else if(strcmp(thelabel,"KDATRD") == 0)
	{
              if(strcmp(ghdr.kdatrd,sac_null.kdatrd)==0){
                *access_error = NOT_DEF;
                return;
              }
	      strncpy(value,ghdr.kdatrd,8);
	      return;
	}
	else if(strcmp(thelabel,"KINST") == 0)
	{
              if(strcmp(ghdr.kinst,sac_null.kinst)==0){
                *access_error = NOT_DEF;
                return;
              }
	      strncpy(value,ghdr.kinst,8);
	      return;
	}
	else
	{
	       strncpy(value,"-12345",6);
	       printf("Warning %s not matched in getfhv\n",thelabel);
	       *access_error = kIO_ERROR;
	       return;
	}	

}

void setkhv_(char *label, char *value, int *access_error, long int llen, long int vlen)
{	
	char thelabel[16], thevalue[16];
	char labeltoken[16], valuetoken[16];
	int i;
	
	/* convert to upper case */
	for(i=0;i<llen;i++)
	{
	    labeltoken[i] = toupper((int) label[i]);
	}
	labeltoken[i] = '\0';
	for(i=0;i<vlen;i++)
	{
	    valuetoken[i] = toupper((int) value[i]);
	}
	valuetoken[i] = '\0';

	trim(labeltoken, thelabel);
	trim(valuetoken, thevalue);
	
	
	if(strcmp(thelabel,"KSTNM") == 0)
	{
	      strncpy(ghdr.kstnm,thevalue,8);
	}
	else if(strcmp(thelabel,"KEVNM") == 0)
	{
	      strncpy(ghdr.kevnm,thevalue,16);
	}
	else if(strcmp(thelabel,"KHOLE") == 0)
	{
	      strncpy(ghdr.khole,thevalue,8);
	}
	else if(strcmp(thelabel,"KO") == 0)
	{
	      strncpy(ghdr.ko,thevalue,8);
	}
	else if(strcmp(thelabel,"KA") == 0)
	{
	      strncpy(ghdr.ka,thevalue,8);
	}
	else if(strcmp(thelabel,"KT0") == 0)
	{
	      strncpy(ghdr.kt0,thevalue,8);
	}
	else if(strcmp(thelabel,"KT1") == 0)
	{
	      strncpy(ghdr.kt1,thevalue,8);
	}
	else if(strcmp(thelabel,"KT2") == 0)
	{
	      strncpy(ghdr.kt2,thevalue,8);
	}
	else if(strcmp(thelabel,"KT3") == 0)
	{
	      strncpy(ghdr.kt3,thevalue,8);
	}
	else if(strcmp(thelabel,"KT4") == 0)
	{
	      strncpy(ghdr.kt4,thevalue,8);
	}
	else if(strcmp(thelabel,"KT5") == 0)
	{
	      strncpy(ghdr.kt5,thevalue,8);
	}
	else if(strcmp(thelabel,"KT6") == 0)
	{
	      strncpy(ghdr.kt6,thevalue,8);
	}
	else if(strcmp(thelabel,"KT7") == 0)
	{
	      strncpy(ghdr.kt7,thevalue,8);
	}
	else if(strcmp(thelabel,"KT8") == 0)
	{
	      strncpy(ghdr.kt8,thevalue,8);
	}
	else if(strcmp(thelabel,"KT9") == 0)
	{
	      strncpy(ghdr.kt9,thevalue,8);
	}
	else if(strcmp(thelabel,"KF") == 0)
	{
	      strncpy(ghdr.kf,thevalue,8);
	}
	else if(strcmp(thelabel,"KUSER0") == 0)
	{
	      strncpy(ghdr.kuser0,thevalue,8);
	}
	else if(strcmp(thelabel,"KUSER1") == 0)
	{
	      strncpy(ghdr.kuser1,thevalue,8);
	}
	else if(strcmp(thelabel,"KUSER2") == 0)
	{
	      strncpy(ghdr.kuser2,thevalue,8);
	}
	else if(strcmp(thelabel,"KCMPNM") == 0)
	{
	      strncpy(ghdr.kcmpnm,thevalue,8);
	}
	else if(strcmp(thelabel,"KNETWK") == 0)
	{
	      strncpy(ghdr.knetwk,thevalue,8);
	}
	else if(strcmp(thelabel,"KDATRD") == 0)
	{
	      strncpy(ghdr.kdatrd,thevalue,8);
	}
	else if(strcmp(thelabel,"KINST") == 0)
	{
	      strncpy(ghdr.kinst,thevalue,8);
	}
	else
	{
	       printf("Warning %s not matched in setkhv\n",thelabel);
               *access_error = kIO_ERROR;
	       return;
	}	

}


void setfhv_(char *label, float *value, int *access_error, long int llen)
{	
	char thelabel[16], token[16];
	int i;
	
	/* convert to upper case */
	for(i=0;i<llen;i++)
	{
	    token[i] = toupper((int) label[i]);
	}
	token[i] = '\0';

	trim(token, thelabel);
	
	
	if(strcmp(thelabel,"B") == 0)
	{
	     ghdr.b = *value ;
	      return;
	}
	else if(strcmp(thelabel,"E") == 0)
	{
	     ghdr.e = *value;
	      return;
	}
	else if(strcmp(thelabel,"DELTA") == 0)
	{
	     ghdr.delta = *value;
	      return;
	}
	else if(strcmp(thelabel,"DEPMIN") == 0)
	{
	     ghdr.depmin = *value;
	      return;
	}
	else if(strcmp(thelabel,"DEPMEN") == 0)
	{
	     ghdr.depmen = *value;
	      return;
	}
	else if(strcmp(thelabel,"DEPMAX") == 0)
	{
	     ghdr.depmax = *value;
	      return;
	}
	else if(strcmp(thelabel,"EVLA") == 0)
	{
	      ghdr.evla = *value;
	      return;
	}
	else if(strcmp(thelabel,"EVLO") == 0)
	{
	     ghdr.evlo = *value;
	      return;
	}
	else if(strcmp(thelabel,"EVDP") == 0)
	{
	     ghdr.evdp = *value;
	      return;
	}
	else if(strcmp(thelabel,"STLA") == 0)
	{
	     ghdr.stla = *value;
	      return;
	}
	else if(strcmp(thelabel,"STLO") == 0)
	{
	     ghdr.stlo = *value;
	      return;
	}
	else if(strcmp(thelabel,"STEL") == 0)
	{
	     ghdr.stel = *value;
	      return;
	}
	else if(strcmp(thelabel,"CMPAZ") == 0)
	{
	     ghdr.cmpaz = *value;
	      return;
	}
	else if(strcmp(thelabel,"CMPINC") == 0)
	{
	     ghdr.cmpinc = *value;
	      return;
	}
	else if(strcmp(thelabel,"DIST") == 0)
	{
	     ghdr.dist = *value;
	      return;
	}
	else if(strcmp(thelabel,"GCARC") == 0)
	{
	     ghdr.gcarc = *value;
	      return;
	}
	else if(strcmp(thelabel,"AZ") == 0)
	{
	     ghdr.az = *value;
	      return;
	}
	else if(strcmp(thelabel,"BAZ") == 0)
	{
	     ghdr.baz = *value;
	      return;
	}
	else if(strcmp(thelabel,"O") == 0)
	{
	     ghdr.o = *value;
	      return;
	}
	else if(strcmp(thelabel,"A") == 0)
	{
	     ghdr.a = *value;
	      return;
	}
	else if(strcmp(thelabel,"F") == 0)
	{
	     ghdr.f = *value;
	      return;
	}
	else if(strcmp(thelabel,"T0") == 0)
	{
	     ghdr.t0 = *value;
	      return;
	}
	else if(strcmp(thelabel,"T1") == 0)
	{
	     ghdr.t1 = *value;
	      return;
	}
	else if(strcmp(thelabel,"T2") == 0)
	{
	     ghdr.t2 = *value;
	      return;
	}
	else if(strcmp(thelabel,"T3") == 0)
	{
	     ghdr.t3 = *value;
	      return;
	}
	else if(strcmp(thelabel,"T4") == 0)
	{
	     ghdr.t4 = *value;
	      return;
	}
	else if(strcmp(thelabel,"T5") == 0)
	{
	     ghdr.t5 = *value;
	      return;
	}
	else if(strcmp(thelabel,"T6") == 0)
	{
	     ghdr.t6 = *value;
	      return;
	}
	else if(strcmp(thelabel,"T7") == 0)
	{
	     ghdr.t7 = *value;
	      return;
	}
	else if(strcmp(thelabel,"T8") == 0)
	{
	     ghdr.t8 = *value;
	      return;
	}
	else if(strcmp(thelabel,"T9") == 0)
	{
	     ghdr.t9 = *value;
	      return;
	}
	else if(strcmp(thelabel,"USER0") == 0)
	{
	     ghdr.user0 = *value;
	      return;
	}
	else if(strcmp(thelabel,"USER1") == 0)
	{
	     ghdr.user1 = *value;
	      return;
	}
	else if(strcmp(thelabel,"USER2") == 0)
	{
	     ghdr.user2 = *value;
	      return;
	}
	else if(strcmp(thelabel,"USER3") == 0)
	{
	     ghdr.user3 = *value;
	      return;
	}
	else if(strcmp(thelabel,"USER4") == 0)
	{
	     ghdr.user4 = *value;
	      return;
	}
	else if(strcmp(thelabel,"USER5") == 0)
	{
	     ghdr.user5 = *value;
	      return;
	}
	else if(strcmp(thelabel,"USER6") == 0)
	{
	     ghdr.user6 = *value;
	      return;
	}
	else if(strcmp(thelabel,"USER7") == 0)
	{
	     ghdr.user7 = *value;
	      return;
	}
	else if(strcmp(thelabel,"USER8") == 0)
	{
	     ghdr.user8 = *value;
	      return;
	}
	else if(strcmp(thelabel,"USER9") == 0)
	{
	     ghdr.user9 = *value;
	      return;
	}
	else if(strcmp(thelabel,"XMAXIMUM") == 0)
	{
	     ghdr.xmaximum = *value;
	      return;
	}
	else if(strcmp(thelabel,"XMINIMUM") == 0)
	{
	     ghdr.xminimum = *value;
	      return;
	}
	else if(strcmp(thelabel,"YMAXIMUM") == 0)
	{
	     ghdr.ymaximum = *value;
	      return;
	}
	else if(strcmp(thelabel,"YMINIMUM") == 0)
	{
	     ghdr.yminimum = *value;
	      return;
	}
	else
	{
	       printf("Warning %s not matched in setfhv\n",thelabel);
               *access_error = kIO_ERROR;
	       return;
	}	

}


void getnhv_(char *label, int *value, int *access_error, long int llen)
{	
	char thelabel[16], token[16];
	int i;
	
        *access_error = kNO_ERROR;
	/* convert to upper case */
	for(i=0;i<llen;i++)
	{
	    token[i] = toupper((int) label[i]);
	}
	token[i] = '\0';

	trim(token, thelabel);
	
	
	if(strcmp(thelabel,"NZYEAR") == 0)
	{
             if(ghdr.nzyear == sac_null.nzyear){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.nzyear;
	      return;
	}
	else if(strcmp(thelabel,"NZJDAY") == 0)
	{
             if(ghdr.nzjday == sac_null.nzjday){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.nzjday;
	      return;
	}
	else if(strcmp(thelabel,"NZHOUR") == 0)
	{
             if(ghdr.nzhour == sac_null.nzhour){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.nzhour;
	      return;
	}
	else if(strcmp(thelabel,"NZMIN") == 0)
	{
             if(ghdr.nzmin == sac_null.nzmin){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.nzmin;
	      return;
	}
	else if(strcmp(thelabel,"NZSEC") == 0)
	{
             if(ghdr.nzsec == sac_null.nzsec){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.nzsec;
	      return;
	}
	else if(strcmp(thelabel,"NZMSEC") == 0)
	{
             if(ghdr.nzmsec == sac_null.nzmsec){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.nzmsec;
	      return;
	}
	else if(strcmp(thelabel,"NVHDR") == 0)
	{
             if(ghdr.nvhdr == sac_null.nvhdr){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.nvhdr;
	      return;
	}
	else if(strcmp(thelabel,"NPTS") == 0)
	{
             if(ghdr.npts == sac_null.npts){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.npts;
	      return;
	}
	else if(strcmp(thelabel,"NXSIZE") == 0)
	{
             if(ghdr.nxsize == sac_null.nxsize){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.nxsize;
	      return;
	}
	else if(strcmp(thelabel,"NYSIZE") == 0)
	{
             if(ghdr.nysize == sac_null.nysize){
                *access_error = NOT_DEF;
                return;
             }
	     *value = ghdr.nysize;
	      return;
	}
	else
	{
	       printf("Warning %s not matched in getnhv\n",thelabel);
	       *value = -12345;
               *access_error = kIO_ERROR;
	       return;
	}	

}


void setnhv_(char *label, int *value, int *access_error, long int llen)
{	
	char thelabel[16], token[16];
	int i;
	
	/* convert to upper case */
	for(i=0;i<llen;i++)
	{
	    token[i] = toupper((int) label[i]);
	}
	token[i] = '\0';

	trim(token, thelabel);
	
	
	if(strcmp(thelabel,"NZYEAR") == 0)
	{
	     ghdr.nzyear = *value;
	      return;
	}
	else if(strcmp(thelabel,"NZJDAY") == 0)
	{
	     ghdr.nzjday = *value;
	      return;
	}
	else if(strcmp(thelabel,"NZHOUR") == 0)
	{
	     ghdr.nzhour = *value;
	      return;
	}
	else if(strcmp(thelabel,"NZMIN") == 0)
	{
	     ghdr.nzmin = *value;
	      return;
	}
	else if(strcmp(thelabel,"NZSEC") == 0)
	{
	     ghdr.nzsec = *value;
	      return;
	}
	else if(strcmp(thelabel,"NZMSEC") == 0)
	{
	     ghdr.nzmsec = *value;
	      return;
	}
	else if(strcmp(thelabel,"NVHDR") == 0)
	{
	     ghdr.nvhdr = *value;
	      return;
	}
	else if(strcmp(thelabel,"NPTS") == 0)
	{
	     ghdr.npts = *value;
	      return;
	}
	else if(strcmp(thelabel,"NXSIZE") == 0)
	{
	     ghdr.nxsize = *value;
	      return;
	}
	else if(strcmp(thelabel,"NYSIZE") == 0)
	{
	     ghdr.nysize = *value;
	      return;
	}
	else
	{
	       printf("Warning: %s not found in setnhv\n",thelabel);
	       return;
	}	

}


void getihv_(char *label, char *value, int *access_error, long int llen, long int vlen)
{
	char thelabel[16];
	char labeltok[16];
	int i;
	
	/* convert to upper case */
	for(i=0;i<llen;i++)
	{
	    labeltok[i] = toupper((int) label[i]);
	}
	labeltok[i] = '\0';

	trim(labeltok, thelabel);
	

	if(strcmp(thelabel,"IFTYPE") == 0)
	{
		switch (ghdr.iftype) {
			case ITIME:
	      			strncpy(value,"ITIME",5);
				break;
			case IRLIM:
	      			strncpy(value,"IRLIM",5);
				break;
			case IAMPH:
	      			strncpy(value,"IAMPH",5);
				break;
			case IXY:
	      			strncpy(value,"IXY",3);
				break;
			case IXYZ:
	      			strncpy(value,"IXYZ",4);
				break;
			default:
				printf("Unknown IFTYPE %ld in getihv.\n",ghdr.iftype);
				*access_error = NOT_EXIST;
				break;
		}
		return;
	}
	else if (strcmp(thelabel,"IDEP") == 0)
	{
	switch (ghdr.idep) {
			case IUNKN:
	      			strncpy(value,"IUNKN",5);
				break;
			case IDISP:
	      			strncpy(value,"IDISP",5);
				break;
			case IVEL:
	      			strncpy(value,"IVEL",4);
				break;
			case IVOLTS:
	      			strncpy(value,"IVOLTS",6);
				break;
			case IACC:
	      			strncpy(value,"IACC",4);
				break;
			default:
				printf("Unknown IDEP %ld in getihv.\n",ghdr.idep);
				*access_error = NOT_EXIST;
				break;
		}
		return;

	}
else if (strcmp(thelabel,"IZTYPE") == 0)
	{
	switch (ghdr.iztype) {
			case IUNKN:
	      			strncpy(value,"IUNKN",5);
				break;
			case IB:
	      			strncpy(value,"IB",2);
				break;
			case IDAY:
	      			strncpy(value,"IDAY",4);
				break;
			case IO:
	      			strncpy(value,"IO",2);
				break;
			case IA:
	      			strncpy(value,"IA",2);
				break;
			case IT1:
	      			strncpy(value,"IT1",3);
				break;
			case IT2:
	      			strncpy(value,"IT2",3);
				break;
			case IT3:
	      			strncpy(value,"IT3",3);
				break;
			case IT4:
	      			strncpy(value,"IT4",3);
				break;
			case IT5:
	      			strncpy(value,"IT5",3);
				break;
			case IT6:
	      			strncpy(value,"IT6",3);
				break;
			case IT7:
	      			strncpy(value,"IT7",3);
				break;
			case IT8:
	      			strncpy(value,"IT8",3);
				break;
			case IT9:
	      			strncpy(value,"IT9",3);
				break;
			case IT0:
	      			strncpy(value,"IT0",3);
				break;
			default:
				printf("Unknown IZTYPE %ld in getihv.\n",ghdr.iztype);
				*access_error = NOT_EXIST;
				break;
		}
		return;

	}
	else
	{
		printf("Warning: %s not found in getihv\n",thelabel);
		*access_error = NOT_EXIST;
		return;
	}

}

void setihv_(char *label, char *value, int *access_error, long int llen, long int vlen)
{	
	char thelabel[16], thevalue[16];
	char labeltok[16], valuetok[16];
	int i;
	
	/* convert to upper case */
	for(i=0;i<llen;i++)
	{
	    labeltok[i] = toupper((int) label[i]);
	}
	labeltok[i] = '\0';
	
	/* convert to upper case */
	for(i=0;i<vlen;i++)
	{
	    valuetok[i] = toupper((int) value[i]);
	}
	valuetok[i] = '\0';

	trim(labeltok, thelabel);
	trim(valuetok, thevalue);

	
	if(strcmp(thelabel,"IFTYPE") == 0)
	{
		if(strcmp(thevalue,"ITIME") == 0) 
	     		ghdr.iftype = ITIME;
		else if (strcmp(thevalue,"IRLIM") == 0)
			ghdr.iftype = IRLIM;
		else if (strcmp(thevalue,"IAMPH") == 0)
			ghdr.iftype = IAMPH;
		else if (strcmp(thevalue,"IXY") == 0)
			ghdr.iftype = IXY;
		else if (strcmp(thevalue,"IXYZ") == 0)
			ghdr.iftype = IXYZ;
		else{
			ghdr.iftype = 12345;
			printf("Unknown value %s of IFTIME in setihv.\n",thevalue);
		}
	     return;
	}
	else if(strcmp(thelabel,"IDEP") == 0)
	{
		if(strcmp(thevalue,"IUNKN") == 0) 
	     		ghdr.idep = IUNKN;
		else if (strcmp(thevalue,"IDISP") == 0)
			ghdr.idep = IDISP;
		else if (strcmp(thevalue,"IVEL") == 0)
			ghdr.idep = IVEL;
		else if (strcmp(thevalue,"IVOLTS") == 0)
			ghdr.idep = IVOLTS;
		else if (strcmp(thevalue,"IACC") == 0)
			ghdr.idep = IACC;
		else{
			ghdr.idep = 12345;
			printf("Unknown value %s of IDEP in setihv.\n",thevalue);
		}
	     return;
	}
else if(strcmp(thelabel,"IZTYPE") == 0)
	{
		if(strcmp(thevalue,"IUNKN") == 0) 
	     		ghdr.iztype = IUNKN;
		else if (strcmp(thevalue,"IB") == 0)
			ghdr.iztype = IB;
		else if (strcmp(thevalue,"IDAY") == 0)
			ghdr.iztype = IDAY;
		else if (strcmp(thevalue,"IO") == 0)
			ghdr.iztype = IO;
		else if (strcmp(thevalue,"IA") == 0)
			ghdr.iztype = IA;
		else if (strcmp(thevalue,"IT0") == 0)
			ghdr.iztype = IT0;
		else if (strcmp(thevalue,"IT1") == 0)
			ghdr.iztype = IT1;
		else if (strcmp(thevalue,"IT2") == 0)
			ghdr.iztype = IT2;
		else if (strcmp(thevalue,"IT3") == 0)
			ghdr.iztype = IT3;
		else if (strcmp(thevalue,"IT4") == 0)
			ghdr.iztype = IT4;
		else if (strcmp(thevalue,"IT5") == 0)
			ghdr.iztype = IT5;
		else if (strcmp(thevalue,"IT6") == 0)
			ghdr.iztype = IT6;
		else if (strcmp(thevalue,"IT7") == 0)
			ghdr.iztype = IT7;
		else if (strcmp(thevalue,"IT8") == 0)
			ghdr.iztype = IT8;
		else if (strcmp(thevalue,"IT9") == 0)
			ghdr.iztype = IT9;
		else{
			ghdr.iztype = 12345;
			printf("Unknown value %s of IZTYPE in setihv.\n",thevalue);
		}
	     return;
	}

	else
	{
	       printf("Warning: %s not found in setihv\n",thelabel);
	       return;
	}	

}



void setlhv_(char *label, int *value, int *access_error, long int llen)
{	
	char thelabel[16];
	int i;
	
	/* convert to upper case */
	for(i=0;i<llen;i++)
	{
	    thelabel[i] = toupper((int) label[i]);
	}
	thelabel[i] = '\0';
	
	
	if(strcmp(thelabel,"LEVEN") == 0)
	{
	     ghdr.leven = *value;
	     return;
	}
	else if(strcmp(thelabel,"LCALDA") == 0)
	{
	     ghdr.lcalda = *value;
	     return;
	}
	else
	{
	       printf("Warning: %s not found in setlhv\n",thelabel);
	       return;
	}	

}

/* ************************************************************ */
/* Remove the trailing blanks etc. in fortran-passed strings    */
/* Return the length of the resulting string                    */
/* ************************************************************ */

int trim(char *s, char *t)
{
	int i, original_length;

	original_length=strlen(s);

	for(i=0;i<=original_length;i++)
		if(isspace(s[i])){
			t[i]='\0';
			return i;
		}
		else
			t[i]=s[i];

	return original_length;

}

#if 0
/*****************************************************************************/
/* write a sac spectral file ??? */
void 
wrtsac1(FO, dt, ns, ar)
	char            FO[];
	int             ns;
	float           dt, ar[];
{
	FILE           *strm;
	struct SacHeader      hd;
	void            NewSacHeader();
	int             i;
	float           zero;
	zero = 0.;
	NewSacHeader(&hd);
	hd.npts = ns;
	hd.iftype = IRLIM;
	hd.iztype = IB;
	hd.b = zero;
	hd.e = 0.5 / dt;
	hd.internal2 = 0.;
	hd.internal3 = dt;
	hd.internal7 = ns;
	ns = ns / 2;
	hd.delta = hd.e / ns;
	hd.leven = TRUE;
	strm = fopen(FO, "w");
	fwrite(&hd, sizeof(struct SacHeader), 1, strm);
	for (i = 0; i < ns; i++)
		fwrite(&ar[2 * i], sizeof(float), 1, strm);
	fwrite(&ar[1], sizeof(float), 1, strm);
	for (i = ns - 1; i > 0; i--)
		fwrite(&ar[2 * i], sizeof(float), 1, strm);
	fwrite(&zero, sizeof(float), 1, strm);
	for (i = 1; i < ns; i++)
		fwrite(&ar[2 * i + 1], sizeof(float), 1, strm);
	fwrite(&zero, sizeof(float), 1, strm);
	for (i = ns - 1; i > 0; i--) {
		zero = -ar[2 * i + 1];
		fwrite(&zero, sizeof(float), 1, strm);
	}
	fclose(strm);
}

/*****************************************************************************/
/* Write an IMAGE sac file ??? */
void 
wrtsac3(FO, nx, ny, xmin, xmax, ymin, ymax, ar)
	char            FO[];
	int             nx, ny;
	float           xmin, xmax, ymin, ymax, ar[];
{
	FILE           *strm;
	struct SacHeader      hd;
	void            NewSacHeader();
	NewSacHeader(&hd);
	hd.npts = nx * ny;
	hd.iftype = INIV51;
	hd.delta = 1.;
	hd.iztype = IB;
	hd.b = 0.;
	hd.e = hd.npts * hd.delta;
	hd.leven = TRUE;
	hd.unused2 = xmin;
	hd.unused3 = xmax;
	hd.unused4 = ymin;
	hd.unused5 = ymax;
	hd.unused13 = nx;
	hd.unused14 = ny;
	strm = fopen(FO, "w");
	fwrite(&hd, sizeof(struct SacHeader), 1, strm);
	fwrite(ar, hd.npts * sizeof(float), 1, strm);
	fclose(strm);
}
/*****************************************************************************/
#endif
