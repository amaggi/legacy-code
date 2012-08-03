           integer*4 lyrs,lsd,lgrm,lgrm2
	   parameter (lyrs=200, lsd=4,lgrm=2050,lgrm2=lgrm*2)
	   integer*4 maxrecs,maxpee,maxcol,maxrow
	   parameter (maxcol=lyrs*3,maxpee=10)
	   parameter (maxrecs=8,maxrow=maxrecs*lgrm+maxcol+maxpee)
	   integer*4 maxsect
	   parameter (maxsect=maxcol*lgrm)
	   integer*4 maxbrn,maxbrn2
	   parameter (maxbrn=50,maxbrn2=maxbrn*2)
