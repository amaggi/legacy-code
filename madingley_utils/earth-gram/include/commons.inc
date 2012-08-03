	  real*8 sig,del,gam,fmom,sol,tconst
	  real*4 xlength
  	  common/source/sig,del,gam,fmom,sol(6),tconst,xlength
 
          real*4 d0,th,ph,sec
	  integer*4 jy,jd,jh,jm
          common/dhead/d0,th,ph,jy,jd,jh,jm,sec

          real*4 to,dist,t1,p1
	  integer*4 ivel,mdmin,mdmx,nsdp,init
          common/vary/ivel,to,mdmin,mdmx,nsdp,dist,init,t1,p1

          real*4 sdep
	  integer*4 npts,nsrce,jcom,nbran,IRAY,ILOV
          common/same/npts,nsrce,jcom,nbran,sdep(lsd),IRAY,ILOV

          complex resp
          common/resp/resp(lgrm/2+1)

 	  character  sta*4,chn*4,typ*4
	  real*4 ss,dt
	  integer*4 nscan,iy,id,ih,im
          common/myhed/nscan,sta,chn,typ,iy,id,ih,im,ss,dt
