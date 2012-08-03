c
c $Id: slresp.f,v 1.1.1.1 2002/07/12 11:15:20 maggi Exp $
c $Log: slresp.f,v $
c Revision 1.1.1.1  2002/07/12 11:15:20  maggi
c
c
c Revision 1.1  2002/05/23 10:28:35  maggi
c Initial revision
c
c
       complex function slresp(w)
       double precision pi,a0,s0
       parameter(nzero=5,npole=11,pi=3.14159265358979)
       double complex z,p,zero(nzero),pole(npole)
       data a0,s0/0.0276,5D+08/
       data zero/(0.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.)/
       data pole/(-3.770000D-01,  1.830000D-01),
     *           (-3.770000D-01, -1.830000D-01),
     *           (-6.540000D-01,  0.000000D+00),
     *           (-2.320000D-01,  0.000000D+00),
     *           (-2.320000D-01,  0.000000D+00),
     *           (-2.320000D-01,  0.000000D+00),
     *           (-3.280000D-01,  0.000000D+00),
     *           (-3.280000D-01,  0.000000D+00),
     *           (-3.280000D-01,  0.000000D+00),
     *           (-2.140000D-02,  0.000000D+00),
     *           (-2.140000D-02,  0.000000D+00)/



       z=dcmplx(0,w)-zero(1)
       do 10 i=2,nzero
       z=z*(dcmplx(0,w)-zero(i))
 10    continue

       p=dcmplx(0,w)-pole(1)
       do 20 j=2,npole
       p=p*(dcmplx(0,w)-pole(j))
 20    continue

c      gains are counts/meter
c      change polarization to VRT epicentral coordinates
       slresp=-a0*s0*z/p


       return
       end
