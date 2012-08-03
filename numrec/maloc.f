      FUNCTION maloc(len)
      INTEGER maloc,len,NG,MEMLEN
      PARAMETER (NG=5,MEMLEN=13*2**(2*NG)/3+14*2**NG+8*NG-100/3)
C     PARAMETER (NG=5,MEMLEN=17*2**(2*NG)/3+18*2**NG+10*NG-86/3)
      INTEGER mem
      DOUBLE PRECISION z
      COMMON /memory/ z(MEMLEN),mem
      if (mem+len+1.gt.MEMLEN) pause 'insufficient memory in maloc'
      z(mem+1)=len
      maloc=mem+2
      mem=mem+len+1
      return
      END
