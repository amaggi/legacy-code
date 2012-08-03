      FUNCTION irbit2(iseed)
      INTEGER irbit2,iseed,IB1,IB2,IB5,IB18,MASK
      PARAMETER (IB1=1,IB2=2,IB5=16,IB18=131072,MASK=IB1+IB2+IB5)
      if(iand(iseed,IB18).ne.0)then
        iseed=ior(ishft(ieor(iseed,MASK),1),IB1)
        irbit2=1
      else
        iseed=iand(ishft(iseed,1),not(IB1))
        irbit2=0
      endif
      return
      END
