      FUNCTION irbit1(iseed)
      INTEGER irbit1,iseed,IB1,IB2,IB5,IB18
      PARAMETER (IB1=1,IB2=2,IB5=16,IB18=131072)
      LOGICAL newbit
      newbit=iand(iseed,IB18).ne.0
      if(iand(iseed,IB5).ne.0)newbit=.not.newbit
      if(iand(iseed,IB2).ne.0)newbit=.not.newbit
      if(iand(iseed,IB1).ne.0)newbit=.not.newbit
      irbit1=0
      iseed=iand(ishft(iseed,1),not(IB1))
      if(newbit)then
        irbit1=1
        iseed=ior(iseed,IB1)
      endif
      return
      END
