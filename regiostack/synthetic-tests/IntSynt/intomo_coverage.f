      program intomocov
      include 'param_V3.h'
      real late,lone,lats,lons
      integer nbtra

      dimension late(NBTRAMAX), lone(NBTRAMAX)
      dimension lats(NBTRAMAX), lons(NBTRAMAX)


      call lecintomogeom(lats,lons,late,lone,nbtra)
      call write_coverage(lats,lons,late,lone,nbtra)

      end
