C                                                                               
C  Copyright 1990  Regents of the University of California                      
C                                                                               
C                                                                               
C  Author:  Dave Harris                                                         
C                                                                               
C           Lawrence Livermore National Laboratory                              
C           L-205                                                               
C           P.O. Box 808                                                        
C           Livermore, CA  94550                                                
C           USA                                                                 
C                                                                               
C           (415) 423-0617                                                      
C                                                    LPTBP                      
C                                                                               
C  Subroutine to convert an prototype lowpass filter to a bandpass filter via   
C    the analog polynomial transformation.  The lowpass filter is               
C    described in terms of its poles and zeros (as input to this routine).      
C    The output consists of the parameters for second order sections.           
C                                                                               
C  Input Arguments:                                                             
C  ----------------                                                             
C                                                                               
C    P                       Array containing poles                             
C                                                                               
C    Z                       Array containing zeros                             
C                                                                               
C    RTYPE                   Character array containing type information        
C                              (SP) single real pole  or                        
C                              (CP) complex conjugate pole pair  or             
C                              (CPZ) complex conjugate pole/zero pairs          
C                                                                               
C    DCVALUE                 Zero frequency value of filter                     
C                                                                               
C    NSECTS                  Number of second-order sections upon input         
C                                                                               
C    FL                      Low-frequency cutoff                               
C                                                                               
C    FH                      High-frequency cutoff                              
C                                                                               
C  Output Arguments:                                                            
C  -----------------                                                            
C                                                                               
C    SN                      Numerator polynomials for second order             
C                              sections.                                        
C                                                                               
C    SD                      Denominator polynomials for second order           
C                              sections.                                        
C                                                                               
C    NSECTS                  Number of second order sections upon output        
C                              This subroutine doubles the number of            
C                              sections.                                        
C                                                                               
C                                                                               
      SUBROUTINE LPTBP( P, Z, RTYPE, DCVALUE, NSECTS, FL, FH, SN, SD )  

        IMPLICIT NONE

        COMPLEX*16 P(*), Z(*)
        CHARACTER*3 RTYPE(*)                                            
        REAL*8 SN(*), SD(*), DCVALUE,FL,FH 
        INTEGER NSECTS

        COMPLEX*16  CTEMP, P1, P2, Z1, Z2, S, H
        REAL*8 PI,TWOPI,A,B,SCALE
        INTEGER N,I,IPTR

C                                                                               
        PI = 3.14159265                                                 
        TWOPI = 2.*PI                                                   
        A = TWOPI*TWOPI*FL*FH                                           
        B = TWOPI*( FH - FL )                                           

C                                                                               
        N = NSECTS                                                      
        NSECTS = 0                                                      
        IPTR = 1                                                        
        DO    1 I = 1, N                                                
C                                                                               
          IF (    RTYPE(I) .EQ. 'CPZ' ) THEN                            
C                                                                               
            CTEMP = ( B*Z(I) )**2 - 4.*A                                
            CTEMP = CDSQRT( CTEMP )                                      
            Z1 = 0.5*( B*Z(I) + CTEMP )                                 
            Z2 = 0.5*( B*Z(I) - CTEMP )                                 
            CTEMP = ( B*P(I) )**2 - 4.*A                                
            CTEMP = CDSQRT( CTEMP )                                      
            P1 = 0.5*( B*P(I) + CTEMP )                                 
            P2 = 0.5*( B*P(I) - CTEMP )                                 
            SN( IPTR )     = DBLE( Z1 * DCONJG( Z1 ) )                   
            SN( IPTR + 1 ) = -2. * DBLE( Z1 )                           
            SN( IPTR + 2 ) = 1.                                         
            SD( IPTR )     = DBLE( P1 * DCONJG( P1 ) )                   
            SD( IPTR + 1 ) = -2. * DBLE( P1 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
            SN( IPTR )     = DBLE( Z2 * DCONJG( Z2 ) )                   
            SN( IPTR + 1 ) = -2. * DBLE( Z2 )                           
            SN( IPTR + 2 ) = 1.                                         
            SD( IPTR )     = DBLE( P2 * DCONJG( P2 ) )                   
            SD( IPTR + 1 ) = -2. * DBLE( P2 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
C                                                                               
            NSECTS = NSECTS + 2                                         
C                                                                               
          ELSE IF (   RTYPE(I) .EQ. 'CP' ) THEN                         
C                                                                               
            CTEMP = ( B*P(I) )**2 - 4.*A                                
            CTEMP = CDSQRT( CTEMP )                                      
            P1 = 0.5*( B*P(I) + CTEMP )                                 
            P2 = 0.5*( B*P(I) - CTEMP )                                 
            SN( IPTR )     = 0.                                         
            SN( IPTR + 1 ) = B                                          
            SN( IPTR + 2 ) = 0.                                         
            SD( IPTR )     = DBLE( P1 * DCONJG( P1 ) )                   
            SD( IPTR + 1 ) = -2. * DBLE( P1 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
            SN( IPTR )     = 0.                                         
            SN( IPTR + 1 ) = B                                          
            SN( IPTR + 2 ) = 0.                                         
            SD( IPTR )     = DBLE( P2 * DCONJG( P2 ) )                   
            SD( IPTR + 1 ) = -2. * DBLE( P2 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             

C                                                                               
            NSECTS = NSECTS + 2                                         
C                                                                               
          ELSE IF (  RTYPE(I) .EQ. 'SP' ) THEN                          
C                                                                               
            SN( IPTR )     = 0.                                         
            SN( IPTR + 1 ) = B                                          
            SN( IPTR + 2 ) = 0.                                         
            SD( IPTR )     = A                                          
            SD( IPTR + 1 ) = -B*DBLE( P(I) )                            
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
C                                                                               
            NSECTS = NSECTS + 1                                         
C                                                                               
          END IF                                                        
C                                                                               
    1   CONTINUE                                                        
C                                                                               
C  Scaling - use the fact that the bandpass filter amplitude at sqrt( omega_l * 
C            equals the amplitude of the lowpass prototype at d.c.              
C                                                                               
        S = CMPLX( 0., DSQRT(A) )                                        
        H = CMPLX( 1., 0. )                                             
C                                                                               
        IPTR = 1                                                        
        DO    2 I = 1, NSECTS                                           
          H = H * ( ( SN(IPTR+2)*S + SN(IPTR+1) )*S + SN(IPTR) )        
     &      / ( ( SD(IPTR+2)*S + SD(IPTR+1) )*S + SD(IPTR) )            
          IPTR = IPTR + 3                                               
    2   CONTINUE                                                        
        SCALE = DCVALUE / SQRT( REAL( H ) * CONJG( H ) )                
                                                                        
        SN(1) = SN(1) * SCALE                                           
        SN(2) = SN(2) * SCALE                                           
        SN(3) = SN(3) * SCALE                                           
C                                                                               
      RETURN                                                            
      END                                                               
