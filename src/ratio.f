C **************  R A T I O  *******************                                
C                                                                               
      SUBROUTINE RATIO(V,I,J,K,II)                                              

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION V(60,30)                                                         
      V(I,II)=-999.                                                             
      IF (V(K,II) .NE. 0.0) V(I,II)=V(J,II)/V(K,II)                             
      RETURN                                                                    
      END                                                                       
C
