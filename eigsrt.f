      SUBROUTINE EIGSRT(D,V,N,NP)                                               
c given the eigenvalues D and eigenvectors V as output from jacobi.f this
c routine sorts the eigenvalues into ascending order, and rearranges the 
c columns of V correspondingly. The method is straight insertion.
c (c) numerical recipes
c sab 11.12.93 changed to sorting into descending order

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION D(NP),V(NP,NP)                                                  
      DO 13 I=1,N-1                                                             
        K=I                                                                     
        P=D(I)                                                                  
        DO 11 J=I+1,N                                                           
          IF(D(J).le.P)THEN                                                     
            K=J                                                                 
            P=D(J)                                                              
          ENDIF                                                                 
11      CONTINUE                                                                
        IF(K.NE.I)THEN                                                          
          D(K)=D(I)                                                             
          D(I)=P                                                                
          DO 12 J=1,N                                                           
            P=V(J,I)                                                            
            V(J,I)=V(J,K)                                                       
            V(J,K)=P                                                            
12        CONTINUE                                                              
        ENDIF                                                                   
13    CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
