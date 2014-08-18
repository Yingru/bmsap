C **************   A N G L E  *********************************                 
C                                                                               
       SUBROUTINE ANGLE(IEIG,THET,PHI)                                          
C                                                                               
C      EVALUATES POLAR AND AZIMUTHAL ANGLE                                      
C      OF FLOW EIGENVECTOR WITH LENGTH XL                                       
C ** NOTE ** THE IEIG INDEX REFERS TO THE DECENDING NUMBERING                   
C            SCHEME, WHILE THE VECTORS IN VEC ARE STORED USING                  
C            ASCENDING EIGENVALUE ORDERING.                                     
C-------------------------------------------------------------------            
C 6/9/83 EIGENVECTORS ARE NORMALISED WITH A POSITIVE Z COMPONENT                
C                                                                               

       IMPLICIT REAL*8(A-H,O-Z)

       COMMON /MATX/ F(3,3),VEC(3,3)                                              
C      THIS IS FLOW MATRIX                                                      
       DATA RAD/57.29578/,EPS/1.E-15/                                           
       J=4-IEIG                                                                 
C                                                                               
       IF (VEC(3,J) .LT. 0.) THEN                                               
         DO 500 I=1,3                                                           
           VEC(I,J)=-VEC(I,J)                                                   
 500     CONTINUE                                                               
       END IF                                                                   
C                                                                               
       THET=RAD*ACOS(VEC(3,J))                                                  
       X=VEC(1,J)                                                               
       Y=VEC(2,J)                                                               
       IF ( (ABS(X) .LT. EPS) .AND. (ABS(Y) .LT. EPS)) THEN                     
         PHI=0.                                                                 
       ELSE                                                                     
         PHI=RAD*ATAN2(Y,X)                                                     
       END IF                                                                   
       RETURN                                                                   
       END                                                                      
C
