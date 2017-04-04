C ****************   C U B I C  ***********************************             
C                                                                               
      SUBROUTINE CUBIC(Z1,Z2,Z3,PER,IERR)                                       
C                                                                               
C     SOLVES FOR EIGENVALUES OF 3X3 SYMMETRIC MATRIX F                          
C     Z1.GE.Z2.GE.Z3 IN ORDER                                                   
C ** NOTE ** THE EIGENVALUES FROM EIGRS ARE ASCENDING ON OUTPUT                 
C                                                                               

      IMPLICIT REAL*8(A-H,O-Z)

      COMMON /MATX/ F(3,3),VEC(3,3)                                            
C                                                                               
      DIMENSION EIG(3),WORK(9)                                                  
C                                                                               
      IERR=0                                                                    
c sab 06.04.92
c sab 21.07.93
c keine IMSL mehr noetig
c Modifikation wg. Aenderungenen in der IMSL-Bibliothek
c      CALL EIGRS(F,3,2,EIG,VEC,3,WORK,IERR)                                    
c      IF (IERR .NE. 0) RETURN 
c      call evcsf(3,f,3,eig,vec,3)                                             
       call jacobi(f,3,3,eig,vec,ndumrot)
       call eigsrt(eig,vec,3,3)
c      PER=MAX(PER,WORK(1))                                                     
c      per=episf(3,3,f,3,eig,vec)  
      Z1=EIG(3)                                                                 
      Z2=EIG(2)                                                                 
      Z3=EIG(1)                                                                 
      RETURN                                                                    
      END                                                                       
C
