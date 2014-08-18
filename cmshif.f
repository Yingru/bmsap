C******************   C M S H I F T   *********************************         
C                                                                               
      SUBROUTINE CMSHIF(PL1,PL2,PL3,XM,E,N,BX,BY,BZ,G,nmax)                          
C                                                                               
C                                                                               
C FIND TOTAL ENERGY-MOMENTUM 4-VECTOR AND LORENTZ TRANSFORM                     
C THE MOMENTA (IN-PLACE) TO THE CENTER-OF-MOMENTUM FRAME.                       
C                                                                               
c     DIMENSION PL1(10000),PL2(10000),PL3(10000),E(10000),XM(10000)

      IMPLICIT REAL*8(A-H,O-Z)

      dimension pl1(nmax),pl2(nmax),pl3(nmax)
      dimension e(nmax),xm(nmax)
C                                                                               
      PX=0.                                                                     
      PY=0.                                                                     
      PZ=0.                                                                     
      ET=0.                                                                     
C                                                                               
C ...  COMPUTE TOTAL 4-MOMENTUM                                                 
C                                                                               
      DO 1000 I=1,N                                                             
        PX=PX+PL1(I)                                                            
        PY=PY+PL2(I)                                                            
        PZ=PZ+PL3(I)                                                            
        ET=ET+E(I)                                                              
 1000  CONTINUE                                                                 
C                                                                               
C ...  FIND BETA-C.M. AND GAMMA                                                 
C                                                                               
      BX=PX/ET                                                                  
      BY=PY/ET                                                                  
      BZ=PZ/ET                                                                  
      G=1./SQR(1.-BX*BX-BY*BY-BZ*BZ)                                            
      GG=1./(1.+G)                                                              
      GBX=G*BX                                                                  
      GBY=G*BY                                                                  
      GBZ=G*BZ                                                                  
C                                                                               
C ...  TRANSFORM                                                                
C                                                                               
      DO 2000 I=1,N                                                             
        PX=PL1(I)                                                               
        PY=PL2(I)                                                               
        PZ=PL3(I)                                                               
        TEMP=(PX*GBX+PY*GBY+PZ*GBZ)*GG-E(I)                                     
        PL1(I)=PX+GBX*TEMP                                                      
        PL2(I)=PY+GBY*TEMP                                                      
        PL3(I)=PZ+GBZ*TEMP                                                      
        E(I)=E(I)/G                                                             
 2000  CONTINUE                                                                 
C                                                                               
      RETURN                                                                    
      END                                                                       
C




