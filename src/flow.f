C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
C                                                                               
        SUBROUTINE FLOW(PL1,PL2,PL3,XM,LCOUNT,PER,IERR,EK,VAR,II,
     &                  nmax,rapcm)               
c        SUBROUTINE FLOW(PL1,PL2,PL3,XM,LCOUNT,PER,IERR,EK,II)
C                                                                               
C       PL1 IS CM MOMENTUM ALONG X DIRECTION                                    
C       PL2     ...              Y    ...                                       
C       PL3     ...              Z    ...                                       
C       XM  IS  MASS OF FRAGMENT                                                
C       LCOUNT IS TOTAL NUMBER OF FRAGMENTS
c       ek is E_tot
c       ii is the groupcounter
c       per is useless
c       ierr is errorcode for eigenvalue calculation

      IMPLICIT REAL*8(A-H,O-Z) 

c	include 'stndrdb.f'
        dimension pl1(nmax),pl2(nmax),pl3(nmax),xm(nmax)
        dimension ek(nmax)
c      DIMENSION PL1(10000),PL2(10000),PL3(10000),XM(10000),EK(10000)                 
        COMMON /MATX/ F(3,3),VEC(3,3)                                          
        COMMON /FLSTAB/ QAVSQ(4),QAVSQM(4)                                      
        DIMENSION VAR(60,30)                                                     
        DIMENSION U(3),V(3),W(3),G(3,3),E(3)                                    
C                                                                               
      save

C                                                                               
      PER=0.                                                                    
      IF (LCOUNT .GT. 0) THEN                                                   
        CALL CMSHIF(PL1,PL2,PL3,XM,EK,LCOUNT,BX,BY,BZ,GAM,nmax)                      
        VAR(1,II)=BX                                                            
        VAR(2,II)=BY                                                            
        VAR(3,II)=BZ                                                            
        VAR(46,II)=GAM                                                          
      ELSE                                                                      
c        WRITE(30,*)0.,0.,0.                                                     
        VAR(1,II)=0.                                                            
        VAR(2,II)=0.                                                            
        VAR(3,II)=0.                                                            
        VAR(46,II)=0.0                                                          
      END IF                                                                    
C                                                                               
C     IF (LCOUNT .EQ. 0) THEN                                                   
C       WRITE(30,*)(0.0,I=1,6)                                                  
C       DO 500 I=1,3                                                            
C         WRITE(30,*)(0.0,J=1,9)                                                
C500    CONTINUE                                                                
C       RETURN                                                                  
C     END IF                                                                    
C                                                                               
        RAD=57.29578                                                            
        IA=ABS(LCOUNT)                                                          
        IiCOUNT=0                                                                
        XIA=FLOAT(IA)                                                           
        VAR(60,II)=XIA                                                          
        P1DAV=0.                                                                
        P1AV=0.                                                                 
        P2AV=0.                                                                 
        P3AV=0.                                                                 
        PAV=0.                                                                  
        PMAGAV=0.                                                               
        P1AVSQ=0.                                                               
        P2AVSQ=0.                                                               
        P3AVSQ=0.                                                               
        PAVSQ=0.                                                                
        ARA3=0.                                                                 
        EKAV=0.                                                                 
        EKAVSQ=0.                                                               
        BXAV=0.                                                                 
        BXAVSQ=0.                                                               
        BYAV=0.                                                                 
        BYAVSQ=0.                                                               
        BZAV=0.                                                                 
        BZAVSQ=0.                                                               
        QXY=0.                                                                  
        QXZ=0.                                                                  
        QYZ=0.                                                                  
        FXX=0                                                                   
        FYY=0                                                                   
        FZZ=0                                                                   
        FXY=0.                                                                  
        FXZ=0.                                                                  
        FYZ=0.                                                                  
        RQXX=0.                                                                 
        RQYY=0.                                                                 
        RQZZ=0.                                                                 
        RQXY=0.                                                                 
        RQXZ=0.                                                                 
        RQYZ=0.                                                                 
C                                                                               
        PXYT=0.0                                                                
        PTSQ=0.0                                                                
        PTAV=0.0                                                                
        QZZAV=0.0
        direct=0.0
        dirpt=0.0
        dirpx=0.0
        dirpy=0.0
C                                                                               
C                                                                               
        DO 200 I=1,IA
           etot=sqrt(xm(i)**2+pl1(i)**2+pl2(i)**2+pl3(i)**2)
           rap=0.5*dlog( (etot+pl3(i)) / (etot-pl3(i)) )
           if(rap.gt.rapcm) then
              dirpt=dirpt+SQR(PL1(I)**2+PL2(I)**2)
              dirpx=dirpx+pl1(i)
              dirpy=dirpy+pl2(i)
           endif
        iICOUNT=iICOUNT+1                                                         
        PMAG=SQR(PL1(I)**2+PL2(I)**2+PL3(I)**2)                                 
        PTRA=SQR(PL1(I)**2+PL2(I)**2)                                           
        P1AV=P1AV+PL1(I)
c        write(6,*) "check1"                                            
        P1DAV=P1DAV+PL1(I)*SIGN(1.,real(PL3(I)))  
c        write(6,*) "check2"                    
        P2AV=P2AV+PL2(I)                                                        
        P3AV=P3AV+PL3(I)                                                        
        PMAGAV=PMAGAV+PMAG                                                      
        P1AVSQ=P1AVSQ+PL1(I)**2                                                 
        P2AVSQ=P2AVSQ+PL2(I)**2                                                 
        P3AVSQ=P3AVSQ+PL3(I)**2                                                 
        EKAV=EKAV+EK(I)                                                         
        EKAVSQ=EKAVSQ+EK(I)**2                                                  
        BXAV=BXAV+PL1(I)/EK(I)                                                  
        BXAVSQ=BXAVSQ+(PL1(I)/EK(I))**2                                         
        BYAV=BYAV+PL2(I)/EK(I)                                                  
        BYAVSQ=BYAVSQ+(PL2(I)/EK(I))**2                                         
        BZAV=BZAV+PL3(I)/EK(I)                                                  
        BZAVSQ=BZAVSQ+(PL3(I)/EK(I))**2                                         
        PTAV=PTAV+PTRA                                                          
        PTSQ=PTSQ+PTRA*PTRA                                                     
        PXYT=PXYT+PL1(I)*PL1(I)-PL2(I)*PL2(I)                                   
        QZZAV=QZZAV+2.*PL3(I)*PL3(I)-PTRA*PTRA                                  
        FXX=FXX+PL1(I)**2/PMAG                                                  
        FYY=FYY+PL2(I)**2/PMAG                                                  
        FZZ=FZZ+PL3(I)**2/PMAG                                                  
        QXY=QXY+PL1(I)*PL2(I)                                                   
        QXZ=QXZ+PL1(I)*PL3(I)                                                   
        QYZ=QYZ+PL2(I)*PL3(I)                                                   
        FXY=FXY+PL1(I)*PL2(I)/PMAG                                              
        FXZ=FXZ+PL1(I)*PL3(I)/PMAG                                              
        FYZ=FYZ+PL2(I)*PL3(I)/PMAG                                              
        RQXX=RQXX+PL1(I)**2/XM(I)                                               
        RQYY=RQYY+PL2(I)**2/XM(I)                                               
        RQZZ=RQZZ+PL3(I)**2/XM(I)                                               
        RQXY=RQXY+PL1(I)*PL2(I)/XM(I)                                           
        RQXZ=RQXZ+PL1(I)*PL3(I)/XM(I)                                           
        RQYZ=RQYZ+PL2(I)*PL3(I)/XM(I)                                           
200     CONTINUE                                                                
        PAVSQ=P1AVSQ+P2AVSQ+P3AVSQ                                              
C                                                                               
C                                                                               
        dirpx=dirpx/xia
        dirpy=dirpy/xia
        dirpt=dirpt/xia
        direct= sqr(dirpx**2+dirpy**2)/max(0.0000001,dirpt)
        VAR(4,II)=P1AV/XIA                                                      
        VAR(47,II)=P1DAV/XIA                                                    
        VAR(5,II)=P2AV/XIA                                                      
        VAR(6,II)=P3AV/XIA                                                      
        PAV=PAV/XIA                                                             
        VAR(49,II)=PMAGAV/XIA                                                   
        VAR(51,II)=EKAV/XIA                                                     
        VAR(53,II)=BXAV/XIA                                                     
        VAR(54,II)=BYAV/XIA                                                     
        VAR(55,II)=BZAV/XIA                                                     
        QXX=P1AVSQ                                                              
        QYY=P2AVSQ                                                              
        QZZ=P3AVSQ                                                              
        RQXX=RQXX/(2.*XIA)                                                      
        RQYY=RQYY/(2.*XIA)                                                      
        RQZZ=RQZZ/(2.*XIA)                                                      
        RQXY=RQXY/(2.*XIA)                                                      
        RQXZ=RQXZ/(2.*XIA)                                                      
        RQYZ=RQYZ/(2.*XIA)                                                      
C                                                                               
        P1AVSQ=P1AVSQ/XIA                                                       
        P2AVSQ=P2AVSQ/XIA                                                       
        P3AVSQ=P3AVSQ/XIA                                                       
        PAVSQ=PAVSQ/XIA                                                         
        EKAVSQ=EKAVSQ/XIA                                                       
        BXAVSQ=BXAVSQ/XIA                                                       
        BYAVSQ=BYAVSQ/XIA                                                       
        BZAVSQ=BZAVSQ/XIA                                                       
        ARA3=P3AVSQ/PAVSQ                                                       
        QAVSQ(II)=2.*P3AVSQ-P1AVSQ-P2AVSQ                                       
        PTAV=PTAV/XIA                                                           
        PXYT=PXYT/PTSQ                                                          
        QZZAV=QZZAV/XIA/PMAGAV                                                     
        PTSQ=SQR(PTSQ/XIA-PTAV*PTAV)                                            
C                                                                               
C                                                                               
        VAR(48,II) = SQR(P1AVSQ-VAR(47,II)**2)                                  
        VAR(7,II)  = SQR(P1AVSQ-VAR(4,II)**2)                                   
        VAR(8,II)  = SQR(P2AVSQ-VAR(5,II)**2)                                   
        VAR(9,II)  = SQR(P3AVSQ-VAR(6,II)**2)                                   
        VAR(50,II) = SQR(PAVSQ-VAR(49,II)**2)                                   
        VAR(52,II) = SQR(EKAVSQ-VAR(51,II)**2)                                  
        VAR(56,II) = SQR(BXAVSQ-VAR(53,II)**2)                                  
        VAR(57,II) = SQR(BYAVSQ-VAR(54,II)**2)                                  
        VAR(58,II) = SQR(BZAVSQ-VAR(55,II)**2)                                  
C       VAR(59,II) = BXAVSQ+BYAVSQ+BZAVSQ-VAR(53,II)*VAR(53,II)                 
C    *              -VAR(54,II)*VAR(54,II)-VAR(55,II)*VAR(55,II)                
        SIGP  = 0.                                                              
C                                                                               
C       NOW SOME CHANGED DEFINITIONS OF VAR(53,54,55,59)                        
C             PREV. MEAN BETAX BETAY BETAZ  SIG2 BETA  NOW:                     
        VAR(53,II)=PTAV                                                         
        VAR(54,II)=PTSQ                                                         
        VAR(55,II)=PXYT                                                         
        VAR(59,II)=QZZAV                                                        
        var(46,ii)=direct
C                                                                               
c      WRITE(30,*)(VAR(I,II),I=4,9),VAR(47,II),VAR(48,II)                        
C                                                                               
C                                                                               
C       NOW DO 3D MOMENTUM FLOW ANALYSIS                                        
C                                                                               
        F(1,1)=FXX                                                              
        F(2,2)=FYY                                                              
        F(3,3)=FZZ                                                              
        F(2,1)=FXY
        f(1,2)=fxy                                                           
        F(1,3)=FXZ
        f(3,1)=fxz                                                              
        F(2,3)=FYZ
        f(3,2)=fyz                                                              
        CALL CUBIC(F1,F2,F3,PER,IERR)                                           
        IF (IERR .NE. 0) RETURN                                                 
        CALL ANGLE(1,TH1,PH1)                                                   
        CALL ANGLE(2,TH2,PH2)                                                   
        CALL ANGLE(3,TH3,PH3)                                                   
C                                                                               
C     NOW DO 3D SPHERICITY FLOW                                                 
C                                                                               
       FACTR=1.                                                                 
       F(1,1)=QXX*FACTR
       F(2,2)=QYY*FACTR
       F(3,3)=QZZ*FACTR
       F(1,2)=QXY*FACTR
       F(2,1)=QXY*FACTR                                                         
       F(1,3)=QXZ*FACTR            
       F(3,1)=QXZ*FACTR            
       F(2,3)=QYZ*FACTR            
       F(3,2)=QYZ*FACTR            
       CALL CUBIC(Q1,Q2,Q3,PER,IERR)                                            
       IF (IERR .NE. 0) RETURN                                                  
       CALL ANGLE(1,QT1,QP1)                                                    
       CALL ANGLE(2,QT2,QP2)                                                    
       CALL ANGLE(3,QT3,QP3)                                                    
C                                                                               
C     NOW DO 3D KINETIC FLOW ANALYSIS                                           
C                                                                               
       F(1,1)=RQXX                                                              
       F(2,2)=RQYY                                                              
       F(3,3)=RQZZ                                                              
       F(1,2)=RQXY                                                             
       F(2,1)=RQXY                                                              
       F(1,3)=RQXZ                                                              
       F(3,1)=RQXZ                                                              
       F(2,3)=RQYZ                                                              
       F(3,2)=RQYZ                                                              
       CALL CUBIC(RQ1,RQ2,RQ3,PER,IERR)                                         
       IF (IERR .NE. 0) RETURN                                                  
       CALL ANGLE(1,RQT1,RQP1)                                                  
       CALL ANGLE(2,RQT2,RQP2)                                                  
       CALL ANGLE(3,RQT3,RQP3)                                                  
c       WRITE (30,*) F1,TH1,PH1,F2,TH2,PH2,F3,TH3,PH3                            
c       WRITE (30,*) Q1,QT1,QP1,Q2,QT2,QP2,Q3,QT3,QP3                            
c       WRITE (30,*) RQ1,RQT1,RQP1,RQ2,RQT2,RQP2,RQ3,RQT3,RQP3                   
C                                                                               
       VAR(10,II)= F1                                                           
       VAR(11,II)= TH1                                                          
       VAR(12,II)= PH1                                                          
       VAR(13,II)= F2                                                           
       VAR(14,II)= TH2                                                          
       VAR(15,II)= PH2                                                          
       VAR(16,II)= F3                                                           
       VAR(17,II)= tH3                                                          
       VAR(18,II)= pH3                                                          
       VAR(21,II)= Q1                                                           
       VAR(22,II)= QT1                                                          
       VAR(23,II)= QP1                                                          
       VAR(24,II)= Q2                                                           
       VAR(25,II)= QT2                                                          
       VAR(26,II)= QP2                                                          
       VAR(27,II)= Q3                                                           
       VAR(28,II)= QT3                                                          
       VAR(29,II)= QP3                                                          
       VAR(32,II)= RQ1                                                          
       VAR(33,II)= RQT1                                                         
       VAR(34,II)= RQP1                                                         
       VAR(35,II)= RQ2                                                          
       VAR(36,II)= RQT2                                                         
       VAR(37,II)= RQP2                                                         
       VAR(38,II)= RQ3                                                          
       VAR(39,II)= RQT3                                                         
       VAR(40,II)= RQP3                                                         
      CALL RATIO(VAR,19,10,16,II)                                               
      CALL RATIO(VAR,20,13,16,II)                                               
      CALL RATIO(VAR,30,21,27,II)                                               
      CALL RATIO(VAR,31,24,27,II)                                               
      CALL RATIO(VAR,41,32,38,II)                                               
      CALL RATIO(VAR,42,35,38,II)                                               
C                                                                               
      VAR(43,II)=COS(0.0174533*VAR(11,II))                                      
      VAR(44,II)=COS(0.0174533*VAR(22,II))                                      
      VAR(45,II)=COS(0.0174533*VAR(33,II))                                      
C                                                                               
       RETURN                                                                   
       END                                                                      
C
