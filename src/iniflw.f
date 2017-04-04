C     ************ I N I F L W ****************************************
C                                      
      SUBROUTINE INIFLW                

      IMPLICIT REAL*8(A-H,O-Z)

      include 'ucoms.f'
      include 'hicoms.f'
      DIMENSION E(NMAX),XM(NMAX),ekin(NMAX),
     *          PL1(NMAX),PL2(NMAX),PL3(NMAX)  
C
c hier rapcm berechnen, wird fuer directivity benoetigt:
C       RAPIDITIES OF PROJECTILE AND TARGET
      EPROJ=SQRT(.9380**2+PPROJ**2)
      ETARG=SQRT(.9380**2+PTARG**2)
      RAPPR=0.5*DLOG((EPROJ+PPROJ)/(EPROJ-PPROJ))
      RAPTA=0.5*DLOG((ETARG+PTARG)/(ETARG-PTARG))
      XPROJ=Ap
      XTARG=At
      RAPCM=(XPROJ*RAPPR + XTARG*RAPTA)/max(1.,(XPROJ+XTARG))


c loop ueber alle gruppen
      do 10 igrp=1,grpsel
         ilc=0
c ab hier do-loop ueber alle Teilchen in gruppe  
            if(ap.eq.at) then                                           
c here symmetric momentum flow
            do 1301 ip=1,nopart(igrp)
               i=grpart(igrp,ip)
               Ekin(i)=abs(SQRT(FMASS(i)*FMASS(i)+PX(i)*PX(i)+PY(i)*
     &                 PY(i)+PZ(i)*PZ(i))-fmass(i))
               PL1(IP)=PX(i)                               
               PL2(IP)=PY(i)  
               PL3(IP)=PZ(i)                               
               XM(IP)=fmass(i) 
               E(IP)=EKIN(i)
               if(refsys.eq.1) then
                  e1=ekin(i)+fmass(i)
                  pl3(ip)=pl3(ip)*gammal-betal*gammal*e1
               end if
 1301       continue
            ipart=nopart(igrp)
            else
c here asymmetric momentum flow
               do 1302 ip=1,nopart(igrp)
               ik=grpart(igrp,ip)
               Ekin(ik)=abs(SQRT(FMASS(ik)*FMASS(ik)+
     &                 PX(ik)*PX(ik)+PY(ik)*
     &                 PY(ik)+PZ(ik)*PZ(ik))-fmass(ik))
               if(pz(ik).gt.0) then
                  ilc=ilc+2
                  PL1(ilc-1)=-PX(ik)                                         
                  PL2(ilc-1)=-PY(ik)                                         
                  PL3(ilc-1)=-PZ(ik)                                         
                  XM(ilc-1)=fmass(ik)                                        
                  E(ilc-1)=EKIN(ik)                                          
                  PL1(ilc)=PX(ik)                                            
                  PL2(ilc)=PY(ik)                                            
                  PL3(ilc)=PZ(ik)                                            
                  XM(ilc)=fmass(ik)                                          
                  E(ilc)=EKIN(ik)
                  if(refsys.eq.1) then
                     e1=ekin(ik)+fmass(ik)
                     pl3(ilc)=pl3(ilc)*gammal-betal*gammal*e1
                     pl3(ilc-1)=pl3(ilc-1)*gammal-betal*gammal*e1
                  end if
               endif
 1302          continue
               ipart=ilc               
            endif
c lcount bestimmen : anzahl der Teilchen in gruppe
         IF(ipart.GT.3)THEN
            CALL FLOW(PL1,PL2,PL3,XM,ipart,PER,IRET,E,VAR,
     &                igrp,nmax,rapcm)              
         ELSE                                                              
            DO 1310 J=1,60                                                  
               VAR(J,igrp)=0.0                                                 
 1310       CONTINUE
         endif
 10   continue
c ende                                                            
      RETURN                                                                    
      END




















