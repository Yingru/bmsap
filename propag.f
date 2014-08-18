      subroutine propag(ntfr1,ntfr2,x,y,z,px,py,pz,xnew,ynew,znew)

      IMPLICIT REAL*8(A-H,O-Z)

c  von Christiane Lesny
c  17.11.1994

      parameter (m=938.)
   
c Differenz der Zeiten
      ndeltat = ntfr1-ntfr2

c Propagation der Orte      
       xnew = x + px/m*float(ndeltat)
       ynew = y + py/m*float(ndeltat)
       znew = z + pz/m*float(ndeltat)

      return
      end
