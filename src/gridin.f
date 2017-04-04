      subroutine gridin
c sab initialize grids and book observables into them

      IMPLICIT REAL*8(A-H,O-Z)

      include 'ucoms.f'
      include 'hicoms.f'

c sab calculate grid-intervalls
c 
      do 111 i=1,zaxis
         xbins(i)=max(1,min(nbgr,xbins(i)))
         if(xbins(i).gt.1) then
            dxbin(i)=(xmax(i)-xmin(i))/(xbins(i)-1)
         else
            dxbin(i)=(xmax(i)-xmin(i))
            xmin(i)=xmin(i)+dxbin(i)*0.5d0
         endif
         if(ybins(i).gt.1) then
            dybin(i)=(ymax(i)-ymin(i))/(ybins(i)-1)
         else
            dybin(i)=(ymax(i)-ymin(i))
         endif
 111  continue
c           start calculation with grid                                        
c sab i1grdm,i2grdm muessen noch aus bbb herausgeklaut werden
        do 6351 i=1,30                                                      
           do 6352 j=1,nbgr
              do 222 k=1,nbgr
                zxygrd(i,j,k)=0.0d0                                             
                dnsgrd(i,j,k)=0
                dnsgrd(i,j,k)=0.0d0
                errgrd(i,j,k)=0.0d0
 222        continue
 6352      continue                                                           
 6351   continue                                                           
         return
         end
