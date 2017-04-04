      subroutine frag
c     this subroutine determines fragments according to a coalescence model
c     
c     author  : Steffen A. Bass
c     date    : 21.02.94
c     revision: 0.1 beta
c     
c     based on Clusters, (c) by Jens Konopka /Univ. Frankfurt
c

      IMPLICIT REAL*8(A-H,O-Z)
      include 'ucoms.f'
      include 'hicoms.f'

      dimension xn(nmax),yn(nmax),zn(nmax),                      
     &          pxn(nmax),pyn(nmax),pzn(nmax),  
     &          emn(nmax),itypn(nmax),ichargen(nmax),icharge(nmax)
     &          lastcln(nmax),nclcntn(nmax),
     &          iclarr(nmax,nmax),number(nmax)
      np=nnucl
c     itotpan is the new all-particle counter
      itotpan=0
c initalization, generate charge-vector
c
c     iparts(k) will later give the number of nucl. in a cluster -> mass
c     iclarr gives the linkes between particles
         do 80 k=1,np
            iparts(k) = 1
            icharge(k)=charge(k)
            do 80 l=1,np
               iclarr(k,l) = 0
  80     continue

c ***
c *** distance criterion
c ***

         do 400 k=1,np
            do 410 l=k+1,np
               dr2 = (rx(l)-rx(k))**2 + (ry(l)-ry(k))**2 + 
     +              (rz(l)-rz(k))**2
c               dp2 = (px(l)-px(k))**2 + (py(l)-py(k))**2 +
c     +               (pz(l)-pz(k))**2
               if ( dr2 .le. 9.00 ) then

c ***
c *** particle k and l belong to the same cluster
c ***

                  iclarr(k,l) = 1
                  iclarr(l,k) = 1
               end if
  410       continue
  400    continue


c ***
c *** counting of clusters number, size and charge
c ***
         do 120 k=1,np
            if ( iparts(k) .lt. 1 ) then
               goto 120
            end if
  777       continue
            iwh = 0
            do 130 l=k+1,np
               if ( iclarr(k,l) .eq. 1 ) then
                  iclarr(l,k) = 0
                  iparts(l) = 0
                  do 140 m=1,np
                     if ( iclarr(l,m) .eq. 1 ) then
                        iparts(m) = 0
                        iclarr(k,m) = 1
                        iclarr(m,k) = 0
                        iclarr(l,m) = 0
                        iclarr(m,l) = 0
                        if ( m .le. l ) iwh = 1
                     end if
  140             continue
               end if
  130       continue
            if ( iwh .eq. 1 ) goto 777
            do 144 l=k+1,np
               if ( iclarr(k,l) .eq. 1 ) then
                  iparts(k) = iparts(k)+1
                  icharge(k) = icharge(k)+icharge(l)
               end if
  144       continue
  120    continue

c final loop: cluster forming...
         nsingle=0
         nfrag=0
         do 150 k=1,np
            if ( iparts(k) .eq. 0 ) then
               goto 150
            else if ( iparts(k) .eq. 1 ) then
               itotpan=itotpan+1
               nsingle=nsingle+1
               xn(itotpan)=rx(k)
               yn(itotpan)=ry(k)
               zn(itotpan)=rz(k)
               pxn(itotpan)=px(k)
               pyn(itotpan)=py(k)
               pzn(itotpan)=pz(k)
               emn(itotpan)=fmass(k)
               itypn(itotpan)=ityp(k)
               ichargen(itotpan)=icharge(k)
               lastcln(itotpan)=lstcoll(k)
               nclcntn(itotpan)=ncoll(k)
               goto 150
            endif
c
c *** Numerierungsprozedur fuer ein Cluster
            itotpan=itotpan+1
            nfrag=nfrag+1
            index = 1
            number(index) = k
            do 160 l=k+1,np
               if ( iclarr(k,l) .eq. 1 ) then
                  index = index+1
                  number(index) = l
               end if
  160       continue

c ***
c *** cluster center of mass coordinates and momentum
c ***
c hier muss spaeter mal eine Abfrage rein, ob deuteronen
c in frag oder deut berechnet werden.
            xcm = 0.
            ycm = 0.
            zcm = 0.
            pxcm = 0.
            pycm = 0.
            pzcm = 0.
            totem = 0.
            ntotcl = 0
            rhom = 0.0
            frti=0
            mfdelt=0.0
            npimi=0
            npize=0
            npipl=0
            do 170 k1=1,iparts(k)
               j=number(k1)

               xcm = xcm + rx(j)*fmass(j)
               ycm = ycm + ry(j)*fmass(j)
               zcm = zcm + rz(j)*fmass(j)
               pxcm = pxcm + px(j)
               pycm = pycm + py(j)
               pzcm = pzcm + pz(j)
               totem = totem + fmass(j)

               ntotcl = ntotcl + ncoll(j)

  170       continue

c            div = real(iparts(k))
            xcm = xcm/totem
            ycm = ycm/totem
            zcm = zcm/totem
c these are momenta per nucleon
c            pxcm = pxcm/div
c            pycm = pycm/div
c            pzcm = pzcm/div

            xn(itotpan)=xcm
            yn(itotpan)=ycm
            zn(itotpan)=zcm
            pxn(itotpan)=pxcm
            pyn(itotpan)=pycm
            pzn(itotpan)=pzcm
            emn(itotpan)=totem
            itypn(itotpan)=-1
            ichargen(itotpan)=icharge(k)
            lastcln(itotpan)=-99
            nclcntn(itotpan)=ntotcl

            
c            do 180 k1=1,iparts(k)
c               j=number(k1)
c               x(j) = x(j) - xcm   
c               y(j) = y(j) - ycm   
c               z(j) = z(j) - zcm   
c               px(j) = px(j) - pxcm
c               py(j) = py(j) - pycm
c               pz(j) = pz(j) - pzcm
c  180       continue
c hier werden Bindungsenergien berechnet - muss total umgestellt werden

c ***** end of main-loop
  150    continue

c insert pions into new array
         do 200 k=(np+1),npart
            itotpan=itotpan+1
               xn(itotpan)=rx(k)
               yn(itotpan)=ry(k)
               zn(itotpan)=rz(k)
               pxn(itotpan)=px(k)
               pyn(itotpan)=py(k)
               pzn(itotpan)=pz(k)
               emn(itotpan)=fmass(k)
               itypn(itotpan)=ityp(k)
               ichargen(itotpan)=charge(k)
               lastcln(itotpan)=lstcoll(k)
               nclcntn(itotpan)=ncoll(k)
 200        continue
            
c rewrite particle-arrays
            do 300 k=1,itotpan
               rx(k)=xn(k)
               ry(k)=yn(k)
               rz(k)=zn(k)
               px(k)=pxn(k)
               py(k)=pyn(k)
               pz(k)=pzn(k)
               fmass(k)=emn(k)
               ityp(k)=itypn(k)
               charge(k)=ichargen(k)
               lstcoll(k)=lastcln(k)
               ncoll(k)=nclcntn(k)
 300        continue
            npart=itotpan
      return
      end



