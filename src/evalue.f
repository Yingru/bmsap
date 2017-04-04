      function evalue(obs,jgrp)
c     this function returns the value of observable obs of event consisting
c     of particles in group jgrp
c     
c     author  : Steffen A. Bass
c     date    : 17.02.94
c     revision: 1.1
c
      IMPLICIT REAL*8(A-H,O-Z)

      real*8 evalue
      integer obs,jgrp,i,j
      parameter(degrad=180./3.1415927)

c
      include 'ucoms.f'
      include 'hicoms.f'
c
c WARNING:
c for the calculation of eventlike values only the particles
c in the respective group are used, which have passed all
c particle-like cuts referring to them!!!
c this means a particle-loop from 1 to nopart(jgrp)
c and accessing the particles pointed to in
c grpart(jgrp,nopart(jgrp))
c
c now goto the desired observable
c           1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
      goto( 1, 2, 3, 4, 5, 6, 7,99,99,99,11,11,11,11,11,11,11,11,11,11,
     &     11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,
     &     11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,
     &     11,11,11,11,11,11,11,11,11,11,99,99,99,99,99,99,99,99,99,99,
     &     99,82,83,84,85,86,87,88,89,99,99,99,99,99,99,99,99,99,99,
     &     99,99), (obs-100)
      goto 99
c
c here come the observables
c feel free to add your favourite observable;
c in case of changes, please contact the authors 
c for a general revision update.
c
c dN/d(obsvalue) distribution
 1    evalue=1.
      return
c impact parameter
 2    evalue=bimp
      return
c particle multiplicity (# of particles in group)
 3    evalue=nopart(jgrp)
      return
c incident beam energy (lab frame)
 4    evalue=ebeam
      return
c center of mass energy (sqrt(s))
 5    evalue=ecm
      return
c momentum of beam particles (lab-frame)
 6    evalue=pbeam
      return
c momentum of beam particles (cm-frame)
 7    evalue=pcm
      return
c
 11   evalue=var(obs-110,jgrp)
      return
c number of hard collisions
 82   continue
      evalue=0d0
      do 111 i=1,inpart
         evalue=evalue+dble(NINT(lstcoll(I))/MSTV3**2)
 111  continue
      evalue=evalue/2d0
      return
c <rt> 
 83   continue
      if(nopart(jgrp).gt.0) then
      do 112 i=1,nopart(jgrp)
         evalue=evalue+
     &        sqrt(frrx(grpart(jgrp,i))**2+frry(grpart(jgrp,i))**2)
 112  continue
      evalue=evalue/nopart(jgrp)
      else
         evalue=0d0
      endif
      return
c <tf>
 84   continue
      if(nopart(jgrp).gt.0) then
      do 113 i=1,nopart(jgrp)
         evalue=evalue+frr0(grpart(jgrp,i))
 113  continue
      evalue=evalue/nopart(jgrp)
      else
         evalue=0d0
      endif
      return
c <pt>
 85   continue
      if(nopart(jgrp).gt.0) then
      do 114 i=1,nopart(jgrp)
         evalue=evalue+sqrt(px(grpart(jgrp,i))**2+py(grpart(jgrp,i))**2)
 114  continue
      evalue=evalue/nopart(jgrp)
      else
         evalue=0d0
      endif
      return
c E_t of the particle group in the event
 86   continue
      if(nopart(jgrp).gt.0) then
      do 115 i=1,nopart(jgrp)
         evalue=evalue+sqrt(px(grpart(jgrp,i))**2+py(grpart(jgrp,i))**2
     &                      +fmass(grpart(jgrp,i))**2)
 115  continue
c      evalue=evalue/nopart(jgrp)
      else
         evalue=0d0
      endif
      return
c E_tot of the particle group in the event
 87   continue
      if(nopart(jgrp).gt.0) then
         do 116 i=1,nopart(jgrp)
         evalue=evalue+sqrt(px(grpart(jgrp,i))**2+py(grpart(jgrp,i))**2
     &                      +pz(grpart(jgrp,i))**2
     &        +fmass(grpart(jgrp,i))**2)
 116     continue
c     evalue=evalue/nopart(jgrp)
      else
         evalue=0d0
      endif
      return
c v2 for particle group
 88   continue
      evalue=0d0
      if(nopart(jgrp).gt.0) then
         do 117 i=1,nopart(jgrp)
            j=grpart(jgrp,i)
            pt=dsqrt(px(j)**2+py(j)**2)
            if(pt.gt.1d-10) then
               evalue=evalue+(px(j)/pt)**2-(py(j)/pt)**2
            endif
 117  continue
      evalue=evalue/nopart(jgrp)
      endif
      return
c eccentricity via coordinate space
 89   continue
      evalue=0d0
      a1=0d0
      a2=0d0
      if(nopart(jgrp).gt.0) then
         do 118 i=1,nopart(jgrp)
            j=grpart(jgrp,i)
            a1=a1+((ry(j)*ry(j))-(rx(j)*rx(j)))
            a2=a2+((ry(j)*ry(j))+(rx(j)*rx(j)))
 118     continue
         if(abs(a2).gt.1d-10) then
            evalue=a1/a2
         endif
      endif
      return
c dummy observable
 99   write(6,*)'Event-Observable ',obs,' not implemented.'
      write(6,*)'EMERGENCY-STOP!!'
      stop
      end
