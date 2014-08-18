      subroutine obsvalue(obs,inpart,retvec)
c     this function returns the value of observable obs of particle idpart
c     
c     author  : Steffen A. Bass
c     date    : 17.02.94
c     revision: 1.1
c

      IMPLICIT REAL*8(A-H,O-Z)

      integer obs,inpart
      parameter(degrad=180./3.1415927)

c
      include 'ucoms.f'
      include 'hicoms.f'

      double precision pt(nmax),etot(nmax),rap(nmax),phi(nmax),
     &                 virt(nmax)
      integer am(nmax),charge,itmp
      double precision theta(nmax),retvec(nmax),etra(nmax),ppz,rapp
c
c first calculate some global observables (all in nn-system)
c
c calculate projectile rapidity, assuming a nuclear projectile
      en=EMNUC
      ppz=pcm
      entot=dsqrt(en*en+ppz*ppz)
      rapp=0.5*dlog((entot+ppz) / (entot-ppz))

      al=8.66d0
      sal=21.0d0
c
      do 999 i=1,inpart
c
      pt(i)=sqrt(px(i)*px(i)+py(i)*py(i))
c      etot(i)=sqrt(fmass(i)*fmass(i)+px(i)*px(i)
c     &     +py(i)*py(i)+pz(i)*pz(i))
      etot(i)=p0(i)
      rap(i)=0.5*dlog( (etot(i)+pz(i)) / (etot(i)-pz(i)) )
      phi(i)=90.-degrad*atan(px(i)/max(abs(py(i)),1d-8) )
      theta(i)=90.-degrad*atan(pz(i)/max(pt(i),1d-8) )
      am(i)=int((fmass(i)+0.01)/.938d0)
      virt(i)=sign(fmass(i)*fmass(i),fmass(i))
      if(am(i).lt.1) am(i)=1

      if(lstcoll(i).gt.0d0) then
         ncoll=NINT(lstcoll(I))/MSTV3**2
         nsbra=MOD(NINT(lstcoll(I))/MSTV3,MSTV3)
         ntbra=MOD(NINT(lstcoll(I)),MSTV3)
      else
         ncoll=0
         nsbra=0
         ntbra=0
      endif

c
c now goto the desired observable
      obs=abs(obs)
c           1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
      goto( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,
     &     21,22,23,24,25,26,27,28,29,30,99,99,99,99,35,36,37,38,39,40,
     &     41,42,99,99,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,
     &     61,62,63,64,65,66,67,68,69,99,99,99,99,99,99,99,99,99,99,99,
     &     99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,101
     & ), obs
      goto 99
c
c here come the observables
c feel free to add your favourite observable;
c in case of changes, please contact the authors 
c for a general revision update.
c
c dN/d? distribution
 1    retvec(i)=1.d0
      goto 999
c px
 2    retvec(i)=px(i)
      goto 999
c py
 3    retvec(i)=py(i)
      goto 999
c pz
 4    retvec(i)=pz(i)
      goto 999
c pt
 5    retvec(i)=pt(i)
      goto 999
c rapidity
 6    retvec(i)=rap(i)
      goto 999
c E_tot
 7    retvec(i)=etot(i)
      goto 999
c phi
 8    retvec(i)=phi(i)
      goto 999
c theta
 9    retvec(i)=theta(i)
      goto 999
c abs(p)
 10   retvec(i)=sqrt(pt(i)*pt(i)+pz(i)*pz(i))
      goto 999
c impossible Normalization
 11   write(6,*)'wrong column selected for normalization'
      write(6,*)'EMERGENCY-STOP!!!'
      stop
c mass
 12   retvec(i)=fmass(i)
      goto 999
c transverse mass
 13   retvec(i)=sqrt(fmass(i)*fmass(i)+pt(i)*pt(i))
      goto 999
c raw r(i,5) information
 14   retvec(i)=lstcoll(i)
      goto 999
c three times electric charge
 15   retvec(i)=dble(charge(ityp(i)))
      goto 999
c ityp of particle
 16   retvec(i)=dble(ityp(i))
      goto 999
c number of (hard) collisions
 17   retvec(i)=dble(ncoll)
      goto 999
c number of space-like branchings
 18   retvec(i)=dble(nsbra)
      goto 999
c number of time-like branchings
 19   retvec(i)=dble(ntbra)
      goto 999
c flow coefficient v2 (only useful when seen as average over many particles)
c (px/pt)^2-(py/pt)^2 flow coefficient v2
 20   if(pt(i).le.1d-10) then
         retvec(i)=0d0
      else
         retvec(i)=(px(i)/pt(i))**2-(py(i)/pt(i))**2
      endif
      goto 999
c virtuality
 21   retvec(i)=virt(i)
      goto 999
c time coordinate
 22   retvec(i)=r0(i)
      goto 999
c x coordinate
 23   retvec(i)=rx(i)
      goto 999
c y coordinate
 24   retvec(i)=ry(i)
      goto 999
c z coordinate
 25   retvec(i)=rz(i)
      goto 999
 26   retvec(i)=rap(i)/rapp
      goto 999
 27   retvec(i)=px(i)/ppz
      goto 999
 28   retvec(i)=pt(i)/ppz
      goto 999
c Feynman x
 29   retvec(i)=pz(i)/ppz
      goto 999
c qzz
 30   retvec(i)=2*pz(i)*pz(i)-pt(i)*pt(i)
      goto 999
c px/A
 35   retvec(i)=px(i)/dble(am(i))
      goto 999
c py/A
 36   retvec(i)=py(i)/dble(am(i))
      goto 999
c pz/A
 37   retvec(i)=pz(i)/dble(am(i))
      goto 999
c pt/A
 38   retvec(i)=pt(i)/dble(am(i))
      goto 999
c mt-m0
 39   retvec(i)=sqrt(fmass(i)*fmass(i)+pt(i)*pt(i))-fmass(i)
      goto 999
c particle index
 40   retvec(i)=i
      goto 999
c transverse energy
 41   retvec(i)=sqrt(pt(i)*pt(i)+fmass(i)*fmass(i))
      goto 999
c transverse pressure
 42   retvec(i)=pt(i)*pt(i)/sqrt(fmass(i)*fmass(i)+px(i)*px(i)
     &     +py(i)*py(i)+pz(i)*pz(i))
      goto 999
c time at freezeout
 45   retvec(i)=frr0(i)
      goto 999
c rx at freeze out
 46   retvec(i)=frrx(i)
      goto 999
c ry at freeze out
 47   retvec(i)=frry(i)
      goto 999
c rz at freeze out
 48   retvec(i)=frrz(i)
      goto 999
c abs(r) at freeze out
 49   retvec(i)=sqrt(frrx(i)*frrx(i)+frry(i)*frry(i)+
     & frrz(i)*frrz(i))
      goto 999
c  abs(r_t) at freeze out
 50   retvec(i)=sqrt(frrx(i)*frrx(i)+frry(i)*frry(i))
      goto 999
c  abs(r_t) 
 51   retvec(i)=sqrt(rx(i)*rx(i)+ry(i)*ry(i))
      goto 999
c r_trans*p_trans/|p_trans|
 52   retvec(i)=(frrx(i)*px(i)+frry(i)*py(i))
     &     /sqrt(px(i)**2+py(i)**2)
      goto 999
c pseudo-rapidity
 53   retvec(i)=-1.0*log(
     &     tan(acos(pz(i)/sqrt(pt(i)*pt(i)+pz(i)*pz(i)))/2.0))
      goto 999
c log of Feynman x
 54   retvec(i)=dlog10(pz(i)/ppz)
      goto 999
c initial rapidity
 55   if(pz(i).ge.0d0) then
         retvec(i)=rapp+dlog(pz(i)/ppz)
      else
         retvec(i)=-1d0*(rapp+dlog(dabs(pz(i))/ppz))
      endif
      goto 999
c pz*pz
 56   retvec(i)=pz(i)*pz(i)
      goto 999
c pt*pt
 57   retvec(i)=pt(i)*pt(i)
      goto 999
c Thydro: temperature of the medium
 58   retvec(i)=Thydro(i)
      goto 999
c initial transverse momentum
 59   retvec(i)=sqrt(frpx(i)*frpx(i)+frpy(i)*frpy(i))
      goto 999
c initial energy (be careful, fmass may not correspond to initial particle)
 60   retvec(i)=sqrt(frpx(i)*frpx(i)+frpy(i)*frpy(i)
     &       +frpz(i)*frpz(i)+fmass(i)*fmass(i))
      goto 999
c weight of the initial particle
 61   retvec(i)=weight(i)
      goto 999

c flow coefficient v4 (only useful when seen as average over many particles)
 62   if(pt(i).le.1d-10) then
         retvec(i)=0d0
      else
         retvec(i)=1d0-8d0*(px(i)*py(i)/pt(i)/pt(i))**2
      endif
      goto 999

c flow coefficient v4 (only useful when seen as average over many particles)
 63   if(pt(i).le.1d-10) then
         retvec(i)=1d0
      else
         retvec(i)=1d0-(2d0*px(i)*py(i)/(px(i)*px(i)-py(i)*py(i)))**2
      endif
      goto 999

c v_T of cell
 64   retvec(i)=sqrt(c_vx(i)**2+c_vy(i)**2)
      goto 999

c v2 of cell velocity
 65   if(c_vx(i)**2+c_vy(i)**2.lt.1d-10) then
         retvec(i)=0d0
      else
         retvec(i)=(c_vx(i)**2-c_vy(i)**2)/(c_vx(i)**2+c_vy(i)**2)
      endif
      goto 999

 66   retvec(i)=c_vx(i)
      goto 999

 67   retvec(i)=c_vy(i)
      goto 999

 68   retvec(i)=c_vz(i)
      goto 999

 69   retvec(i)=c_vx(i)**2-c_vy(i)**2
      goto 999


c here comes a minimum spanning tree to distinguish between bound/unbound part.
 101  if(retvec(i).lt.0.) goto 999
      retvec(i)=1.
      ii=i+1
      do 888 j=ii,inpart
         deltar2=(rx(i)-rx(j))**2+(ry(i)-ry(j))**2+(rz(i)-rz(j))**2
         if(deltar2.lt.9.) then
            retvec(i)=-1.
            retvec(j)=-1.
            goto 999
         endif
 888  continue
      goto 999
c observable not implemented
 99   write(6,*)'observable ',obs,' not implemented!!!'
      write(6,*)'EMERGENCY-STOP!!'
      stop
 999  continue
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function charge(pdgid)

c this function returns 3 times the electric charge of the particle

      implicit none
      integer pdgid,qcharge(6),anti,i,digits(5),itmp
      
c                   d u  s c  b t
      data qcharge/-1,2,-1,2,-1,2/

      charge=0

c     treat gluon and photon explicitely
      if(pdgid.eq.21.or.pdgid.eq.22) then
         charge=0
         return
      endif

      anti=1
      if(pdgid.lt.0) anti=-1

c     (anti-)quarks
      if(iabs(pdgid).le.6) then
        charge=anti*qcharge(iabs(pdgid))
        return
      endif

c     hadrons:

c     dissect ID into individual digits

      itmp = iabs(pdgid)
      do i = 5,1,-1
         if ( itmp .lt. 10**(i-1)) then
            digits(i) = 0
         else 
            digits(i) = itmp/(10**(i-1))
            itmp = itmp - digits(i)*10**(i-1)
         end if 
      enddo


c     return if trouble:
      if((pdgid.lt.100).or.(pdgid.gt.1000.and.digits(2).eq.0)) then
         write(6,*) 'trouble in charge routine ',pdgid
         stop
      endif

c     baryon or meson:
      if(digits(4).ne.0) then ! baryon
         do 20 i=2,4
            charge=charge+qcharge(digits(i))
 20      continue
      else ! meson
         if(mod(digits(3),2).ne.0) then ! left quark is antiquark
            charge=charge+qcharge(digits(2))
            charge=charge-qcharge(digits(3))
         else ! right quark is antiquark
            charge=charge-qcharge(digits(2))
            charge=charge+qcharge(digits(3))
         endif
      endif

      charge=anti*charge

      return
      end







