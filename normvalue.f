      subroutine normvalue(obs,inpart,flag,retvec)
c     this function returns the value of observable obs of particle idpart
c     
c     author  : Steffen A. Bass
c     date    : 14.09.93
c     revision: 1.0
c
      IMPLICIT REAL*8(A-H,O-Z)

      integer obs,inpart,flag
      parameter(degrad=180./3.1415927)

c
      include 'ucoms.f'
      include 'hicoms.f'

      real*8 pt(nmax),etot(nmax),rap(nmax),phi(nmax)
      real*8 theta(nmax),retvec(nmax)

c
c first calculate some global observables (all in nn-system)
c
      do 999 i=1,inpart
c
         ppz=pcm
         if(obs.lt.100) then
            pt(i)=dsqrt(px(i)*px(i)+py(i)*py(i))
c            etot(i)=dsqrt(fmass(i)*fmass(i)+
c     &           px(i)*px(i)+py(i)*py(i)+pz(i)*pz(i))
            etot(i)=p0(i)
            rap(i)=0.5*dlog( (etot(i)+pz(i)) / (etot(i)-pz(i)) )
            phi(i)=90.-degrad*atan(px(i)/max(abs(py(i)),1e-8) )
            theta(i)=90.-degrad*atan(pz(i)/max(pt(i),1e-8) )
         endif
c
c now goto the desired observable
c           1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
      goto( 1, 1, 1, 1, 5, 1, 7, 1, 9,10, 1, 1,13, 1, 1, 1, 1, 1, 1, 1,
     &      1, 1, 1, 1, 1, 1, 1, 1,29, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     &      1, 1, 1, 1, 1, 1, 1, 1, 1,50, 1, 1, 1,54, 1, 1, 1, 1, 1, 1,
     &      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     &      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     &      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     &      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     &      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     &      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
     &), obs
      goto 1
c
c here come the observables
c feel free to add your favourite observable;
c in case of changes, please contact the authors 
c for a general revision update.
c
c dN/d? distribution
 1    retvec(i)=1.
      goto 999
c pt
 5    if(flag.gt.0) then
         retvec(i)=1./pt(i)
      else
         retvec(i)=1./(pt(i)*pt(i))
      endif
      goto 999
c E_tot
 7    if(flag.gt.0)then
         retvec(i)=1./etot(i)
      else
         retvec(i)=1./(etot(i)*sqrt(pt(i)*pt(i)+pz(i)*pz(i)))
      endif
      goto 999
c theta
 9    retvec(i)=sqrt(pt(i)*pt(i)+pz(i)*pz(i))/pt(i)
      goto 999
c abs(p)
 10   if(flag.gt.0)then
         retvec(i)=1./(pt(i)*pt(i)+pz(i)*pz(i))
      else
         retvec(i)=etot(i)/(pt(i)*pt(i)+pz(i)*pz(i))
      endif
      goto 999
 13   retvec(i)=1./sqrt(fmass(i)*fmass(i)+pt(i)*pt(i))
      goto 999
 29   retvec(i)=pz(i)/ppz
      goto 999
 50   retvec(i)=1./sqrt(frrx(i)*frrx(i)+frry(i)*frry(i))
      goto 999
 54   retvec(i)=1./(pz(i)/ppz)
      goto 999
 999  continue
      return
      end









