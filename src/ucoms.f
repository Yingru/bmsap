cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c standard common block for uout2nt
c     Unit     : all modules
c     Author   : Steffen A. Bass, mofied by L.A.Winckelmann
c     Date     : 05/04/94
c     Revision : 1.0 beta - untested
c
c
      double precision EMNUC
      parameter (EMNUC = 0.938d0)

      integer nmax
      parameter (nmax = 2000000)

c     this parameter MUST match MSTV(3) in BMS
      integer mstv3
      parameter (mstv3 = 1000)

      integer npart, ap,at,zp,zt,ntim,itim,istep
      integer event,iret,evsamp
      double precision  ebeam, bimp,ecm,pbeam,pcm
      double precision betann,betapro,betatar
c 7 integer

      character*8 model_tag, version_tag
      common /names/model_tag,version_tag

      integer seed, ttime, refsys, logscale, addition
      common /upara/ seed,refsys,logscale,addition

      common /sys/ npart,event,ap,zp,at,zt,ntim,itim,iret,istep,evsamp

      common /rsys/bimp,ebeam,ecm,pbeam,pcm,
     +             betapro,betatar,betann
c 2*nmax*nmax logical

      integer ityp(nmax),origin(nmax)
c 6*nmax integer

      double precision eps, er0, pi, rho0
      parameter (eps  = 1.0E-12,
     +           er0  = 1.128379167d0,
     +           pi   = 3.1415926535d0,
     +           rho0 = 0.16d0)


      double precision r0(nmax),rx(nmax), ry(nmax), rz(nmax),
     +     p0(nmax),px(nmax), py(nmax), pz(nmax),
     +     frr0(nmax),frrx(nmax), frry(nmax), frrz(nmax),
     +     fmass(nmax),lstcoll(nmax),frpx(nmax),frpy(nmax),
     +     frpz(nmax),frp0(nmax),Thydro(nmax),weight(nmax),
     +     c_vx(nmax),c_vy(nmax),c_vz(nmax)
      

      common /isys/ ityp, lstcoll, origin
      common /coor/ r0, rx, ry, rz, p0, px, py, pz, fmass,
     &              frr0,frrx,frry,frrz,frpx,frpy,frpz,frp0,
     &              Thydro,weight,c_vx,c_vy,c_vz


      integer count1,count2
      common /sscount/ count1,count2

      double precision tot_wt,tot_obs
      common /average/ tot_wt,tot_obs

c     analysis with weights or not (by Shanshan)
c     only available for 1d analysis now
      integer flag_wt,wt_num,wt_num_MAX
      double precision wt_int,wt_Tab_min
c      parameter (flag_wt = 2) 
          ! 0: don't worry about weight
          ! 1: contain weight in particle list
          ! 2: need to convert initial pT into weight first
      parameter(wt_num_MAX=140)
      double precision pT_wt_table(wt_num_MAX)
      common/wtTable1/flag_wt,wt_num
      common/wtTable2/wt_int,wt_Tab_min
      common/wtTable3/pT_wt_table

