c
c common blocks for ubap:
c

      integer nbgr
      parameter(nbgr=200)

      logical allpart
      character*2 cflag
      character*4 model
      character*78 header
      integer cuts,grpsel,xaxis,zaxis,splts
      integer ecuts
      real*8 rotphi,rottheta,bmean,bbmin,bbmax
      integer cutobs(30),cutgrp(30),xobs(30),ecutobs(30),normflag(30)
      integer xbins(30),yobs(30),ybins(30),scatx(30),ecutgrp(30)
      integer zobs(30),zunit(30),zgrp(30),scaty(30),scatz(30)
      integer scatunit(30),scatgrp(30),angrp(nbgr,30),totev,anev,tstep
      real*8 ecutmin(30),ecutmax(30)
      real*8 cutmin(30),cutmax(30),xmin(30),xmax(30),ymin(30)
      real*8 ymax(30),var(60,30)
      real*8 vavalc(nmax),vavalx(nmax),vavaly(nmax),vavalz(nmax)
      real*8 vawt(nmax)
      real*8 vavalnx(nmax),vavalny(nmax)
      real*8 zxygrd(30,nbgr,nbgr),errgrd(30,nbgr,nbgr)
      real*8 dnsgrd1(30,nbgr,nbgr)
      real*8 dxbin(30),dybin(30)
      integer dnsgrd(30,nbgr,nbgr),wsum(30),wpart(30),evcntr(30,30)
      integer cutarr(30,nmax),nopart(30),grpart(30,nmax)
      integer encntr(30,30)
      common /rots/rotphi,rottheta,bmean,bbmin,bbmax
      common /strings/cflag,header
      common /ecut1/ecuts,ecutobs,anev
      common /ecut2/ecutmin,ecutmax,ecutgrp
      common /fflw/var
      common /cut/cutobs,cutmin,cutmax,cutgrp,cutarr,cuts
      common /parts/nopart,grpart,evcntr,totev,encntr
      common /grps/grpsel,angrp,normflag,allpart
      common /xa/xobs,xmin,xmax,xbins,dxbin,xaxis
      common /ya/yobs,ymin,ymax,ybins,dybin
      common /za/zaxis,zobs,zgrp,zunit,tstep
      common /sp/splts,scatx,scaty,scatz,scatgrp,scatunit
      common /grids/zxygrd,errgrd,dnsgrd,dnsgrd1
      save

