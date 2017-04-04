      subroutine input
c     this routine handles the input of hicap
c     author  : Steffen A. Bass
c     date    : 04.01.96
c     revision: 0.1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT REAL*8(A-H,O-Z)

      character*2 flag
      character*78 inputstr,file9
      integer ffflag,line
      include 'ucoms.f'
      include 'hicoms.f'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c initialize all arrays
      do 10 i=1,30
         cutobs(i)=0
         ecutobs(i)=0
         xobs(i)=0
         yobs(i)=0
         zobs(i)=0
         scatx(i)=0
         scaty(i)=0
         scatz(i)=0
         xbins(i)=0
         ybins(i)=0
         zunit(i)=0
         cutgrp(i)=0
         ecutgrp(i)=0
         zgrp(i)=0
         cutmin(i)=0.0
         ecutmin(i)=0.0
         xmin(i)=0.0
         ymin(i)=0.0
         cutmax(i)=0.0
         ecutmax(i)=0.0
         xmax(i)=0.0
         ymax(i)=0.0
         nopart(i)=0
         normflag(i)=0
         do 20 j=1,30
               angrp(j,i)=0
               evcntr(j,i)=0
  20      continue
         do 10 j=1,nparmx
            cutarr(i,j)=0
            grpart(i,j)=0
 10   continue
c initialize counters
      allpart=.false.
      line=1
      ffflag=0
      cuts=0
      ecuts=0
      grpsel=0
      xaxis=0
      zaxis=0
      splts=0
c initialize the output scale
      logscale=0
      addition=0
c initialize weight method
      flag_wt=0
      wt_num=0
      wt_int=0.5d0
      wt_Tab_min=0.5d0
c open fortran-unit 9 for input
      call getenv('ftn09',file9)
      OPEN(UNIT=9,FILE=file9,STATUS='OLD',FORM='FORMATTED')
c read header
      read(9,99) cflag,header
 99   format(1A2,1A78)
c read input lines
 1    continue
      line=line+1
      read(9,99) flag,inputstr
 3    continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c select action according to flag:
c # : treat line as a comment
      if(flag.eq.'# ') goto 1
c xx: treat line as end of input marker
      if(flag.eq.'xx') goto 2
c ff: input-flags: refsys, timestep  output-flag: logscale,addition
      if(flag.eq.'ff') then
         ffflag=ffflag+1
         read(inputstr,fmt=*,err=88,end=88) refsys,tstep,logscale,
     *       addition 
         if(ffflag.gt.1) then
            write(6,*)'ff-parameters are defined more than once!!'
            write(6,*)'-> last definition is used'
         endif
      elseif(flag.eq.'wt') then
         read(inputstr,fmt=*,err=88,end=88) flag_wt,wt_num,wt_int,
     &                                      wt_Tab_min
c cu: perform a cut
      elseif(flag.eq.'cu') then
         cuts=cuts+1
         read(inputstr,fmt=*,err=88,end=88) cutobs(cuts),
     *       cutmin(cuts),cutmax(cuts),cutgrp(cuts)
c     now check, if cut refers to event-like observables
c     event-like cuts are stored and handled extra
         if(cutobs(cuts).gt.100) then
            ecuts=ecuts+1
            ecutobs(ecuts)=cutobs(cuts)
            ecutmin(ecuts)=cutmin(cuts)
            ecutmax(ecuts)=cutmax(cuts)
            ecutgrp(ecuts)=cutgrp(cuts)
            cuts=cuts-1
         endif
c gr: definitions for analysis group 1-4
      elseif(flag.eq.'gr') then
         grpsel=grpsel+1
         if(grpsel.gt.30) then
            write(6,*)'maximum # of groups exceeded, input scipped!'
            grpsel=grpsel-1
 5          continue
            line=line+1
            read(9,99) flag,inputstr
            if(flag.eq.'# ') goto 5
            if(flag.ne.'tp') goto 3
            goto 1
         endif
c     angrp(grpsel,29) contains the groupflag (in/out of cut etc.)
c     angrp(grpsel,30) contains the number of different itypes in the group
         read(inputstr,fmt=*,err=88,end=88) angrp(grpsel,29),
     *       angrp(grpsel,30)
         if(angrp(grpsel,30).eq.-1) then
            allpart=.true.
            write(6,*)'distinctive particle groups disabled'
            write(6,*)'-> all particles used for analysis'
            write(6,*)'   (use cuts for selection)'
            write(6,*)'-> all other gr/tp/i3 statements ignored!'
            angrp(1,29)=angrp(grpsel,29)
            angrp(1,30)=angrp(grpsel,30)
            goto 1
         endif
c tp: definitions for particle-types in group 1-30
 6       continue
         line=line+1
         read(9,99) flag,inputstr
            if(flag.eq.'# ') goto 6
            if(flag.ne.'tp') then
               write(6,*)'tp-declaration must follow gr-flag!!'
               write(6,*)'-> input aborted on line ',line
               stop
            endif
c     it is mandatory that a tp-flag follows a gr-flag
         read(inputstr,fmt=*,err=88,end=88) 
     *       (angrp(grpsel,i),i=1,angrp(grpsel,30))
c tp-error:
      elseif(flag.eq.'tp')then
            write(6,*)'warning: all tp-lines not following',
     * ' directly a gr line are ignored!'
            write(6,*)'-> scipping tp-definition on line ',line
c xa: definitions for x-axis of grid
      elseif(flag.eq.'xa') then
         xaxis=xaxis+1
         read(inputstr,fmt=*,err=88,end=88) xobs(xaxis),
     *       xmin(xaxis),xmax(xaxis),xbins(xaxis)
         if(xbins(xaxis).gt.NBGR)then
            write(6,*)'Warning: NBGR exceeded on line ',line
            xbins(xaxis)=NBGR
         endif
c ya: definitions for y-axis of grid
      elseif(flag.eq.'ya') then
         if(xaxis.eq.0) then
            write(6,*)'fatal error: xa must be defined',
     * ' before ya!'
            stop
         endif   
         read(inputstr,fmt=*,err=88,end=88) yobs(xaxis),
     *       ymin(xaxis),ymax(xaxis),ybins(xaxis)
         if(ybins(xaxis).gt.NBGR)then
            write(6,*)'Warning: NBGR exceeded on line ',line
            ybins(xaxis)=NBGR
         endif
c za: definitions for z-axis of grid
      elseif(flag.eq.'za') then
         ix=zaxis
         zaxis=zaxis+1
         read(inputstr,fmt=*,err=88,end=88) zobs(zaxis),
     *       zgrp(zaxis),zunit(zaxis)
         if(abs(zobs(zaxis)).eq.11) then
            normflag(zaxis)=zobs(zaxis)
            zobs(zaxis)=1
         endif   
         if(xaxis.eq.0) then
            write(6,*)'fatal error: xa must be defined',
     * ' before za!'
            stop
         elseif(ix.eq.xaxis) then
c     use old x- and yaxis
            xaxis=xaxis+1
            xobs(xaxis)=xobs(ix)
            xmin(xaxis)=xmin(ix)
            xmax(xaxis)=xmax(ix)
            xbins(xaxis)=xbins(ix)
            yobs(xaxis)=yobs(ix)
            ymin(xaxis)=ymin(ix)
            ymax(xaxis)=ymax(ix)
            ybins(xaxis)=ybins(ix)
         endif
c sp: definitions for scatter-plots
      elseif(flag.eq.'sp') then
         splts=splts+1
         read(inputstr,fmt=*,err=88,end=88) scatx(splts),
     *   scaty(splts),scatz(splts),scatgrp(splts),scatunit(splts)
      else
         write(6,*)'undefined opcode in input-file on line',line
         stop
      endif
      goto 1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 2    continue
c here some validity checks of the input should be performed
      if(model.eq.'baeh') then
         write(6,*)'fatal error: no model specified for input!!'
         stop
      elseif(ffflag.lt.1) then
         write(6,*)'fatal error: no ff-definitions in input!!'
         stop
      elseif((zaxis.eq.0).and.(splts.eq.0)) then
         write(6,*)'fatal error: no analysis selected!'
c hier Aenderungen wenn GEANT-output gewuenscht
         stop
      elseif(grpsel.lt.1) then
         write(6,*)'fatal error: no particle groups selected!'
         stop
      endif
c now print the selected analysis
      write(6,*) '************ This is BMS_AP 1.0 beta  ************'
      write(6,*)
      write(6,*) '(c) by S.A. Bass 2002'
      write(6,*) 
      write(6,*) header
      write(6,*)
      write(6,*)
      write(6,*)'Program parameters:'
      write(6,207) refsys,tstep,logscale
      write(6,209) flag_wt,wt_num,wt_int,wt_Tab_min
c      write(6,208) rotphi,rottheta
 207  format(1X,'Refsys:',I2,' Timestep:',I3,' Logscale:',I2)
 209  format(1X,'Weight Type:',I2,' wt_num:',I4,' wt_int:',F6.3,
     &       ' min_pT:',F6.3)
 
c 208  format(1X,'Phi-Rot:',F8.2,' Theta-Rot:',F8.2)
      write(6,*)
      write(6,*)'selected groups are: '
      if(allpart) then
         write(6,*)' all particles in one group'
         write(6,*)' -> particle selection via cuts'
      else
         do 104 i=1,grpsel
            write(6,206) i,angrp(i,29),
     *           (angrp(i,j),j=1,angrp(i,30))
 104     continue
 206     format(1X,'Grp:',I2,' Flag:',I3,' Itypes:',28I5)
      endif
      write(6,*)
      if(cuts.eq.0) then
         write(6,*) 'no cuts selected'
      else
         write(6,*) 'selected cuts: '
         do 101 i=1,cuts
            write(6,201) cutobs(i),cutmin(i),cutmax(i),cutgrp(i)
 101     continue
      endif   
 201  format(1X,'Obs:',I3,' Min:',F8.2,' Max:',F8.2,' Grp:',I3)
      write(6,*)
      if(ecuts.eq.0) then
         write(6,*) 'no event-like cuts selected'
      else
         write(6,*) 'selected event-cuts: '
         do 302 i=1,ecuts
            write(6,201) ecutobs(i),ecutmin(i),ecutmax(i),ecutgrp(i)
 302     continue
      endif   
      write(6,*)
      if(zaxis.eq.0)then
         write(6,*) 'no 1D histos / 2D grids selected'
      else
         write(6,*) 'selected 1D histos /2D grids:'
         write(6,*)
         do 102 i=1,zaxis
            write(6,202) xobs(i),xmin(i),xmax(i),xbins(i)
            if(ybins(i).eq.0) goto 11
            write(6,203) yobs(i),ymin(i),ymax(i),ybins(i)
 11         write(6,204) zobs(i),zgrp(i),normflag(i),zunit(i)
            write(6,*)
 102     continue
      endif
 202     format(1X,'Xobs:',I3,' Min:',F8.2,' Max:',F8.2,' Bins:',I3)
 203     format(1X,'Yobs:',I3,' Min:',F8.2,' Max:',F8.2,' Bins:',I3)
 204     format(1X,'Zobs:',I3,' Grp:',I2,' Norm:',I3,' Unit:',I3)
      write(6,*)   
      if(splts.eq.0)then
         write(6,*)'no scatter-plots selected'
      else
         write(6,*)'selected scatter-plots:'
         write(6,*)
         do 103 i=1,splts
            write(6,205) scatx(i),scaty(i),scatz(i),
     * scatgrp(i),scatunit(i)
            write(6,*)
 103     continue
 205     format(1X,'Xobs:',I3,' Yobs:',I3,' Zobs:',
     * I3,' Grp:',I2,' Unit:',I3)
      endif
      write(6,*)
      write(6,*) 'end of input'
      write(6,331) logscale
 331    format(1X,'logscale in input file:', I2)
c 
      return
c error-exit
 88   write(6,*) 'syntax-error in input-file on line ',line
      stop
      end







