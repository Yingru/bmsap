cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program bms_ap
c     
c     this program is a multi-purpose analysis program for modified BMS
c     OSCAR output
c     For details in running the program please consult the manual bms_ap.tex.
c     authors: 
c     S.A. Bass     bass@phy.duke.edu
c
c     date 03/21/2002
c     revision 1.0 beta
c
c     this version is not yet fully tested!!!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT NONE

      double precision evalue,axval,azval
      include 'ucoms.f'
      include 'hicoms.f'
      character*5 envname
      character*15 ostr
      character*80 fnam
      integer nsteps, iii
      integer dum_wt,ipn,igrp,igpart,iun
      integer ixbin,iybin,i,j,isp,inf,iiflag
      integer iana,ibookflag,icuts,idpart
      integer iecuts,ievflag,ig,igsp
      double precision aweight,acv,ayval
      double precision xx,yy,zz
      integer wt_table

c first read input-file 
      call input
c now initialize grids for analysis
      call gridin
      totev=0
      anev=0
      iret=99

      dum_wt=61

      tot_wt=0.d0
      tot_obs=0.d0

      count1=0 ! count how many particle have been analyzed
      count2=0 ! count the # of selected particles

c debug
      if(flag_wt.ne.0.and.flag_wt.ne.1.and.flag_wt.ne.2) then
         write(6,*) "Wrong value for flag_wt ..."
         write(6,*) "Terminating ..."
         stop
      endif

      if(flag_wt.eq.2) then
         if(wt_num.gt.wt_num_MAX) then
            write(6,*) "wt_num exceeds maximum value allowed."
            write(6,*) "Terminating ..."
            stop
         endif
         wt_table=0
         call read_wt_table(wt_table)
         if(wt_table.ne.0) then
            write(6,*) "error occurs during reading weight table ..."
            write(6,*) "terminating ..."
            stop
         endif
      endif

c read in the ration angle (may be participant plane, may be the event plane)
c added by Yingru
      cos_pp = 1d0
      sin_pp = 0d0
      call read_rotation_angle
      write(6,*) cos_pp, sin_pp

c define one group in case of allpart
      if(allpart) grpsel=1
c totev: total events; anev: analysed events

c read file header:
 3    continue

      istep=1
      call read_osc_event
      if(iret.eq.0) then
         write(6,332) logscale
 332     format('logscale in main code:', I2)        
         write(6,*) '*** end of event-input ***'
         call output
         stop
      elseif(iret.lt.0) then
         write(6,*) '**** error during input in read_osc_event ***'
         write(6,*) '**** last event was skipped, terminating...***'
         call output
         stop
      end if

      nsteps=ntim

c      write(6,*) ' steps ',tstep,nsteps

      if(tstep.gt.1) then ! the first timestep is already in ucoms
         do 33 iii=2,(tstep-1)*evsamp
c            write(6,*) ' read next ',iii
            istep=iii
            call gettimestep
            if(iret.eq.0) then
               write(6,*) '*** end of event-input ***'
               call output
               stop
            elseif(iret.eq.-1) then
             write(6,*) '**** error during event-input in uqmdnxtev ***'
             write(6,*) '**** last event was scipped, terminating...***'
             call output
             stop
            end if
 33      continue
      endif

      do 933 iii=(tstep-1)*evsamp+1,tstep*evsamp
         istep=iii
         if(istep.ne.1) call gettimestep
         if(iret.eq.0) then
            write(6,*) '*** end of event-input ***'
            call output
            stop
         elseif(iret.eq.-1) then
           write(6,*) '**** error during event-input in uqmdnxtev ***'
           write(6,*) '**** last event was scipped, terminating...***'
           call output
           stop
         end if

         totev=totev+1
c calculate some additional info
         write(6,*)'processing event # ',totev,' timestep # ',tstep
         write(6,*)'true event # ',event,' timestep # ',itim
         call reftrans
         bmean=bmean+bimp

         if(npart.gt.nmax) then
            write(6,*)'warning: too many particles in event ',totev
            write(6,*)'terminating...'
c            write(6,*)'event is skipped !!!'
            stop
            goto 3
         endif

         count1=count1+npart
 
c determine, if particles are in cuts (without regard of grouping)
         if(cuts.eq.0) goto 13
         do 10 icuts=1,cuts
            call obsvalue(cutobs(icuts),npart,vavalc)
            do 10 idpart=1,npart
               if(vavalc(idpart).lt.cutmin(icuts)) then
                  cutarr(icuts,idpart)=-1
               elseif(vavalc(idpart).gt.cutmax(icuts)) then
                  cutarr(icuts,idpart)=1
               else
                  cutarr(icuts,idpart)=0
               endif
 10      continue
 13      continue
c now fill the analysis-groups 
         do 20 igrp=1,grpsel
            nopart(igrp)=0
            do 30 idpart=1,npart
c ibookflag makes certain, that every particle is booked only
c once into a analysis group
               ibookflag=-1
c first find out, if ityp(idpart) is in group igrp
               iiflag=0
               if(allpart) then
                  iiflag=1
               else
                  do 35 igpart=1,angrp(igrp,30)
                     if((ityp(idpart)).eq.angrp(igrp,igpart)) 
     &                   iiflag=iiflag+1
 35               continue
               endif
               if(iiflag.eq.0) goto 30
c book particles if no cuts are selected
               if(cuts.eq.0) then
                  ibookflag=0
               else
                  ibookflag=0
                  do 40 icuts=1,cuts
c book particles which fall into a cut
                     if((cutarr(icuts,idpart).eq.angrp(igrp,29)).and.
     *           ((cutgrp(icuts).eq.igrp).or.(cutgrp(icuts).eq.0)))then
                     ibookflag=ibookflag+1
c book particles with grpflag=2
                     elseif((angrp(igrp,29).eq.2).and.
     *                  (abs(cutarr(icuts,idpart)).eq.1).and.
     *           ((cutgrp(icuts).eq.igrp).or.(cutgrp(icuts).eq.0)))then
                     ibookflag=ibookflag+1
c now particles should be booked which do not fall into a cut
                     elseif((cutgrp(icuts).ne.igrp).and.
     *                     (cutgrp(icuts).ne.0))then
                     ibookflag=ibookflag+1
                     else
                       goto 30
                     endif
 40               continue
c end of icuts loop
               endif
c now determine whether particle meets all cuts
c and book them into a pointer array grpart for further use
               if(ibookflag.eq.cuts)then
                  nopart(igrp)=nopart(igrp)+1
                  grpart(igrp,nopart(igrp))=idpart
               else
                  write(6,*) 'ibookflag - error'
                  stop
               endif
               
               count2=count2+1

 30         continue
c end of idpart loop
 20      continue
c end of igrp loop
c
c calculate all event-like observables
         call iniflw
c
c now perform event-like cuts
         if(ecuts.gt.0)then
            do 50 iecuts=1,ecuts
               acv=evalue(ecutobs(iecuts),ecutgrp(iecuts))
               if((acv.lt.ecutmin(iecuts)).or.
     &            (acv.gt.ecutmax(iecuts))) goto 31
c if event is outside of cut, skip analysis and get next event
 50         continue
         endif
c
         anev=anev+1
c now make histos and grids
c
         do 100 iana=1,zaxis
            if(yobs(iana).eq.0) then
c make a 1d-histo for all groups
              do 110 igrp=1,grpsel
c if eventobs. then do histo without particle loop
                 if(xobs(iana).gt.100) then
c first check, if both obs. refer to events
c dn/dobs has different code for event-obs
                    if(zobs(iana).eq.1) zobs(iana)=101
                    if(zobs(iana).lt.100) then
                       write(6,*)'warning: no mixing of event-
     & observables with particle observables!!'
                       goto 100
                    endif
c sab
c                    axval=evalue(xobs(iana),zgrp(iana))
                    axval=evalue(xobs(iana),igrp)
c so far no normalization for eventlike observables
c               inf=normflag(iana)
c               if(inf.ne.0) then
c                  azval=evalue(zobs(iana),igrp)*
c     &                  normvalue(xobs(iana),zgrp(iana),inf)
c               else
                    azval=evalue(zobs(iana),igrp)
c               endif
                    
c add weight for analysis (by Shanshan)
                    if(flag_wt.ne.0) then
                       write(6,*) "weight not available here yet ..."
                       write(6,*) "terminating ..."
                       stop
                    endif

                    ixbin=int((axval-xmin(iana))/dxbin(iana)+1.5)
                    if((ixbin.lt.1).or.(ixbin.gt.xbins(iana))) goto 110
c     count contributing events
                    evcntr(iana,igrp)=evcntr(iana,igrp)+1
                    encntr(iana,igrp)=encntr(iana,igrp)+1
                    zxygrd(iana,ixbin,igrp)=zxygrd(iana,ixbin,igrp)+
     &                   azval
                    dnsgrd(iana,ixbin,igrp)=dnsgrd(iana,ixbin,igrp)+1
                    errgrd(iana,ixbin,igrp)=errgrd(iana,ixbin,igrp)+
     &                   azval*azval
                
              elseif(xobs(iana).lt.100) then 
c     particle observables:
c     loop over all particles in group
                 ievflag=0
                 inf=normflag(iana)
                 call obsvalue(xobs(iana),npart,vavalx)
                 call obsvalue(zobs(iana),npart,vavalz)
                 call obsvalue(dum_wt,npart,vawt)
 
                 if(inf.ne.0) then
                    call normvalue(xobs(iana),npart,inf,vavalnx)
                 endif
                 do 120 igpart=1,nopart(igrp)
c     get particle-pointer
                    idpart=grpart(igrp,igpart)
c     get x and z values of observble
                    axval=vavalx(idpart)
                    azval=vavalz(idpart)
                    aweight=vawt(idpart)
c                    tot_wt=tot_wt+aweight
                    if(inf.ne.0) then
                       azval=azval*vavalnx(idpart)
                    endif
c     check whether particle is in interval
                    ixbin=int((axval-xmin(iana))/dxbin(iana)+1.5)
                    if((ixbin.lt.1).or.(ixbin.gt.xbins(iana))) goto 120
c     if at least one particle of event contributes then count event
                    ievflag=1
                    encntr(iana,igrp)=encntr(iana,igrp)+1
                    tot_wt=tot_wt+aweight

c     add weights, modified by Shanshan
                    if(flag_wt.ne.0) then
                       zxygrd(iana,ixbin,igrp)=zxygrd(iana,ixbin,igrp)
     &                      + azval*aweight
                       dnsgrd1(iana,ixbin,igrp)=dnsgrd1(iana,ixbin,igrp)
     &                      + aweight
                       dnsgrd(iana,ixbin,igrp)=dnsgrd(iana,ixbin,igrp)
     &                      + 1
                       errgrd(iana,ixbin,igrp)=errgrd(iana,ixbin,igrp)
     &                      + azval*azval*aweight
                    else
                       zxygrd(iana,ixbin,igrp)=zxygrd(iana,ixbin,igrp)
     &                      + azval
                       dnsgrd(iana,ixbin,igrp)=dnsgrd(iana,ixbin,igrp)
     &                      + 1
                       errgrd(iana,ixbin,igrp)=errgrd(iana,ixbin,igrp)
     &                      + azval*azval
                    endif
 120             continue

c    end of particle loop
                 if(ievflag.eq.1) evcntr(iana,igrp)=evcntr(iana,igrp)+1 
              endif
 110       continue
c          end of grp loop
           elseif(yobs(iana).ne.0) then
c          make a 2d-grid
              igrp=zgrp(iana)
c             if eventobs. then do histo without particle loop
              if((xobs(iana).gt.100).and.(yobs(iana).gt.100)) then
c             first check, if both obs. refer to events
c             dn/dobs has different code for event-obs
              if(zobs(iana).eq.1) zobs(iana)=101
                 if(zobs(iana).lt.100) then
                    write(6,*)'warning: no mixing of event-
     & observables with particle observables!!'
                    goto 100
                 endif
               axval=evalue(xobs(iana),igrp)
               ayval=evalue(yobs(iana),igrp)
               inf=normflag(iana)
c so far no normalization for event-like observables
c               if(inf.ne.0) then
c                  azval=evalue(zobs(iana),igrp)*
c     &                  normvalue(xobs(iana),igrp,inf)*
c     &                  normvalue(yobs(iana),igrp,inf)
c               else
                  azval=evalue(zobs(iana),igrp)
c               endif
               ixbin=int((axval-xmin(iana))/dxbin(iana)+1.5)
               iybin=int((ayval-ymin(iana))/dybin(iana)+1.5)
               if((ixbin.lt.1).or.(ixbin.gt.xbins(iana))) goto 130
               if((iybin.lt.1).or.(iybin.gt.ybins(iana))) goto 130
               evcntr(iana,igrp)=evcntr(iana,igrp)+1
               encntr(iana,igrp)=encntr(iana,igrp)+1
              zxygrd(iana,ixbin,iybin)=zxygrd(iana,ixbin,iybin)+azval
               dnsgrd(iana,ixbin,iybin)=dnsgrd(iana,ixbin,iybin)+1
               errgrd(iana,ixbin,iybin)=errgrd(iana,ixbin,iybin)+
     &              azval*azval
c             now x is eventobs and y is particleobs, zobs must be particleobs
              elseif((xobs(iana).gt.100).and.(yobs(iana).lt.100))then
c                zobs has to be particleobs:
                 if(zobs(iana).gt.100) then
                    write(6,*)'wrong zobs for x- and yobs!!'
                    goto 100
                 endif
               axval=evalue(xobs(iana),igrp)
               ixbin=int((axval-xmin(iana))/dxbin(iana)+1.5)
               if((ixbin.lt.1).or.(ixbin.gt.xbins(iana))) goto 130
               ievflag=0
               inf=normflag(iana)
               call obsvalue(yobs(iana),npart,vavaly)
               call obsvalue(zobs(iana),npart,vavalz)
               if(inf.ne.0) then
                  call normvalue(yobs(iana),npart,inf,vavalny)
               endif
               do 160 igpart=1,nopart(igrp)
c               get particle-pointer
                idpart=grpart(igrp,igpart)
c               get y and z values of observble
                ayval=vavaly(idpart)
                azval=vavalz(idpart)
                if(inf.ne.0) then
                  azval=azval*vavalny(idpart)
                endif
c               check whether particle is in interval
                iybin=int((ayval-ymin(iana))/dybin(iana)+1.5)
                if((iybin.lt.1).or.(iybin.gt.ybins(iana))) goto 160
                ievflag=1
                encntr(iana,igrp)=encntr(iana,igrp)+1
              zxygrd(iana,ixbin,iybin)=zxygrd(iana,ixbin,iybin)+azval
                dnsgrd(iana,ixbin,iybin)=dnsgrd(iana,ixbin,iybin)+1
                errgrd(iana,ixbin,iybin)=errgrd(iana,ixbin,iybin)+
     &              azval*azval
 160           continue
c              end of particle loop
               if(ievflag.eq.1) evcntr(iana,igrp)=evcntr(iana,igrp)+1
c             now y is eventobs and x is particleobs, zobs must be particleobs
              elseif((yobs(iana).gt.100).and.(xobs(iana).lt.100))then
c                zobs has to be particleobs:
                 if(zobs(iana).gt.100) then
                    write(6,*)'wrong zobs for x- and yobs!!'
                    goto 100
                 endif
               ayval=evalue(yobs(iana),igrp)
               iybin=int((ayval-ymin(iana))/dybin(iana)+1.5)
               if((iybin.lt.1).or.(iybin.gt.ybins(iana))) goto 130
               ievflag=0
               inf=normflag(iana)
               call obsvalue(xobs(iana),npart,vavalx)
               call obsvalue(zobs(iana),npart,vavalz)
               if(inf.ne.0) then
                  call normvalue(xobs(iana),npart,inf,vavalnx)
               endif
               do 170 igpart=1,nopart(igrp)
c               get particle-pointer
                idpart=grpart(igrp,igpart)
c               get x and z values of observble
                axval=vavalx(idpart)
                azval=vavalz(idpart)
                if(inf.ne.0) then
                  azval=azval*vavalnx(idpart)
                endif
c               check whether particle is in interval
                ixbin=int((axval-xmin(iana))/dxbin(iana)+1.5)
                if((ixbin.lt.1).or.(ixbin.gt.xbins(iana))) goto 170
                ievflag=1
                encntr(iana,igrp)=encntr(iana,igrp)+1
              zxygrd(iana,ixbin,iybin)=zxygrd(iana,ixbin,iybin)+azval
                dnsgrd(iana,ixbin,iybin)=dnsgrd(iana,ixbin,iybin)+1
                errgrd(iana,ixbin,iybin)=errgrd(iana,ixbin,iybin)+
     &              azval*azval
 170           continue
c              end of particle loop
               if(ievflag.eq.1) evcntr(iana,igrp)=evcntr(iana,igrp)+1
              else 
c             particle observables:
c             loop over all particles in group
              ievflag=0
              inf=normflag(iana)
              call obsvalue(xobs(iana),npart,vavalx)
              call obsvalue(yobs(iana),npart,vavaly)
              call obsvalue(zobs(iana),npart,vavalz)
              if(inf.ne.0) then
                 call normvalue(xobs(iana),npart,inf,vavalnx)
                 call normvalue(yobs(iana),npart,inf,vavalny)
              endif
              do 140 igpart=1,nopart(igrp)
c              get particle-pointer
               idpart=grpart(igrp,igpart)
c              get x and z values of observble
               axval=vavalx(idpart)
               ayval=vavaly(idpart)
               azval=vavalz(idpart)
               if(inf.ne.0) then
                  azval=azval*vavalnx(idpart)*vavalny(idpart)
               endif
c              check whether particle is in interval
               ixbin=int((axval-xmin(iana))/dxbin(iana)+1.5)
               iybin=int((ayval-ymin(iana))/dybin(iana)+1.5)
               if((ixbin.lt.1).or.(ixbin.gt.xbins(iana))) goto 140
               if((iybin.lt.1).or.(iybin.gt.ybins(iana))) goto 140
               ievflag=1
               encntr(iana,igrp)=encntr(iana,igrp)+1
              zxygrd(iana,ixbin,iybin)=zxygrd(iana,ixbin,iybin)+azval
               dnsgrd(iana,ixbin,iybin)=dnsgrd(iana,ixbin,iybin)+1
               errgrd(iana,ixbin,iybin)=errgrd(iana,ixbin,iybin)+
     &              azval*azval
 140          continue
c             end of particle loop
              if(ievflag.eq.1) evcntr(iana,igrp)=evcntr(iana,igrp)+1
              endif
 130       continue
c          jump label
         endif
c        end of 1d/2d histo selection
 100  continue
c     end of zaxis loop
c output of scatterplots
      if(splts.eq.0)goto 31
      do 850 isp=1,splts
         write(envname,851)'ftn',scatunit(isp)
 851     format(a3,i2)
         call getenv(envname,fnam)
         iun=scatunit(isp)
         open(unit=iun,file=fnam,status='new',form='formatted')
c        now write file-header
         write(iun,501)cflag
         write(iun,502)cflag
         write(iun,503)cflag
         write(iun,504)cflag,header
         write(iun,503)cflag
         write(iun,505)cflag,totev,model_tag,version_tag
         write(iun,503)cflag
         write(iun,506)cflag
         write(iun,89)cflag,refsys,tstep
 501     format(a2,'************* BMS_AP 1.0 beta ***********')
 502     format(a2,'(c) S. A. Bass')
 503     format(a2)
 504     format(a2,a78)
 505     format(a2,i8,' events calulated by ',2a8)
 506     format(a2,'HICAP parameters:')
 89      format(a2,' Refsys:',I2,' Timestep:',I3)
c        print out selected particle-cuts
         if(cuts.eq.0) then
            write(iun,507)cflag
 507        format(a2,'no particle-cuts selected')
         else
            write(iun,508)cflag
 508        format(a2,'selected particle-cuts:')
            do 101 i=1,cuts
               if((scaty(isp).ne.0).and.((cutgrp(i).ne.0).or.
     &            (scatz(isp).ne.cutgrp(i)))) goto 101
               write(iun,201)cflag,cutobs(i),cutmin(i),cutmax(i),
     &                       cutgrp(i)
 201           format(a2,'Obs:',I3,' Min:',F8.2,' Max:',F8.2,
     &                ' Grp:',I3)

 101        continue
         endif
c        now print out event-cuts
         if(ecuts.eq.0) then
            write(iun,509)cflag
 509        format(a2,'no event-cuts selected')
         else
            write(iun,510)cflag
 510        format(a2,'selected event-cuts:')
            write(iun,511)cflag,anev
 511        format(a2,'# of events in event-cuts: ',I8)
            do 102 i=1,ecuts
               write(iun,201)cflag,ecutobs(i),ecutmin(i),ecutmax(i),
     &                       ecutgrp(i)
 102        continue
         endif
         write(iun,503)cflag
         write(iun,512)cflag
 512     format(a2,'selected analysis:')
         call obstring(scatx(isp),ostr)
         write(iun,203)cflag,ostr
 203     format(a2,'X-Obs: ',a15)
         call obstring(scaty(isp),ostr)
         write(iun,204)cflag,ostr
 204     format(a2,'Y-Obs: ',a15)
         call obstring(scatz(isp),ostr)
         write(iun,205)cflag,ostr
 205     format(a2,'Z-Obs: ',a15)
         ig=scatgrp(isp)
         write(iun,515)cflag,ig
 515     format(a2,'selected analysis-group: ',i3)
         write(iun,516)cflag,ig
 516     format(a2,'group ',I3,' is defined as:')
         write(iun,206)cflag,ig,angrp(ig,29),
     *        (angrp(ig,j),j=1,angrp(ig,30))
 206     format(a2,'Grp:',I2,' Flag:',I3,' Itypes:',28I3)
         if((scatx(igsp).gt.99).and.(scaty(isp).gt.99).and.
     *        (scatz(isp).gt.99)) then
            xx=evalue(scatx(isp),ig)
            yy=evalue(scaty(isp),ig)
            zz=evalue(scatz(isp),ig)
            write(iun,517)xx,yy,zz
 517        format(3E10.3)
         else
            if(scatx(isp).lt.99) then 
               call obsvalue(scatx(isp),npart,vavalx)
            else
               xx=evalue(scatx(isp),ig)
            endif
            if(scaty(isp).lt.99) then 
               call obsvalue(scaty(isp),npart,vavaly)
            else
               yy=evalue(scaty(isp),ig)
            endif
            if(scatz(isp).lt.99) then 
               call obsvalue(scatz(isp),npart,vavalz)
            else
               zz=evalue(scatz(isp),ig)
            endif
            do 518 ipn=1,nopart(ig)
               idpart=grpart(ig,ipn)
               if(scatx(isp).lt.99) then
                  xx=vavalx(idpart)
               endif
               if(scaty(isp).lt.99) then
                  yy=vavaly(idpart)
               endif
               if(scatz(isp).lt.99) then
                  zz=vavalz(idpart)
               endif
               write(iun,517)xx,yy,zz
 518        continue
         endif
 850  continue
      goto 31
 31   continue

 933  continue ! end reading in evsamp blocks


      if(tstep.lt.nsteps) then
         do 32 iii=tstep*evsamp+1,nsteps*evsamp
cdebug
c           write(6,*) ' hello ',iii,nsteps
            istep=iii
            call gettimestep
            if(iret.eq.0) then
               write(6,*) '*** end of event-input ***'
               call output
               stop
            elseif(iret.eq.-1) then
               write(6,*) '** error during event-input in uqmdnxtev **'
               write(6,*) '** last event was scipped, terminating...**'
               call output
               stop
            end if

 32      continue
      endif

      goto 3
                                                                    
      end

