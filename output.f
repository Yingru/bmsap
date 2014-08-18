      subroutine output
c     this subroutine handles the output of histos/grids in hicap
c
c     author   : Steffen A. Bass
c     date     : 04.01.96
c     revision : 0.9 beta
c

      IMPLICIT REAL*8(A-H,O-Z)

      character*5 envname
      character*15 ostr
      character*80 fnam
      real outy(60),outx
      include 'ucoms.f'
      include 'hicoms.f'
c
      bmean=bmean/totev
c     first make a loop over all possible outputs
      do 10 iana=1,zaxis
c        then get the name of the output-file and open it
         write(envname,88)'ftn',zunit(iana)
 88      format(a3,i2)
         call getenv(envname,fnam)
         if(fnam.eq.' ') then
            write(fnam,87)'bms_ap.',zunit(iana)
 87         format(a7,i2)
         endif            
         iun=zunit(iana)
         open(unit=iun,file=fnam,status='unknown',form='formatted')
c        now write file-header
         write(iun,501)cflag
 501     format(a2,'************* BMS_AP 1.0 beta ***********')
         write(iun,503)cflag
 503     format(a2)
         write(iun,502)cflag
 502     format(a2,'(c) S. A. Bass 2002')
         write(iun,503)cflag
         write(iun,504)cflag,header
 504     format(a2,a78)
         write(iun,503)cflag
         write(iun,505)cflag,totev,model_tag,version_tag
 505     format(a2,i8,' events calulated by ',2a8)
         write(iun,503)cflag
c display the uqmd parameters
         write(iun,601)cflag
 601     format(a2,'model parameters:')
         write(iun,602)cflag,ap,zp,at,zt
 602     format(a2,1X,'Projectile: ',i5,' ',i5,'  Target: ',i5,' ',i5)
         write(iun,603)cflag,bmean
 603     format(a2,1X,'mean impact parameter (fm): ',1f7.2)
         write(iun,604)cflag,ecm
 604     format(a2,1X,' sqrt(s)(GeV): ',e10.4e2)
c
         write(iun,503)cflag
         write(iun,506)cflag
 506     format(a2,'BMS_AP parameters:')
         write(iun,89)cflag,refsys,tstep,logscale,addition
 89      format(a2,' Refsys:',I2,' Timestep:',I3,' Logscale:',I2,
     &          ' Addition:',I2)
         write(iun,90)cflag,flag_wt,wt_num,wt_int,wt_Tab_min
 90      format(a2,' Weight Type:',I2,' wt_num:',I4,' wt_int:',F6.3,
     &       ' min_pT:',F6.3)
c        print out selected particle-cuts
         if(cuts.eq.0) then
            write(iun,507)cflag
 507        format(a2,'no particle-cuts selected')
         else
            write(iun,508)cflag
 508        format(a2,'selected particle-cuts:')
            do 101 i=1,cuts
               if((cutgrp(i).ne.0).and.(yobs(iana).ne.0).and.
     &            (zgrp(iana).ne.cutgrp(i))) goto 101
               call obstring(cutobs(i),ostr)
               write(iun,201)cflag,ostr,cutmin(i),cutmax(i),
     &                       cutgrp(i)
 201           format(a2,'Obs:',a15,'  Min:',d12.4,'  Max:',d12.4,
     &                '  Grp:',I3)

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
               call obstring(ecutobs(i),ostr)
               write(iun,201)cflag,ostr,ecutmin(i),ecutmax(i),
     &                       ecutgrp(i)
 102        continue
         endif
         write(iun,503)cflag
         write(iun,512)cflag
 512     format(a2,'selected analysis:')
c
         if(yobs(iana).eq.0) then
c        this one is a 1D-histo
c        first calculate errorbars
c        sigma^2=<x^2>-<x>^2
          do 20 ixb=1,xbins(iana)
            do 20 ig=1,grpsel
              if(flag_wt.ne.0) ddns1=dnsgrd1(iana,ixb,ig)
              ddns=dnsgrd(iana,ixb,ig)+eps
              if((zobs(iana).eq.1).or.(zobs(iana).eq.11)
     &             .or.(zobs(iana).eq.41).or.(zobs(iana).eq.42)
     &             .or.(zobs(iana).eq.101)) then
c normalize to particle number in each event
c                 if(flag_wt.ne.0) then
c                    zxygrd(iana,ixb,ig)=zxygrd(iana,ixb,ig)/
c     &                   (tot_wt*dxbin(iana))
c                    errgrd(iana,ixb,ig)=zxygrd(iana,ixb,ig)/sqrt(ddns)
c                 else
                    zxygrd(iana,ixb,ig)=zxygrd(iana,ixb,ig)/
     &                   (anev*dxbin(iana))
                    errgrd(iana,ixb,ig)=zxygrd(iana,ixb,ig)/sqrt(ddns)
c                 endif
              else
                 if(flag_wt.ne.0) then
                    zxygrd(iana,ixb,ig)=zxygrd(iana,ixb,ig)/ddns1
                    xsqm=errgrd(iana,ixb,ig)/ddns1
                    xmsq=zxygrd(iana,ixb,ig)
                    xmsq=xmsq*xmsq
                    errgrd(iana,ixb,ig)=sqrt(dabs((xsqm-xmsq)/(ddns-1)))
                 else
                    zxygrd(iana,ixb,ig)=zxygrd(iana,ixb,ig)/ddns
                    xsqm=errgrd(iana,ixb,ig)/ddns
                    xmsq=zxygrd(iana,ixb,ig)
                    xmsq=xmsq*xmsq
                    errgrd(iana,ixb,ig)=sqrt(dabs((xsqm-xmsq)/(ddns-1)))
                 endif
              endif

 20       continue

c         now prepare output:
            call obstring(xobs(iana),ostr)
            write(iun,202)cflag,ostr,xmin(iana),xmax(iana),
     &                    xbins(iana)
 202        format(a2,'X-Obs: ',a15,' Min:',F8.2,' Max:',F8.2,
     &                ' Bins:',I3)
            call obstring(zobs(iana),ostr)
            if(normflag(iana).gt.0) then
               write(iun,401)cflag,ostr
            elseif(normflag(iana).lt.0) then
               write(iun,402)cflag,ostr
            else
               write(iun,203)cflag,ostr
            endif
 401        format(a2,'Z-Obs: ',a15,' kanonical Normalization')
 402        format(a2,'Z-Obs: ',a15,' alternate Normalization')
 203        format(a2,'Z-Obs: ',a15)
            write(iun,513)cflag
 513        format(a2,'all analysis-groups selected')
            write(iun,503)cflag
            write(iun,514)cflag
 514        format(a2,'groups are defined as:')
            if(allpart) then
               write(iun,901)cflag
 901           format(a2,'all particles in one group:')
               write(iun,902)cflag
 902           format(a2,'-> selection via cuts')
            else
               do 103 i=1,grpsel
                  write(iun,206)cflag,i,angrp(i,29),
     *                 (angrp(i,j),j=1,angrp(i,30))
 103           continue
 206           format(a2,'Grp:',I2,' Flag:',I3,' Itypes:',28I7)
            endif
            write(iun,503)cflag
c           now prepare output of histos: x, y1,dy1,y2,dy2...
            write(iun,301)cflag,(ig,ig=1,grpsel)
 301        format(a2,'Grp:        ',30(I3,'      Errorbar:  '))
            write(iun,302)cflag,(encntr(iana,ig),ig=1,grpsel)
 302        format(a2,'Entries   ',30(I8,12X))
            write(iun,303)cflag,(evcntr(iana,ig),ig=1,grpsel)
 303        format(a2,'Events    ',30(I8,12X))
            write(iun,520)cflag,tot_wt/anev
 520        format(a2,'tot_wt    ',E10.3)

           outx=xmin(iana)
c +(dxbin(iana)/2)
            igg=2*grpsel
            do 70 ix=1,xbins(iana)
c             now fill vector for one output line
              igy=0
              do 60 ig=1,igg,2
                igy=igy+1
c  convert the error into logscale if necessary
                if(logscale.eq.1) then
                  outy(ig)=log(zxygrd(iana,ix,igy))
                  outy(ig+1)=errgrd(iana,ix,igy)/zxygrd(iana,ix,igy)
                else
                  outy(ig)=zxygrd(iana,ix,igy)
                  outy(ig+1)=errgrd(iana,ix,igy)
                endif
c  addition is for the convenience of error analysis
                if(addition.eq.1) then
                  outy(ig)=outy(ig)+outy(ig+1)
                elseif(addition.eq.(-1)) then
                  outy(ig)=outy(ig)-outy(ig+1)
                endif
 60           continue
c             write line of output
              write(iun,304) outx,(outy(ig),ig=1,igg)
 304          format(2X,61(1X,E10.3))
              outx=outx+dxbin(iana)
 70         continue
         else
c        here comes a 2D grid
c        first calculate errorbars
c        sigma^2=<x^2>-<x>^2
          do 30 ixb=1,xbins(iana)
            do 30 iyb=1,ybins(iana)
              ddns=dnsgrd(iana,ixb,iyb)+1e-8
              if((zobs(iana).eq.1).or.(zobs(iana).eq.101)) then
c               zxygrd(iana,ixb,iyb)=zxygrd(iana,ixb,iyb)/
c     &          (evcntr(iana,zgrp(iana))*dxbin(iana)*dybin(iana))
               zxygrd(iana,ixb,iyb)=zxygrd(iana,ixb,iyb)/
     &                              (anev*dxbin(iana)*dybin(iana))
               errgrd(iana,ixb,iyb)=zxygrd(iana,ixb,iyb)/sqrt(ddns)
              else
               xsqm=errgrd(iana,ixb,iyb)/ddns
               xmsq=zxygrd(iana,ixb,iyb)/ddns
               xmsq=xmsq*xmsq
               zxygrd(iana,ixb,iyb)=zxygrd(iana,ixb,iyb)/ddns
               errgrd(iana,ixb,iyb)=sqrt(abs((xsqm-xmsq)/(ddns-1)))
              endif
 30       continue
c         now prepare output:
            call obstring(xobs(iana),ostr)
            write(iun,202)cflag,ostr,xmin(iana),xmax(iana),
     &                    xbins(iana)
            call obstring(yobs(iana),ostr)
            write(iun,204)cflag,ostr,ymin(iana),ymax(iana),
     &                    ybins(iana)
 204        format(a2,'Y-Obs: ',a15,' Min:',F8.2,' Max:',F8.2,
     &                ' Bins:',I3)
            call obstring(zobs(iana),ostr)
            if(normflag(iana).gt.0) then
               write(iun,401)cflag,ostr
            elseif(normflag(iana).lt.0) then
               write(iun,402)cflag,ostr
            else
               write(iun,203)cflag,ostr
            endif
            ig=zgrp(iana)
            write(iun,515)cflag,ig
 515        format(a2,'selected analysis-group: ',i3)
            write(iun,516)cflag,ig
 516        format(a2,'group ',I3,' is defined as:')
            write(iun,206)cflag,ig,angrp(ig,29),
     *       (angrp(ig,j),j=1,angrp(ig,30))
            write(iun,517)cflag,encntr(iana,ig)
 517        format(a2,'# of entries: ',I9)
            write(iun,518)cflag,evcntr(iana,ig)
 518        format(a2,'# of events : ',I8)

c           now print out the grid
            do 40 jy=1,ybins(iana)
c  convert the error into logscale if necessary
              if(logscale.eq.1) then
                write(iun,91) (log(zxygrd(iana,jx,jy)),jx=1,xbins(iana))
              else
                write(iun,91) (zxygrd(iana,jx,jy),jx=1,xbins(iana))
              endif
 91           format(100(2X,E10.3))
 40         continue
c           now the errorgrid has to be printed out:
c           right now the errorgrid is printed as a comment to the normal grid
c           otherwise one also could specify another output-file
c           and use the grid as a contour-overlay
            write(iun,503)cflag
            write(iun,519)cflag
 519        format(a2,'Grid with error-values:')
            do 45 jy=1,ybins(iana)
c  convert the error into logscale if necessary
            if(logscale.eq.1) then
              write(iun,95)cflag,(errgrd(iana,jx,jy)/zxygrd(iana,jx,jy),
     *         jx=1,xbins(iana))
            else
               write(iun,95)cflag,(errgrd(iana,jx,jy),jx=1,xbins(iana))
            endif
 95         format(a2,100(1X,E10.3))
 45         continue
         endif
 10   continue
      write(6,333)  logscale
 333   format(1X,'logscale in output file:', I2)
                                                             
      return
      end
