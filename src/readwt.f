      subroutine read_wt_table(wt_flag)

      implicit none

      include 'ucoms.f'

      double precision dummy_pT,wtsum,dsigmaOVERdy
      integer wt_flag,i
      character*77 file10

      wtsum=0d0

      file10='      '
      call getenv('ftn10',file10)
      if (file10(1:4).ne.'    ') then
         OPEN(UNIT=10,FILE=file10,STATUS='OLD',FORM='FORMATTED')
      else
         OPEN(UNIT=10,FILE='pT*_wt.dat',STATUS='OLD',FORM='FORMATTED')
      endif

      i=0
      do while(i.le.wt_num)
         read(unit=10,fmt=*,err=2198,end=2197) dummy_pT,pT_wt_table(i)
         pT_wt_table(i)=dummy_pT*pT_wt_table(i)
         wtsum=wtsum+pT_wt_table(i)
         i=i+1
      enddo

c this was a d\sigma/d^2p_T/d\eta table
c we change it to dN/dp table normalized to 1
c by default, the bin size is 0.5~GeV, might be changed due to request

c calculate cross section of HQ production
      dsigmaOVERdy=wtsum*2d0*PI*wt_int
      write(6,*) "dsigma/dy (y=0): ",dsigmaOVERdy

      i=0
      do while(i.le.wt_num)
         pT_wt_table(i)=pT_wt_table(i)/wtsum/wt_int
         i=i+1
      enddo

      wt_flag=0
      write(6,*) "weight table has been read in successfully :)"
      return

 2197 continue
      wt_flag=1
      write(6,*) 'ERROR: EOF reached in weight table'
      write(6,*) 'terminating ...'
      return

 2198  continue
       wt_flag=-1
       write(6,*) 'READ-ERROR in weight table'
       write(6,*) 'terminating ...'
       return

       end





! this is a subroutine read in the rotate angle
       subroutine read_rotation_angle

       implicit none


       include 'ucoms.f'
       character*77 file24

       file24='    '
       call getenv('ftn24', file24)
       if (file24(1:4) .ne. '    ') then
         OPEN(UNIT=24, FILE=file24, STATUS='OLD', FORM='FORMATTED')
         READ(UNIT=24, fmt=*, err=2200, end=2200) cos_pp, sin_pp
       else
         WRITE(6,*) "no inputs to rotate ..."
       endif
 
       write(6,*) "rotation angle set in successfully :)"
       return

 2200  continue
       write(6,*) "no rotation this time..."
       return 

       end
       

       
