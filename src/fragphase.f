      subroutine frag
c     this subroutine determines fragments according to a coalescence model
c     
c     author  : Christiane Lesny
c     date    : 20.07.94
c     
c     based on Clusters, (c) by Jens Konopka /Univ. Frankfurt
c                            and Steffen Bass /GSI Darmstadt
c
      IMPLICIT REAL*8(A-H,O-Z)

      include 'stndrdb.f'
      logical inclus
      dimension xn(nparmx),yn(nparmx),zn(nparmx),                      
     &          pxn(nparmx),pyn(nparmx),pzn(nparmx),  
     &          emn(nparmx),itypn(nparmx),ichargen(nparmx),
     &          lastcln(nparmx),nclcntn(nparmx),rhomaxn(nparmx),
     &          ndeltn(nparmx),mdeltn(nparmx),nfrtimen(nparmx),
     &          iclarr(nparmx,nparmx),number(nparmx),
     &          kcount(nparmx),inclus(nparmx)
      parameter (pclus2=16900.)
      np=nnucl
c     itotpan is the new all-particle counter
      itotpan=0
c initalization, generate charge-vector
c
c     iparts(k) will later give the number of nucl. in a cluster -> mass
c     iclarr gives the linkes between particles
         do 80 k=1,np
            iparts(k) = 1
c     inclus(k) gibt die Teilchen an, die schon in einem Cluster enthalten
c     sind
            inclus(k)=.false.
            if ( ityp(k) .eq. 2 ) then
               icharge(k) = 1
            else
               icharge(k) = 0
            end if
            do 80 l=1,np
               iclarr(k,l) = 0
  80     continue

c *******************************************************************
c Schwerpunktcoalescence von Christiane Lesny
c momentum coalescence

c Setze das Quadrat von pmin auf unendlich
    5    pmin2=1E20
c Suche pmin = minimaler Impulsabstand zweier Nukleonen

         do 10 k=1,np
          if (inclus(k).eqv..true.) goto 10
          
          do 20 l=k+1,np
           if (inclus(l).eqv..true.) goto 20
 
             dp2=(px(l)-px(k))*(px(l)-px(k))+(py(l)-py(k))*(py(l)-py(k))
     +          +(pz(l)-pz(k))*(pz(l)-pz(k))
                       
             if (dp2.lt.pmin2) then
                pmin2=dp2
                lmin=l
                kmin=k
             endif

   20     continue
   10    continue
 
c Alle Cluster gefunden; der Rest sind Einzelteilchen->Ende der Routine
         
         if (pmin2.gt.(4.*pclus2)) goto 110

c Berechne Schwerpunkt der zwei Teilchen
      
         psptx=(px(lmin)+px(kmin))*0.5
         pspty=(py(lmin)+py(kmin))*0.5
         psptz=(pz(lmin)+pz(kmin))*0.5
       
c Wie viele = icountfrag und welche = kcount Teilchen bilden ein Cluster

   50    icountfrag = 0

c Schlage um Schwerpunkt Kreis und suche alle Teilchen im Kreis
  
         do 30 k=1,np
          if (inclus(k).eqv..true.) goto 30
     
          dppspt2 =(px(k)-psptx)*(px(k)-psptx)+(py(k)-pspty)
     +            *(py(k)-pspty)+(pz(k)-psptz)*(pz(k)-psptz)   
          
          if (dppspt2.le.pclus2) then
            icountfrag=icountfrag+1
            kcount(icountfrag)=k
          endif        

   30    continue

c Berechne neuen Schwerpunkt dieser Teilchen
  
         pnewx=0.
         pnewy=0.
         pnewz=0.

         do 70 k=1,icountfrag
        
           icluspart=kcount(k)

           pnewx=pnewx+px(icluspart)
           pnewy=pnewy+py(icluspart)
           pnewz=pnewz+pz(icluspart)
 
   70    continue

           pnewx=pnewx/dble(icountfrag)
           pnewy=pnewy/dble(icountfrag)
           pnewz=pnewz/dble(icountfrag)
           
c Wenn der neue Schwerpunkt /= alter Schwerpunkt 
c -> noch einmal neu clustern

           if ((abs(psptx-pnewx).gt.1E-4).or.(abs(pspty-pnewy).gt.1E-4)
     +         .or.(abs(psptz-pnewz).gt.1E-4)) then
              
              psptx=pnewx
              pspty=pnewy
              psptz=pnewz
  
              goto 50
           endif

c Wenn Cluster O.K. Links setzen
       
           do 60 k=1,icountfrag
       
               icluspart=kcount(k)
              if (k.eq.icountfrag) then
               icluspart2=kcount(1)
              else
               icluspart2=kcount(k+1)
              endif

                iclarr(icluspart,icluspart2)=1
                iclarr(icluspart2,icluspart)=1
                inclus(icluspart)=.true.
  60       continue

c Naechstes Cluster -> Beginn von vorne
           goto 5

  110      continue    
c Impulsraumkoaleszenz abgeschlossen
c Ortsraumkoaleszenz bei Freezoutzeit 
         
           do 500 k=1,np
             do 510 l=k+1,np

c Wenn die beiden Nukleonen im Ortsraum ein Cluster bilden, dann
              if (iclarr(k,l).eq.1) then
               
c wird das Nukleon mit der kleineren Freezeoutzeit auf einer Geraden
c zur Freezeoutzeit des zweiten Nukleons propagiert

               if (nfrtime(k) .gt. nfrtime(l)) then
                 call propag (nfrtime(k),nfrtime(l),xf(l),yf(l),zf(l),
     &                        pxf(l),pyf(l),pzf(l),xet2,yet2,zet2)
        
                 xet1=xf(k)
                 yet1=yf(k)
                 zet1=zf(k)
              
               else
                 call propag (nfrtime(l),nfrtime(k),xf(k),yf(k),zf(k),   
     &                        pxf(k),pyf(k),pzf(k),xet2,yet2,zet2)
     
                 xet1=xf(l)
                 yet1=yf(l)
                 zet1=zf(l)

               endif

c Abstand im Ortraum
              dret2=(xet1-xet2)*(xet1-xet2)+(yet1-yet2)*(yet1-yet2)
     &             +(zet1-zet2)*(zet1-zet2)

c Ist der Abstand gr"osser als 3 fm, l"ose die Verbindung wieder
              if (dret2 .gt.9.00) then
                  
                   iclarr(k,l)=0
                   iclarr(l,k)=0
              endif
              endif
 510        continue
 500       continue
  
c           do 67 k=1,np
c              write(6,*)(iclarr(k,l),l=1,np)
c 67        continue   
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
               xn(itotpan)=x(k)
               yn(itotpan)=y(k)
               zn(itotpan)=z(k)
               pxn(itotpan)=px(k)
               pyn(itotpan)=py(k)
               pzn(itotpan)=pz(k)
               emn(itotpan)=em(k)
               itypn(itotpan)=ityp(k)
               ichargen(itotpan)=icharge(k)
               lastcln(itotpan)=lastcl(k)
               nclcntn(itotpan)=nclcnt(k)
               rhomaxn(itotpan)=rhomax(k)
               ndeltn(itotpan)=ndelt(k)
               mdeltn(itotpan)=mdelt(k)
               nfrtimen(itotpan)=nfrtime(k)
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

               xcm = xcm + x(j)*em(j)
               ycm = ycm + y(j)*em(j)
               zcm = zcm + z(j)*em(j)
               pxcm = pxcm + px(j)
               pycm = pycm + py(j)
               pzcm = pzcm + pz(j)
               totem = totem + em(j)

               ntotcl = ntotcl + nclcnt(j)
               if(rhomax(j).ge.rhom) rhom=rhomax(j)
               if(nfrtime(j).ge.frti) then 
                 frti=nfrtime(j)
                 mfdelt=mdelt(j)
               endif
               npimi=npimi+mod(ndelt(j),100)
               npize=npize+mod(ndelt(j)/100,100)
               npipl=npipl+mod(ndelt(j)/10000,100)

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
            rhomaxn(itotpan)=rhom
            ndeltn(itotpan)=npimi+100*npize+10000*npipl
            mdeltn(itotpan)=mfdelt
            nfrtimen(itotpan)=frti

            
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
         do 200 k=(np+1),npart(1)
            itotpan=itotpan+1
               xn(itotpan)=x(k)
               yn(itotpan)=y(k)
               zn(itotpan)=z(k)
               pxn(itotpan)=px(k)
               pyn(itotpan)=py(k)
               pzn(itotpan)=pz(k)
               emn(itotpan)=em(k)
               itypn(itotpan)=ityp(k)
               if (ityp(k).eq.7) then
                  ichargen(itotpan)=-1
               elseif(ityp(k).eq.8) then
                  ichargen(itotpan)=0
               else
                  ichargen(itotpan)=1
               endif
               lastcln(itotpan)=lastcl(k)
               nclcntn(itotpan)=nclcnt(k)
               rhomaxn(itotpan)=rhomax(k)
               ndeltn(itotpan)=ndelt(k)
               mdeltn(itotpan)=mdelt(k)
               nfrtimen(itotpan)=nfrtime(k)
 200        continue
            
c rewrite particle-arrays
            do 300 k=1,itotpan
               x(k)=xn(k)
               y(k)=yn(k)
               z(k)=zn(k)
               px(k)=pxn(k)
               py(k)=pyn(k)
               pz(k)=pzn(k)
               em(k)=emn(k)
               ityp(k)=itypn(k)
               icharge(k)=ichargen(k)
               lastcl(k)=lastcln(k)
               nclcnt(k)=nclcntn(k)
               rhomax(k)=rhomaxn(k)
               ndelt(k)=ndeltn(k)
               mdelt(k)=mdeltn(k)
               nfrtime(k)=nfrtimen(k)
 300        continue
            npart(1)=itotpan
      return
      end



