C *************************************
C *       PROGRAMME  placement_H2O.f  *
C *************************************
C * Calcul de la distribution des tailles de pores

C*
C* Usage : placement_H2O -filepos1 'fichier des positions'  -filepos2 'Fichier H2O'
C*

      module routine
        integer, parameter :: e= 8
      end module

      PROGRAM DISTRI_ANGLE
      use routine
      real(kind=e) :: rbohr, zero, pi
      integer :: im0
      PARAMETER (rbohr=0.529177d0)
      PARAMETER (pi=3.14159265d0)
      PARAMETER (zero=0.0d0)
      PARAMETER (IM0=150000)
      character(len=100), dimension(:), allocatable :: buffl
      character(len=116) :: aa2
      character(len=8) :: chaine8
      character(len=2), dimension(9) :: ctyp
      character(len=2) :: c2
      character(len=60) :: nomfile, nomfilepos1, nomfilepos2
      character(len=180) :: buff
      real(kind=e), dimension(IM0,3) :: xp1, vp1, fp1
      real(kind=e), dimension(IM0,3) :: xp2, vp2, fp2
      real(kind=e) :: zl1, zl2, zl3, x1, x2, x3, v1, v2, v3, f1, f2, f3
      real(kind=e) :: ao1, ao2, ao3, ah11, ah12, ah13, ah21, ah22, ah23
      real(kind=e) :: zls1, zls2, zls3
      real(kind=e) :: dro, drh1, drh2
      real(kind=e) :: xo, yo, zo, xh1, yh1, zh1, xh2, yh2, zh2
      real(kind=e) :: dxo, dyo, dzo, dxh1, dyh1, dzh1, dxh2, dyh2, dzh2
      real(kind=e) :: rcut1, xx1
      integer :: narg, istep, natom1, natom2, iform, indice
      integer :: im1, im2, ipoint
      integer :: i, iti, iconfig, ndis, imcon, ncut
      integer :: ifile, ifilepos1, ifilepos2
      integer, dimension(9) :: na
      integer, dimension(IM0) :: ityp1, ityp2
      rcut1=0.6d0    ! cut_off_changing_location 
      rcut12=rcut1*rcut1

      narg=iargc()
      allocate(buffl(narg))  !Tableau des options
      do i=1,narg
      call getarg(i,buff)
      buffl(i) = buff
      enddo
      ifilepos1 = 0
      ifilepos2 = 0

      do i=1,narg
      if (buffl(i) == "-filepos1") then
      ifilepos1=1
      nomfilepos1 = buffl(i+1)
      endif
      if (buffl(i) == "-filepos2") then
      ifilepos2=1
      nomfilepos2 = buffl(i+1)
      endif
      enddo
      deallocate(buffl)  !Tableau des options

      if (ifilepos1.eq.0) then
      write(6,*) 'Le premier fichier des positions n est pas defini'
      stop 
      endif
      if (ifilepos2.eq.0) then
      write(6,*) 'Le second fichier des positions n est pas defini'
      stop 
      endif

      write(6,'(24HFichier des positions : ,1X,A60)') nomfilepos1
!      write(6,*) 'Taille : ',ZL1,ZL2,ZL3 !TEST
      OPEN (60,file=nomfilepos1,form='formatted',iostat=iend)
      rewind 60
      write(6,'(24HFichier des positions : ,1X,A60)') nomfilepos2
      OPEN (70,file=nomfilepos2,form='formatted',iostat=iend)
      rewind 70

      CTYP(1)='Si'
      CTYP(2)='B '
      CTYP(3)='O '
      CTYP(4)='Na'
      CTYP(5)='Al'
      CTYP(6)='Ca'
      CTYP(7)='Zr'
      CTYP(8)='OW'
      CTYP(9)='HW'

      DO I=1,9 ! Initialisation du tableau NA()
      NA(I)=0
      ENDDO

      read (60,168,end=999) AA2
      do while (.true.)
      read (60,168,end=999) AA2
  168 format (A116)
        read(AA2,*) iform,imcon,natom1,istep
! Lecture des dimensions
          read (60,168,end=999) AA2
          read(AA2,*) ZL1,x2,x3
          ZLS1=ZL1/2.d0
          read (60,168,end=999) AA2
          read(AA2,*) x1,ZL2,x3
          ZLS2=ZL2/2.d0
          read (60,168,end=999) AA2
          read(AA2,*) x1,x2,ZL3
          ZLS3=ZL3/2.d0
! Fin de la lecture des dimensions
! Lecture des atomes
          do i=1,natom1
          read (60,168,end=999) AA2
          read(AA2,*) C2,indice
            do j=1,9
            if (C2.eq.CTYP(j)) then
            ityp1(indice)=j
            na(j)=na(j)+1
            endif
            enddo
          read (60,168,end=999) AA2
          read(AA2,*) x1,x2,x3
          XP1(indice,1)=x1
          XP1(indice,2)=x2
          XP1(indice,3)=x3
            if (iform.ge.1) then
            read (60,168,end=999) AA2
            read(AA2,*) v1,v2,v3
            VP1(indice,1)=v1
            VP1(indice,2)=v2
            VP1(indice,3)=v3
              if (iform.ge.2) then
              read (60,168,end=999) AA2
              read(AA2,*) f1,f2,f3
              FP1(indice,1)=f1
              FP1(indice,2)=f2
              FP1(indice,3)=f3
              endif
            endif
          enddo
      enddo
          IM1 = NA(1)+NA(2)+NA(3)+NA(4)+NA(5)+NA(6)+NA(7)
! Fin de la lecture des atomes
  999 CONTINUE

! Lecture de la seconde configuration
      read (70,168,end=998) AA2
      do while (.true.)
      read (70,168,end=998) AA2
        read(AA2,*) iform,imcon,natom2,xx1
! Lecture des dimensions
          read (70,168,end=998) AA2
          read(AA2,*) ZL1,x2,x3
          ZLS1=ZL1/2.d0
          read (70,168,end=998) AA2
          read(AA2,*) x1,ZL2,x3
          ZLS2=ZL2/2.d0
          read (70,168,end=998) AA2
          read(AA2,*) x1,x2,ZL3
          ZLS3=ZL3/2.d0
! Fin de la lecture des dimensions
! Lecture des atomes
          do i=1,natom2
          read (70,168,end=998) AA2
          read(AA2,*) C2,indice
            do j=1,9
            if (C2.eq.CTYP(j)) then
            ityp2(indice)=j
            na(j)=na(j)+1
            endif
            enddo
          read (70,168,end=998) AA2
          read(AA2,*) x1,x2,x3
          XP2(indice,1)=x1
          XP2(indice,2)=x2
          XP2(indice,3)=x3
            if (iform.ge.1) then
            read (70,168,end=998) AA2
            read(AA2,*) v1,v2,v3
            VP2(indice,1)=v1
            VP2(indice,2)=v2
            VP2(indice,3)=v3
              if (iform.ge.2) then
              read (70,168,end=998) AA2
              read(AA2,*) f1,f2,f3
              FP2(indice,1)=f1
              FP2(indice,2)=f2
              FP2(indice,3)=f3
              endif
            endif
          enddo
      enddo
      IM2 = NA(8)+NA(9)
  998 CONTINUE
      OPEN (70,file='CONFIG_GLASS_WATER',form='formatted')
      WRITE(70,*) 'H2O_CONFIG File'
      WRITE(70,*) '        0         3'
      write(70,'(3f20.15)') ZL1,ZERO,ZERO
      write(70,'(3f20.15)') ZERO,ZL2,ZERO
      write(70,'(3f20.15)') ZERO,ZERO,ZL3

! Recherche des H2O (configuration 2) loin des atomes (configuration 1)
      ncut=0
      do i=1,natom2,3
         xo=xp2(i,1)
         yo=xp2(i,2)
         zo=xp2(i,3)
         ao1=-DSIGN(ZL1,xo)
         ao2=-DSIGN(ZL2,yo)
         ao3=-DSIGN(ZL3,zo)
         xh1=xp2(i+1,1)
         yh1=xp2(i+1,2)
         zh1=xp2(i+1,3)
         ah11=-DSIGN(ZL1,xh1)
         ah12=-DSIGN(ZL2,yh1)
         ah13=-DSIGN(ZL3,zh1)
         xh2=xp2(i+2,1)
         yh2=xp2(i+2,2)
         zh2=xp2(i+2,3)
         ah21=-DSIGN(ZL1,xh2)
         ah22=-DSIGN(ZL2,yh2)
         ah23=-DSIGN(ZL3,zh2)

         ipoint=1
         do j=1,natom1 !Boucle sur les atomes (configuration 1)
         if (ipoint.eq.1) then
         dxo=xo-xp1(j,1)
         dyo=yo-xp1(j,2)
         dzo=zo-xp1(j,3)
         if(dabs(dxo).gt.zls1) dxo=dxo+ao1 
         if(dabs(dyo).gt.zls2) dyo=dyo+ao2 
         if(dabs(dzo).gt.zls3) dzo=dzo+ao3 
         dro=dxo*dxo+dyo*dyo+dzo*dzo

         dxh1=xh1-xp1(j,1)
         dyh1=yh1-xp1(j,2)
         dzh1=zh1-xp1(j,3)
         if(dabs(dxh1).gt.zls1) dxh1=dxh1+ah11 
         if(dabs(dyh1).gt.zls2) dyh1=dyh1+ah12 
         if(dabs(dzh1).gt.zls3) dzh1=dzh1+ah13 
         drh1=dxh1*dxh1+dyh1*dyh1+dzh1*dzh1

         dxh2=xh2-xp1(j,1)
         dyh2=yh2-xp1(j,2)
         dzh2=zh2-xp1(j,3)
         if(dabs(dxh2).gt.zls1) dxh2=dxh2+ah21 
         if(dabs(dyh2).gt.zls2) dyh2=dyh2+ah22 
         if(dabs(dzh2).gt.zls3) dzh2=dzh2+ah23 
         drh2=dxh2*dxh2+dyh2*dyh2+dzh2*dzh2

         if (dro.lt.rcut12.or.drh1.lt.rcut12.or.drh2.lt.rcut12) then
         ipoint=0
         endif
         endif
         enddo
      if (ipoint.eq.1) then
      ncut = ncut +1
      write(70,'(a2,I5)') CTYP(ityp2(i)),ncut
      write(70,'(3f20.15)') xo,yo,zo
      ncut = ncut +1
      write(70,'(a2,I5)') CTYP(ityp2(i+1)),ncut
      write(70,'(3f20.15)') xh1,yh1,zh1
      ncut = ncut +1
      write(70,'(a2,I5)') CTYP(ityp2(i+2)),ncut
      write(70,'(3f20.15)') xh2,yh2,zh2
      endif

      enddo 


         do j=1,natom1 !Boucle sur les atomes (configuration 1)
         write(70, '(a2,I5)') CTYP(ityp1(j)),j+ncut
         write(70, '(3f20.15)') xp1(j,1),xp1(j,2),xp1(j,3)
         enddo
      write(6,*) 'Number of selected H2O:',ncut/3


      
      STOP
      END

