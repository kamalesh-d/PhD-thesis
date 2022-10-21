! *****************************************
! *       PROGRAMME  config_aleatoire.f   *
! *****************************************
!  Placement aleatoire de molecules H2O
! To change. 1. PARAMETER(IM=*) * is total number of atoms
! 2. No of O, H
! 3. Size of the simulation box

      PROGRAM PREPAR
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(IM=1173)    !     total_number_of_water_molecules
      DIMENSION NA(10)
      DIMENSION XP(IM,3),VP(IM,3),FP(IM,3)
      DIMENSION CMASS(10),CHARGE(10)
      DIMENSION IOCCUP(IM)
	DIMENSION ITYP(IM)
	CHARACTER*6 filename
	CHARACTER*3 CTYP(10)
	data bohr / 0.5291772D0/
	data timau / 2.41889D-17/
	NA(1)=0    !Nombre de Si
	NA(2)=0     !Nombre de B
	NA(3)=391   !Nombre de O             number_of_oxygen_atoms
	NA(4)=782   !Nombre de Na (ici ce sont des H)        number_of_hydrogen_atoms
	NA(5)=0     !Nombre de Zr
	NA(6)=0     !Nombre de Al
	NA(7)=0     !Nombre de Ca
	NA(8)=0
	NA(9)=0
	NA(10)=0
        ZL_x=18.0095376377d0     !Arete de la boite de simulation  x_coordinate_to_be_modified
        ZL_y=18.0095376377d0        !   y_coordinate_to_be_modified
	ZL_z=36.01907527540d0       !   z_coordinate_to_be_modified
        RSEP2=(1.80d0)**2 !Distance de separation minimale entre atome (en Ang)

        ANGHOH=1.8239 !Angle en radian de HOH
        FAD1=dcos(ANGHOH)
        LIM_O=NA(1)+NA(2)+1 !Indice du premier O
        LIM_H=NA(1)+NA(2)+NA(3)+1 !Indice du premier H
        NH2O = NA(4)/2 !Nombre de H2O a placer


        CTYP(1)='Si'
        CTYP(2)='B '
        CTYP(3)='OW'
        CTYP(4)='HW'
        CTYP(5)='Zr'
        CTYP(6)='Al'
        CTYP(7)='Ca'
        CMASS(1)=28.085
        CMASS(2)=10.81
        CMASS(3)=15.999
        CMASS(4)=1.008
        CMASS(5)=91.22
        CMASS(6)=26.981
        CMASS(7)=40.08
        CHARGE(1)=4.0
        CHARGE(2)=3.0
        CHARGE(3)=-2.0
        CHARGE(4)=1.0
        CHARGE(5)=4.0
        CHARGE(6)=3.0
        CHARGE(7)=2.0
	IALEA=1

	istpgb=0
	tstep=1.d-15/timau
	timel=0.d0
	imagr=0
        zero=0.d0

	ZLS2_x=ZL_x/2.0d0
	ZLS2_y=ZL_y/2.0d0
	ZLS2_z=ZL_z/2.0d0
      NA11=NA(1)+1
      NA12=NA(1)+NA(2)
      NA21=NA12+1
      NA13=NA12+NA(3)
      NA31=NA13+1
      NA14=NA13+NA(4)
      NA41=NA14+1
      NA15=NA14+NA(5)
      NA51=NA15+1
      NA16=NA15+NA(6)
      NA61=NA16+1
      NA17=NA16+NA(7)
      NA71=NA17+1
      NA18=NA17+NA(8)
      NA81=NA18+1
      NA19=NA18+NA(9)
      NA91=NA19+1
      DO 40 I=1,NA(1)
   40 ITYP(I)=1             ! Silicium
      DO 41 I=NA11,NA12
   41 ITYP(I)=2             ! Bore
      DO 42 I=NA21,NA13
   42 ITYP(I)=3             ! Oxygene
      DO 43 I=NA31,NA14
   43 ITYP(I)=4             ! Sodium
      DO 44 I=NA41,NA15
   44 ITYP(I)=5             ! Zirconium
      DO 45 I=NA51,NA16
   45 ITYP(I)=6             ! Aluminium
      DO 46 I=NA61,NA17
   46 ITYP(I)=7             ! Calcium
	DO 47 I=NA71,NA18
   47 ITYP(I)=8             ! Lithium
	DO 48 I=NA81,NA19
   48 ITYP(I)=9             ! Cesium
	DO 49 I=NA91,IM
   49 ITYP(I)=10            ! Uranium
      DO I=1,IM
      IOCCUP(I)=0
      VP(I,1)=0.d0
      VP(I,2)=0.d0
      VP(I,3)=0.d0
      FP(I,1)=0.d0
      FP(I,2)=0.d0
      FP(I,3)=0.d0
      ENDDO

	IF (IALEA.EQ.0) THEN
	READ(97) XP
	ELSE
       DO K=1,7001 !PREPAR1
       Z1=rand(0)
       ENDDO
  195 CONTINUE
      IF (MOD(LIM_O,1).EQ.0) WRITE(6,*) 'LIM_O=',LIM_O !Test
! Placement des H2O
      icount=0 
      Z1=rand(0) !Choix d'une direction aleatoire
      Z2=rand(0)
      Z3=rand(0)
      ZZ2=DSQRT(Z1**2+Z2**2+Z3**2)
      Z1=Z1/ZZ2
      Z2=Z2/ZZ2
      Z3=Z3/ZZ2

  196 CONTINUE
      icount=icount+1 
      ZZ=rand(0) !Choix d'une coordonnée Z
      AA = FAD1-Z3*ZZ
      BB = 1.d0 - ZZ*ZZ
      A1 = 1.d0+Z2*Z2/Z1/Z1
      B1 = -2.d0*AA*Z2/Z1/Z1
      C1 = AA*AA/Z1/Z1-BB
      DELTA = B1*B1-4.d0*A1*C1

      IF (DELTA.LE.0.0.and.icount.LE.20) THEN
      GOTO 196
      ENDIF
      IF (DELTA.LE.0.0.and.icount.GT.20) THEN
      GOTO 195
      ENDIF

      IF (DELTA.GT.0.0) THEN
      YY1=(-B1-DSQRT(DELTA))/2.d0/A1
      YY2=(-B1+DSQRT(DELTA))/2.d0/A1
      XRAND=rand(0)
      IF (XRAND.GE.0.5d0) YY=YY2
      IF (XRAND.LT.0.5d0) YY=YY1
      ENDIF
      XX=(AA-Z2*YY)/Z1
      ZZ=(FAD1-Z1*XX-Z2*YY)/Z3
!     WRITE(6,*) 'X:',Z1,Z2,Z3
!     WRITE(6,*) 'XX:',XX,YY,ZZ
      ZPP1=XX
      ZPP2=YY
      ZPP3=ZZ

! Introduire la translation de H2O 
      RZ1=rand(0) !Choix d'une position aleatoire
      RZ2=rand(0)
      RZ3=rand(0)
      POSX1=RZ1*ZL_x-ZLS2_x
      POSX2=RZ2*ZL_y-ZLS2_y
      POSX3=RZ3*ZL_z-ZLS2_z
      Z1=Z1+POSX1
      Z2=Z2+POSX2
      Z3=Z3+POSX3
      IF(Z1.GT.ZLS2_x) Z1=Z1-ZL_x
      IF(Z2.GT.ZLS2_y) Z2=Z2-ZL_y
      IF(Z3.GT.ZLS2_z) Z3=Z3-ZL_z
      IF(Z1.LT.-ZLS2_x) Z1=Z1+ZL_x
      IF(Z2.LT.-ZLS2_y) Z2=Z2+ZL_y
      IF(Z3.LT.-ZLS2_z) Z3=Z3+ZL_z
      ZPP1=ZPP1+POSX1
      ZPP2=ZPP2+POSX2
      ZPP3=ZPP3+POSX3
      IF(ZPP1.GT.ZLS2_x) ZPP1=ZPP1-ZL_x
      IF(ZPP2.GT.ZLS2_y) ZPP2=ZPP2-ZL_y
      IF(ZPP3.GT.ZLS2_z) ZPP3=ZPP3-ZL_z
      IF(ZPP1.LT.-ZLS2_x) ZPP1=ZPP1+ZL_x
      IF(ZPP2.LT.-ZLS2_y) ZPP2=ZPP2+ZL_y
      IF(ZPP3.LT.-ZLS2_z) ZPP3=ZPP3+ZL_z

      NPOIN1=1
      A1=-DSIGN(ZL_x,POSX1) 
      A2=-DSIGN(ZL_y,POSX2)
      A3=-DSIGN(ZL_z,POSX3)
      DO J=1,IM 
      IF (IOCCUP(J).EQ.1) THEN
      C1=POSX1-XP(J,1)
      C2=POSX2-XP(J,2)
      C3=POSX3-XP(J,3)
      IF(DABS(C1).GT.ZLS2_x) C1=C1+A1
      IF(DABS(C2).GT.ZLS2_y) C2=C2+A2
      IF(DABS(C3).GT.ZLS2_z) C3=C3+A3
      R2=C1*C1+C2*C2+C3*C3
      IF (R2.LT.RSEP2) THEN
      NPOIN1=0
      ENDIF
      ENDIF
      ENDDO
      NPOIN2=1
      A1=-DSIGN(ZL_x,Z1)
      A2=-DSIGN(ZL_y,Z2)
      A3=-DSIGN(ZL_z,Z3)
      DO J=1,IM
      IF (IOCCUP(J).EQ.1) THEN
      C1=Z1-XP(J,1)
      C2=Z2-XP(J,2)
      C3=Z3-XP(J,3)
      IF(DABS(C1).GT.ZLS2_x) C1=C1+A1
      IF(DABS(C2).GT.ZLS2_y) C2=C2+A2
      IF(DABS(C3).GT.ZLS2_z) C3=C3+A3
      R2=C1*C1+C2*C2+C3*C3
      IF (R2.LT.RSEP2) THEN
      NPOIN2=0
      ENDIF
      ENDIF
      ENDDO
      NPOIN3=1
      A1=-DSIGN(ZL_x,ZPP1)
      A2=-DSIGN(ZL_y,ZPP2)
      A3=-DSIGN(ZL_z,ZPP3)
      DO J=1,IM
      IF (IOCCUP(J).EQ.1) THEN
      C1=ZPP1-XP(J,1)
      C2=ZPP2-XP(J,2)
      C3=ZPP3-XP(J,3)
      IF(DABS(C1).GT.ZLS2_x) C1=C1+A1
      IF(DABS(C2).GT.ZLS2_y) C2=C2+A2
      IF(DABS(C3).GT.ZLS2_z) C3=C3+A3
      R2=C1*C1+C2*C2+C3*C3
      IF (R2.LT.RSEP2) THEN
      NPOIN3=0
      ENDIF
      ENDIF
      ENDDO
      IF (NPOIN1+NPOIN2+NPOIN3.EQ.3) THEN !Retenir la molecule H2O
      XP(LIM_O,1)=POSX1
      XP(LIM_O,2)=POSX2
      XP(LIM_O,3)=POSX3
      IOCCUP(LIM_O)=1
      LIM_O=LIM_O+1
      XP(LIM_H,1)=Z1
      XP(LIM_H,2)=Z2
      XP(LIM_H,3)=Z3
      IOCCUP(LIM_H)=1
      LIM_H=LIM_H+1
      XP(LIM_H,1)=ZPP1
      XP(LIM_H,2)=ZPP2
      XP(LIM_H,3)=ZPP3
      IOCCUP(LIM_H)=1
      LIM_H=LIM_H+1
      ENDIF

      IF (LIM_O-1.NE.NA(1)+NA(2)+NH2O) GOTO 195 

! Placement des autres atomes
      I=0
      ITIR=0
      IATOM=IM-NH2O*3 !Nombre d atomes encore a placer
      IF (IATOM.NE.0) THEN
      WRITE(6,*) 'Placement des autres atomes',IATOM !TEST
  197 CONTINUE
      IF (MOD(I,1).EQ.0) WRITE(6,*) 'I=',I !Test
      Z1=rand(0)
      Z2=rand(0)
      Z3=rand(0)
      IF (Z1.GE.1.0.OR.Z2.GE.1.0.OR.Z3.GE.1.0) GOTO 197
      IF (Z1.LE.0.0.OR.Z2.LE.0.0.OR.Z3.LE.0.0) GOTO 197
      X1=Z1*ZL_x-ZLS2_x
      X2=Z2*ZL_y-ZLS2_y
      X3=Z3*ZL_z-ZLS2_z
      A1=-DSIGN(ZL_x,X1)
      A2=-DSIGN(ZL_y,X2)
      A3=-DSIGN(ZL_z,X3)
      NPOIN=1

      DO J=1,IM
      IF (IOCCUP(J).EQ.1) THEN
      C1=X1-XP(J,1)
      C2=X2-XP(J,2)
      C3=X3-XP(J,3)
      IF(DABS(C1).GT.ZLS2_x) C1=C1+A1
      IF(DABS(C2).GT.ZLS2_y) C2=C2+A2
      IF(DABS(C3).GT.ZLS2_z) C3=C3+A3
      R2=C1*C1+C2*C2+C3*C3
      IF (R2.LT.RSEP2) THEN
      NPOIN=0
      ENDIF
      ENDIF
      ENDDO

  199 CONTINUE
      IF (NPOIN.EQ.1) THEN
      ITIR=ITIR+1
  194 CONTINUE
      I=I+1
      IF (IOCCUP(I).EQ.1) GOTO 194
      XP(I,1)=X1
      XP(I,2)=X2
      XP(I,3)=X3
      IOCCUP(I)=1
      ENDIF
      IF (ITIR.NE.IATOM) GOTO 197
      ENDIF

      ENDIF !Fin du IF sur IALEA

! Le tableau XP() contient l ensemble des positions atomiques
! Ecriture du fichier CONFIG pour DLPOLY
      OPEN (70,file='CONFIG_H2O',form='formatted')
      WRITE(70,*) 'H2O_CONFIG File'
      WRITE(70,*) '        0         3'
      write(70,'(3f20.15)') ZL_x,ZERO,ZERO
      write(70,'(3f20.15)') ZERO,ZL_y,ZERO
      write(70,'(3f20.15)') ZERO,ZERO,ZL_z

      DO I=1,IM
      write(70,'(a2,11X,I5,3F12.6)') CTYP(ITYP(I)),I,CMASS(ITYP(I)),
     XCHARGE(ITYP(I)),ZERO
      write(70,'(3f20.15)') XP(I,1),XP(I,2),XP(I,3) 
      write(70,'(3f20.15)') VP(I,1),VP(I,2),VP(I,3) 
      write(70,'(3f20.15)') FP(I,1),FP(I,2),FP(I,3) 
      ENDDO


	STOP
	END

