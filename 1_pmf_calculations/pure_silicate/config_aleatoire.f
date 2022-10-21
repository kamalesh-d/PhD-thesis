! *****************************************
! *       PROGRAMME  config_aleatoire.f   *
! *****************************************
!  Preparation d une structure de 400 atomes aleatoires

      PROGRAM PREPAR
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(IM=390)
      DIMENSION NA(10)
      DIMENSION XP(IM,3),VP(IM,3)
	DIMENSION ITYP(IM)
	CHARACTER*6 filename
	CHARACTER*3 CTYP(10)
	data bohr / 0.5291772D0/
	data timau / 2.41889D-17/
	NA(1)=130   !Nombre de Si
	NA(2)=0     !Nombre de B
	NA(3)=260   !Nombre de O
	NA(4)=0     !Nombre de Na
	NA(5)=0     !Nombre de Zr
	NA(6)=0     !Nombre de Al
	NA(7)=0     !Nombre de Ca
	NA(8)=0
	NA(9)=0
	NA(10)=0
        ZL=18.06d0    !Arete de la boite de simulation
        RSEP2=(1.60d0)**2 !Distance de separation minimale entre atome (en Ang)


        CTYP(1)='Si'
        CTYP(2)='B '
        CTYP(3)='O '
        CTYP(4)='Na'
        CTYP(5)='Zr'
        CTYP(6)='Al'
        CTYP(7)='Ca'
	IALEA=1

	istpgb=0
	tstep=1.d-15/timau
	timel=0.d0
	imagr=0
        zero=0.d0

	ZLS2=ZL/2.0d0
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
      VP(I,1)=0.d0
      VP(I,2)=0.d0
      VP(I,3)=0.d0
      ENDDO

	IF (IALEA.EQ.0) THEN
	READ(97) XP
	ELSE
       I=0
!      DO K=1,5233 !PREPAR1
!      DO K=1,13220 !PREPAR2
!      DO K=1,18444 !PREPAR3
       DO K=1,6554  !number_to_be_changed_in_config_aleatoire
       Z1=rand(0)
       ENDDO
  197 CONTINUE
      IF (MOD(I,1).EQ.0) WRITE(6,*) 'I=',I !Test
      Z1=rand(0)
      Z2=rand(0)
      Z3=rand(0)
      IF (Z1.GE.1.0.OR.Z2.GE.1.0.OR.Z3.GE.1.0) GOTO 197
      IF (Z1.LE.0.0.OR.Z2.LE.0.0.OR.Z3.LE.0.0) GOTO 197
      X1=Z1*ZL-ZLS2
      X2=Z2*ZL-ZLS2
      X3=Z3*ZL-ZLS2
      A1=-DSIGN(ZL,X1)
      A2=-DSIGN(ZL,X2)
      A3=-DSIGN(ZL,X3)
      NPOIN=1

      DO J=1,I-1
      C1=X1-XP(J,1)
      C2=X2-XP(J,2)
      C3=X3-XP(J,3)
      IF(DABS(C1).GT.ZLS2) C1=C1+A1
      IF(DABS(C2).GT.ZLS2) C2=C2+A2
      IF(DABS(C3).GT.ZLS2) C3=C3+A3
      R2=C1*C1+C2*C2+C3*C3
      IF (R2.LT.RSEP2) THEN
      NPOIN=0
      ENDIF
      ENDDO

  199 CONTINUE
      IF (NPOIN.EQ.1) THEN
	I=I+1
	XP(I,1)=X1
      XP(I,2)=X2
      XP(I,3)=X3
      ENDIF
      IF (I.NE.IM) GOTO 197
	WRITE(97) XP
	ENDIF

!     DO I=1,IM !Les positions sont remises entre O et ZL=arete de la boite
!     XP(I,1)=XP(I,1)+ZLS2
!     XP(I,2)=XP(I,2)+ZLS2
!     XP(I,3)=XP(I,3)+ZLS2
!     ENDDO

! Le tableau XP() contient l ensemble des positions atomiques
! Ecriture du fichier CONFIG pour DLPOLY
      OPEN (70,file='CONFIG_CJ1',form='formatted')
      WRITE(70,*) 'CJ1 CONFIG File'
      WRITE(70,*) '         1         1'
      write(70,'(3f20.15)') ZL,ZERO,ZERO
      write(70,'(3f20.15)') ZERO,ZL,ZERO
      write(70,'(3f20.15)') ZERO,ZERO,ZL

      DO I=1,IM
      write(70,'(a2,11X,I5)') CTYP(ITYP(I)),I
      write(70,'(3f20.15)') XP(I,1),XP(I,2),XP(I,3) 
      write(70,'(3f20.15)') VP(I,1),VP(I,2),VP(I,3) 
      ENDDO


	STOP
	END

