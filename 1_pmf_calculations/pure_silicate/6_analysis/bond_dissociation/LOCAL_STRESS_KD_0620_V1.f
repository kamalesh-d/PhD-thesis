!Usage : executable -filepos HISTORY -filefield FILE_FIELD -filecontrol CONTROL
! -idepla CONFIG REVCON       
!Introduction of atom type Al - 03/07/2019 (MODIF_1)
!Print tau for all the types of atoms  - 03/07/2019 (MODIF_2)

! Juillet 2019 : Correction du code pour rajouter Al. Il faut utiliser ce code a l'avenir. Il y avait aussi un probleme de debordement de tableaux sur le calcul des distributions.

! Juin 2020 : Nouveau programme LOCAL_STRESS_KD_0620_V1.f : Programme basé sur LOCAL_STRESS_AJ_V2.f90 (code utilisé par Amreen avec des potentiels dont les paramètres sont donnés dans FIELD)
! Le programme LOCAL_STRESS_KD_0620.f utilise les forces tabulées dans TABLE
! La prise en compte de termes à trois corps de type XOY a également été rajoutée (cas lcas=2)
! Unités des termes de paires et des termes à trois corps : eV/Ang
! Des bugs ont été corrigés par rapport à la première version du code, et le nouveau code s'appelle V1.


      PROGRAM STRESS_kasi

      INTEGER, PARAMETER :: e = 8
      REAL*8, PARAMETER :: avog = 6.02d23
      REAL*8, PARAMETER :: elec_volt = 1.602d-19
      CHARACTER(len=100), DIMENSION(:), ALLOCATABLE :: buffl
      CHARACTER(len=116) :: AA4, AA5
      CHARACTER(len=30) :: nomfilepos,nomfilefield,nomfilecontrol
      CHARACTER(len=30) :: nomfileconfig,nomfilerevcon
      INTEGER :: idepla
      CHARACTER(len=2), DIMENSION(7) ::ctyp    
      CHARACTER(len=2)::C2, C21, C22, C31, C32, C33
      CHARACTER(len=3) :: CHAINE3
      CHARACTER(len=4) :: CHAINE4
      CHARACTER(len=6) :: CHAIN6
      CHARACTER(len=7) :: CHAINE7
      CHARACTER(len=8) :: CHAINE8
      CHARACTER(len=9) :: CHAINE9
      CHARACTER(len=11) :: CHAIN11
      CHARACTER(len=14) :: CPAIRE1
      CHARACTER(len=41) :: fileout
      CHARACTER(len=116) :: AA1, AA2, AA3
      CHARACTER(len=100) :: buff 
      INTEGER, DIMENSION(45000) :: ityp
      INTEGER, DIMENSION(7) :: na
      INTEGER, DIMENSION(9) :: lcas
      INTEGER :: narg, i, ifilepos, j ,k, n
      INTEGER :: iconfig, istep, natom, iform, indice
      INTEGER :: imolecule, nmolecule,ipairs,itrips, iti
      INTEGER :: f, ifilefield, ifilecontrol, ifilehistory
      INTEGER :: ik, ipo, ill
      REAL(kind=e), DIMENSION(45000,9) :: stresl2, stresl3
      REAL(kind=e), DIMENSION(45000,9):: Total_stress
      REAL(kind=e), DIMENSION(45000,9) :: Ke
      REAL(kind=e), DIMENSION(45000,9):: Total_stress_KE
      REAL(kind=e), DIMENSION(45000,9):: Tstrsa ! Average of stress tensor
      REAL(kind=e), DIMENSION(45000,9):: Tstrsa_KE !Average of stree tensor+KE
      REAL(kind=e), DIMENSION(45000,3) :: xp, vp, fp 
      REAL(kind=e), DIMENSION(45000,3) :: xp1, vp1, fp1 
      REAL(kind=e), DIMENSION(45000,3) :: xp2, vp2, fp2
      REAL(kind=e), DIMENSION(15000,28) :: vvdw, gvdw
      REAL(kind=e), DIMENSION(100,7):: countyp    !MODIF_1
      REAL(kind=e), DIMENSION(100,7):: countyp_KE !MODIF_1
      REAL(kind=e), DIMENSION(100,7):: countyp_t    !MODIF_1
      REAL(kind=e), DIMENSION(100,7):: countyp_t_KE !MODIF_1
      REAL(kind=e), DIMENSION(45000) :: dis
      REAL(kind=e), DIMENSION(45000) :: hyd_pres, sig_at, tau_at
      REAL(kind=e), DIMENSION(45000) :: p_local 
      REAL(kind=e), DIMENSION(45000) :: hyd_pres_KE
      REAL(kind=e), DIMENSION(45000) :: sig_at_KE
      REAL(kind=e), DIMENSION(45000) :: tau_at_KE
      REAL(kind=e), DIMENSION(45000) :: p_local_KE
      !REAL(kind=e), DIMENSION(1000) :: p_ltot_KE
      !REAL(kind=e), DIMENSION(1000) :: p_ltot 
      REAL(kind=e), DIMENSION(1000) :: count, lim_inf, lim_sup
      REAL(kind=e), DIMENSION(1000) :: counter, lim_inft, lim_supt
      REAL(kind=e), DIMENSION(1000) :: lim_inf_KE, lim_sup_KE
      REAL(kind=e), DIMENSION(1000) :: lim_inft_KE,lim_supt_KE
      REAL(kind=e), DIMENSION(1000) :: counter_KE,count_KE
      REAL(kind=e), DIMENSION(28):: A, rho, C 
      REAL(kind=e), DIMENSION(9):: Tstrsa_tot, Tstrsa_tot_KE 
      REAL(kind=e), DIMENSION(9):: LAMBDA,CTHEta0
      REAL(kind=e), DIMENSION(9):: GAMA1,GAMA2
      REAL(kind=e), DIMENSION(9):: RC1tb,RC2tb
      REAL(kind=e), DIMENSION(7) :: q, m   !MODIF_1
      REAL(kind=e) :: en, ABUCK,rhoBUCK,CBUCK
      REAL(kind=e) :: sig_tot, sig_tot_KE
      REAL(kind=e) :: Lstwb,CTstwb,G1stwb,G2stwb,Rc1stwb,Rc2stwb
      REAL(kind=e) :: li_KE, lit_KE
      REAL(kind=e) :: p_min_KE, p_max_KE,t_mn_KE,t_mx_KE
      REAL(kind=e) :: p_tot_KE, t_tot_KE
      REAL(kind=e) :: LX, LY, LZ
      REAL(kind=e) :: LX1, LY1, LZ1
      REAL(kind=e) :: LX2, LY2, LZ2
      REAL(kind=e) :: LXS1,LYS1,LZS1
      REAL(kind=e) :: LXS2,LYS2,LZS2
      REAL(kind=e) :: mi, qi,omi
      REAL(kind=e) :: RC1, RC2, li, lit
      REAL(kind=e) :: p_min, p_max,t_mn,t_mx
      REAL(kind=e) :: p_tot, t_tot
      REAL(kind=e) :: x1,x2,x3,v1,v2,v3,f1,f2,f3
      REAL(kind=e) :: y1,y2,y3
      REAL(kind=e) :: x4,x5,x6,v4,v5,v6,f4,f5,f6
      REAL(kind=e) :: x7,x8,x9,v7,v8,v9,f7,f8,f9
      REAL(kind=e) :: u1,u2,u3,u4,ff1,ff2,ff3,ff4
      REAL(kind=e) :: rcut,delta
      REAL(kind=e) :: disx,disy,disz

! Lecture des options

      narg=iargc()
      allocate(buffl(narg))  !Tableau des options
  
      DO i=1,narg                             !#1
      call getarg(i,buff)
      buffl(i) = buff
      END DO                                  !#1
        ifilepos = 0
        ifilefield=0
        ifilecontrol=0
        ifilecontrol=0
        ifileconfig=0
        idepla=0

      DO i=1,narg                             !#2
      IF (buffl(i) == "-filepos") THEN        !#3
      ifilepos=1
      nomfilepos = buffl(i+1)
      END IF                                  !#3

      IF (buffl(i) == "-filefield") THEN      !#4
      ifilefield=1
      nomfilefield = buffl(i+1)
      END IF                                  !#3

      IF (buffl(i) == "-filecontrol") THEN 
      ifilecontrol=1
      nomfilecontrol = buffl(i+1)
      END IF

      IF (buffl(i) == "-idepla") THEN
      idepla=1
      nomfileconfig = buffl(i+1)
      nomfilerevcon = buffl(i+2)
      END IF

      END DO                                  !#2

      deallocate(buffl)                       !Tableau des options

      IF (ifilepos.eq.0) THEN                 !#4
      WRITE(6,*) 'Le fichier des positions n est pas defini'
      STOP
      END IF                                  !#4
      
      IF (ifilefield.eq.0) THEN               !#5
      WRITE(6,*) 'Le fichier field n est pas defini'
      STOP 
      END IF                                  !#5

      IF (ifilecontrol.eq.0) THEN
      WRITE(6,*) 'Le fichier control n est pas defini'
      STOP
      END IF

      WRITE(6,'(24HFichier des positions : ,1X,A30)') nomfilepos
      WRITE(6,'(16HFichier FIELD : ,1X,A30)') nomfilefield
      WRITE(6,'(18HFichier CONTROL : ,1X,A30)') nomfilecontrol

      IF(idepla.eq.1) THEN
      WRITE(6,'(12HFirst file: ,1X,A30)') nomfileconfig
      WRITE(6,'(13HSecond file: ,1X,A30)') nomfilerevcon
      END IF

      fileout=TRIM(nomfilepos)//'_stress_res'
      WRITE(6,'(20HFichier de sortie : ,1X,A31)') fileout
      OPEN (70,file=fileout,form='formatted') !Fichier de sortie

!Initialisation of stress (Averaging stresses on all atoms)
       
        DO i=1,45000
        DO j=1,9
        Tstrsa(i,j)=0.d0
        END DO
        END DO

        DO i=1,45000
        DO j=1,9
        Tstrsa_KE(i,j)=0.d0
        END DO
        END DO


        DO i=1,9
        lcas(i)=0
        lambda(i)=0.d0
        rc1tb(i)=0.d0
        rc2tb(i)=0.d0
        ctheta0(i)=0.d0
        gama1(i)=0.d0
        gama2(i)=0.d0
        END DO

        CTYP(1)='Si'
        CTYP(2)='B '
        CTYP(3)='O '
        CTYP(4)='Na'
        CTYP(5)='Al'
        CTYP(6)='Ca'  
        CTYP(7)='H '  

        iconfig = 0         !Counting configurations for averaging 

! Reading of the empirical potential parameters (File FIELD) 
        
        OPEN (71,file=nomfilefield,form='formatted',iostat=iend)
        REWIND 71
        
        DO I=1,7                            !#8     !MODIF_1
        NA(I)=0
        ENDDO                               !#8

        DO WHILE (.true.)                   !#9

        READ (71,168,end= 998) AA1
        IF (AA1(1:9).eq.'MOLECULES') then   !#10
        READ(AA1,*) CHAINE9,imolecule  ! No of types of molecules 
        
        DO I = 1, imolecule                 !#11 
        READ (71,168,end=998) AA1
  168   format (A116)
        READ (71,168,end=998) AA1
        READ(AA1,*) CHAINE7, nmolecule ! no of atoms of different types
        READ (71,168,end=998) AA1
        READ (71,168,end=998) AA1
        READ(AA1,*) C2,mi,qi           ! Symbol of element, mass, charge

            DO j=1,7                        !#12   !MODIF_1
            IF (C2.eq.CTYP(j)) THEN         !# 13
            m(j)=mi
            q(j)=qi
            na(j)=nmolecule
            END IF                          !# 13
            END DO                          !#12

        READ (71,168,end=998) AA1
        END DO                              !#11
        END IF                              !#10

! begin reading parameters of potential( pair potentials)

        IF (AA1(1:3).eq.'VDW') THEN         !#14 Cas d utilisation de TABLE
        READ(AA1,*) CHAINE3,ipairs 
           DO I = 1, ipairs                    !#15
           READ (71,168,end=998) AA1
           END DO                               !#15
        END IF                               !#14

! begin reading parameters of 3body potential
! O-Si-O:1; Si-O-Si:2; O-B-O:3; O-Al-O:4; Al-O-Al:5; Si-O-Al:6 
! Si-O-H:7; Al-O-H:8; H-O-H:9

        IF (AA1(1:3).eq.'TBP') THEN          !#16
        READ(AA1,*) CHAINE3,itrips 
        DO i = 1, itrips                     !#17
        READ (71,168,end=998) AA1
        READ(AA1,*) C31,C32,C33,CHAINE4
        write(6,*) 'C31:',C31,C32,C33,CHAINE4 !TEST
        if (CHAINE4.eq.'stwb') then
        READ(AA1,*) C31,C32,C33,CHAINE4,Lstwb,CTstwb,G1stwb,
     XG2stwb,Rc1stwb
        Rc2stwb=Rc1stwb
        endif
        if (CHAINE4.eq.'stw1') then
        READ(AA1,*) C31,C32,C33,CHAINE4,Lstwb,CTstwb,G1stwb,
     XG2stwb,Rc1stwb,Rc2stwb
        endif

        IF (C31.eq.'O '.and.C32.eq.'Si'.and.C33.eq.'O ')THEN
        LAMBDA(1)= Lstwb
        CTHETA0(1)=CTstwb
        GAMA1(1)=G1stwb
        GAMA2(1)=G2stwb
        RC1tb(1)=Rc1stwb
        RC2tb(1)=Rc2stwb
        if (na(3).ne.0.and.na(1).ne.0) then
        lcas(1)=1
        write(6,*) 'Terme a trois corps OSiO'
        write(6,*) 'Parametres:',Lstwb,CTstwb,G1stwb,G2stwb,
     XRc1stwb,Rc2stwb
        else
        lcas(1)=0
        endif
        END IF

        IF (C31.eq.'Si'.and.C32.eq.'O '.and.C33.eq.'Si')THEN
        LAMBDA(2)= Lstwb
        CTHETA0(2)=CTstwb
        GAMA1(2)=G1stwb
        GAMA2(2)=G2stwb
        RC1tb(2)=Rc1stwb
        RC2tb(2)=Rc2stwb
        if (na(3).ne.0.and.na(1).ne.0) then
        lcas(2)=1
        write(6,*) 'Terme a trois corps SiOSi'
        write(6,*) 'Parametres:',Lstwb,CTstwb,G1stwb,G2stwb,
     XRc1stwb,Rc2stwb
        else
        lcas(2)=0
        endif
        END IF

        IF (C31.eq.'O '.and.C32.eq.'B '.and.C33.eq.'O ')THEN
        LAMBDA(3)= Lstwb
        CTHETA0(3)=CTstwb
        GAMA1(3)=G1stwb
        GAMA2(3)=G2stwb
        RC1tb(3)=Rc1stwb
        RC2tb(3)=Rc2stwb
        if (na(3).ne.0.and.na(2).ne.0) then
        lcas(3)=1
        write(6,*) 'Terme a trois corps OBO'
        write(6,*) 'Parametres:',Lstwb,CTstwb,G1stwb,G2stwb,
     XRc1stwb,Rc2stwb
        else
        lcas(3)=0
        endif
        END IF

        IF (C31.eq.'O '.and.C32.eq.'Al'.and.C33.eq.'O ')THEN
        LAMBDA(4)= Lstwb
        CTHETA0(4)=CTstwb
        GAMA1(4)=G1stwb
        GAMA2(4)=G2stwb
        RC1tb(4)=Rc1stwb
        RC2tb(4)=Rc2stwb
        if (na(3).ne.0.and.na(5).ne.0) then
        lcas(4)=1
        write(6,*) 'Terme a trois corps OAlO'
        write(6,*) 'Parametres:',Lstwb,CTstwb,G1stwb,G2stwb,
     XRc1stwb,Rc2stwb
        else
        lcas(4)=0
        endif
        END IF

        IF (C31.eq.'Al'.and.C32.eq.'O '.and.C33.eq.'Al')THEN
        LAMBDA(5)= Lstwb
        CTHETA0(5)=CTstwb
        GAMA1(5)=G1stwb
        GAMA2(5)=G2stwb
        RC1tb(5)=Rc1stwb
        RC2tb(5)=Rc2stwb
        if (na(3).ne.0.and.na(5).ne.0) then
        lcas(5)=1
        write(6,*) 'Terme a trois corps AlOAl'
        write(6,*) 'Parametres:',Lstwb,CTstwb,G1stwb,G2stwb,
     XRc1stwb,Rc2stwb
        else
        lcas(5)=0
        endif
        END IF

        IF (C31.eq.'Si'.and.C32.eq.'O '.and.C33.eq.'Al')THEN
        LAMBDA(6)= Lstwb
        CTHETA0(6)=CTstwb
        GAMA1(6)=G1stwb
        GAMA2(6)=G2stwb
        RC1tb(6)=Rc1stwb
        RC2tb(6)=Rc2stwb
        if (na(1).ne.0.and.na(3).ne.0.and.na(5).ne.0) then
        lcas(6)=2 !Type XOY
        write(6,*) 'Terme a trois corps SiOAl'
        write(6,*) 'Parametres:',Lstwb,CTstwb,G1stwb,G2stwb,
     XRc1stwb,Rc2stwb
        else
        lcas(6)=0
        endif
        END IF

        IF (C31.eq.'Si'.and.C32.eq.'O '.and.C33.eq.'H ')THEN
        LAMBDA(7)= Lstwb
        CTHETA0(7)=CTstwb
        GAMA1(7)=G1stwb
        GAMA2(7)=G2stwb
        RC1tb(7)=Rc1stwb
        RC2tb(7)=Rc2stwb
        if (na(1).ne.0.and.na(3).ne.0.and.na(7).ne.0) then
        lcas(7)=2 !Type XOY
        write(6,*) 'Terme a trois corps SiOH'
        write(6,*) 'Parametres:',Lstwb,CTstwb,G1stwb,G2stwb,
     XRc1stwb,Rc2stwb
        else
        lcas(7)=0
        endif
        END IF

        IF (C31.eq.'Al'.and.C32.eq.'O '.and.C33.eq.'H ')THEN
        LAMBDA(8)= Lstwb
        CTHETA0(8)=CTstwb
        GAMA1(8)=G1stwb
        GAMA2(8)=G2stwb
        RC1tb(8)=Rc1stwb
        RC2tb(8)=Rc2stwb
        if (na(5).ne.0.and.na(3).ne.0.and.na(7).ne.0) then
        lcas(8)=2 !Type XOY
        write(6,*) 'Terme a trois corps AlOH'
        write(6,*) 'Parametres:',Lstwb,CTstwb,G1stwb,G2stwb,
     XRc1stwb,Rc2stwb
        else
        lcas(8)=0
        endif
        END IF

        IF (C31.eq.'H '.and.C32.eq.'O '.and.C33.eq.'H ')THEN
        LAMBDA(9)= Lstwb
        CTHETA0(9)=CTstwb
        GAMA1(9)=G1stwb
        GAMA2(9)=G2stwb
        RC1tb(9)=Rc1stwb
        RC2tb(9)=Rc2stwb
        if (na(3).ne.0.and.na(7).ne.0) then
        lcas(9)=1
        write(6,*) 'Terme a trois corps HOH'
        write(6,*) 'Parametres:',Lstwb,CTstwb,G1stwb,G2stwb,
     XRc1stwb,Rc2stwb
        else
        lcas(9)=0
        endif
        END IF

        END DO                              !#17
        END IF                              !#16
        END DO                              !#9
 998   CONTINUE

! Lecture de TABLE
        OPEN (85,file='TABLE',form='formatted',iostat=iend)   
        REWIND 85

        read(85,*,end=996) 

        read(85,*,end=996) delta, rcut, ntot
        do ik=1,ipairs
        read(85,1007,end=996) CPAIRE1
 1007   format (A14)
        ipo=0
        if (CPAIRE1.eq.' Si      Si   ') ipo=1
        if (CPAIRE1.eq.' Si       B   ') ipo=2
        if (CPAIRE1.eq.' B       Si   ') ipo=2
        if (CPAIRE1.eq.' Si       O   ') ipo=3
        if (CPAIRE1.eq.' O       Si   ') ipo=3
        if (CPAIRE1.eq.' Si      Na   ') ipo=4
        if (CPAIRE1.eq.' Na      Si   ') ipo=4
        if (CPAIRE1.eq.' Si      Al   ') ipo=5
        if (CPAIRE1.eq.' Al      Si   ') ipo=5
        if (CPAIRE1.eq.' Si      Ca   ') ipo=6
        if (CPAIRE1.eq.' Ca      Si   ') ipo=6
        if (CPAIRE1.eq.' Si       H   ') ipo=7
        if (CPAIRE1.eq.' H       Si   ') ipo=7
        if (CPAIRE1.eq.' B        B   ') ipo=8
        if (CPAIRE1.eq.' B        O   ') ipo=9
        if (CPAIRE1.eq.' O        B   ') ipo=9
        if (CPAIRE1.eq.' B       Na   ') ipo=10
        if (CPAIRE1.eq.' Na       B   ') ipo=10
        if (CPAIRE1.eq.' B       Al   ') ipo=11
        if (CPAIRE1.eq.' Al       B   ') ipo=11
        if (CPAIRE1.eq.' B       Ca   ') ipo=12
        if (CPAIRE1.eq.' Ca       B   ') ipo=12
        if (CPAIRE1.eq.' B        H   ') ipo=13
        if (CPAIRE1.eq.' H        B   ') ipo=13
        if (CPAIRE1.eq.' O        O   ') ipo=14
        if (CPAIRE1.eq.' O       Na   ') ipo=15
        if (CPAIRE1.eq.' Na       O   ') ipo=15
        if (CPAIRE1.eq.' O       Al   ') ipo=16
        if (CPAIRE1.eq.' Al       O   ') ipo=16
        if (CPAIRE1.eq.' O       Ca   ') ipo=17
        if (CPAIRE1.eq.' Ca       O   ') ipo=17
        if (CPAIRE1.eq.' O        H   ') ipo=18
        if (CPAIRE1.eq.' H        O   ') ipo=18
        if (CPAIRE1.eq.' Na      Na   ') ipo=19
        if (CPAIRE1.eq.' Na      Al   ') ipo=20
        if (CPAIRE1.eq.' Al      Na   ') ipo=20
        if (CPAIRE1.eq.' Na      Ca   ') ipo=21
        if (CPAIRE1.eq.' Ca      Na   ') ipo=21
        if (CPAIRE1.eq.' Na       H   ') ipo=22
        if (CPAIRE1.eq.' H       Na   ') ipo=22
        if (CPAIRE1.eq.' Al      Al   ') ipo=23
        if (CPAIRE1.eq.' Al      Ca   ') ipo=24
        if (CPAIRE1.eq.' Ca      Al   ') ipo=24
        if (CPAIRE1.eq.' Al       H   ') ipo=25
        if (CPAIRE1.eq.' H       Al   ') ipo=25
        if (CPAIRE1.eq.' Ca      Ca   ') ipo=26
        if (CPAIRE1.eq.' Ca       H   ') ipo=27
        if (CPAIRE1.eq.' H       Ca   ') ipo=27
        if (CPAIRE1.eq.' H        H   ') ipo=28
        if (ipo.eq.0) then
        write(6,*) 'Probleme avec IPO'
        write(6,*) 'CPAIRE1:',CPAIRE1
        STOP
        else
        write(6,*) 'IPO: ',ipo
        endif
        ill=0
        do j=1,ntot/4
        read(85,1001,end=996) u1,u2,u3,u4
        ill=ill+1
        vvdw(ill,ipo)=u1
        ill=ill+1
        vvdw(ill,ipo)=u2
        ill=ill+1
        vvdw(ill,ipo)=u3
        ill=ill+1
        vvdw(ill,ipo)=u4
        enddo
        ill=0
        do j=1,ntot/4
        read(85,1001,end=996) ff1,ff2,ff3,ff4
        ill=ill+1
        gvdw(ill,ipo)=ff1
        ill=ill+1
        gvdw(ill,ipo)=ff2
        ill=ill+1
        gvdw(ill,ipo)=ff3
        ill=ill+1
        gvdw(ill,ipo)=ff4
        enddo
        enddo ! Fin de la boucle sur ipairs
 1001   format (4(e15.8,1X))

 996    CONTINUE

! Lecture du fichier CONTROL For Cut off radious

        OPEN (80,file=nomfilecontrol,form='formatted',iostat=iend)   
        REWIND 80

        DO WHILE (.true.)
        READ (80,168,end= 997) AA1
        
        IF (AA1(1:6).eq.'cutoff') THEN       ! first cutoff radius
        READ(AA1,*) CHAIN6,RC1 
        END  IF
    
        IF (cutoff.lt.10.d0) THEN
        cutoff=10.d0
        END IF
        IF (AA1(1:11).eq.'rvdw cutoff') THEN ! second cutoff radius 
        READ(AA1,*) CHAINE4,CHAIN6,RC2
        END IF
    
       END DO
 997   CONTINUE

! Lecture du fichier HISTORY 
      
      OPEN (60,file=nomfilepos,form='formatted',iostat=iend)
      REWIND 60
      DO WHILE(.true.)                     !#18
      READ (60,168,end=999) AA2

        IF (AA2(1:8).eq.'timestep') THEN   !#19
        iconfig = iconfig + 1
        READ(AA2,*) CHAINE8,istep,natom,iform ! inside the line timestep

! Read the dimensions of simulation box
          READ (60,168,end=999) AA2
          READ(AA2,*) LX,y2,y3
          READ (60,168,end=999) AA2
          READ(AA2,*) y1,LY,y3
          READ (60,168,end=999) AA2
          READ(AA2,*) y1,y2,LZ
! Fin de la lecture des dimensions

! Read the coordinates of atom
          DO i=1,natom                      !#20
          READ (60,168,end=999) AA2
          READ(AA2,*) C2,indice,x1,x2
           
            DO j=1,7
            IF (C2.eq.CTYP(j)) THEN 
            ityp(indice)=j
            END IF 
            END DO

          READ (60,168,end=999) AA2 ! Position coordinates
          READ(AA2,*) x1,x2,x3 
          XP(indice,1)=x1
          XP(indice,2)=x2
          XP(indice,3)=x3
            IF (iform.ge.1) THEN             !#21
            READ (60,168,end=999) AA2 ! Velocity coordinates
            READ(AA2,*) v1,v2,v3
            VP(indice,1)=v1
            VP(indice,2)=v2
            VP(indice,3)=v3
              IF (iform.ge.2) THEN           !#22
              READ (60,168,end=999) AA2 !Force coordinates
              READ(AA2,*) f1,f2,f3
              FP(indice,1)=f1
              FP(indice,2)=f2
              FP(indice,3)=f3
              END IF                          !#22
              END IF                          !#21
              END DO                          !#20
! Fin de la lecture des atomes pour une configuration

! Calculate the atomic stresses for all configurations and structure 
      CALL FORCE (vvdw,gvdw,ke,stresl2,stresl3,LX,LY,LZ,XP,VP,
     Xnatom,NA,lambda,gama1,gama2,en,ctheta0,Rc1tb,Rc2tb,RC1,RC2,
     Xityp,A,C,rho,q,m,lcas,delta,rcut)

       omi=dfloat(natom)/LX/LY/LZ            !1/omega value

      ! Calculating total stress tensor over all atoms without and with KE

       DO i=1,natom
       DO j=1,9
       Tstrsa(i,j)=Tstrsa(i,j)+0.5*omi*stresl2(i,j)+
     Xomi*stresl3(i,j)
       END DO
       END DO
    
       DO i=1,natom
       DO j=1,9
       Tstrsa_KE(i,j)=Tstrsa_KE(i,j)+0.5*omi*stresl2(i,j)+
     Xomi*stresl3(i,j)+omi*ke(i,j)*10.00d0/avog/elec_volt
       END DO
       END DO

        END IF                                !#19
        END DO                                !Fin de la boucle #18

 999  CONTINUE
! Averaging total stress over all configurations    

        DO i=1,natom
        DO j=1,9
        Tstrsa(i,j)=Tstrsa(i,j)/iconfig
        END DO
        END DO

        DO i=1,natom
        DO j=1,9
        Tstrsa_KE(i,j)=Tstrsa_KE(i,j)/iconfig
        END DO
        END DO

!Calculating Hydrostatic pressure: P(i)& shear: Z(i) 

        DO i=1,natom

        hyd_pres(i)=(Tstrsa(i,1)+Tstrsa(i,2)+Tstrsa(i,3))/3.d0

        sig_at(i)=Tstrsa(i,1)*Tstrsa(i,2)-(Tstrsa(i,6)**2.d0)
     X+Tstrsa(i,2)*Tstrsa(i,3)-(Tstrsa(i,4)**2.d0)+ 
     XTstrsa(i,1)*Tstrsa(i,3)-(Tstrsa(i,5)**2.d0)

        tau_at(i)=DSQRT(3.d0*(hyd_pres(i)**2.d0)-sig_at(i))

        END DO

         DO i=1,natom

        hyd_pres_KE(i)=(Tstrsa_KE(i,1)+Tstrsa_KE(i,2)+
     XTstrsa_KE(i,3))/3.d0

        sig_at_KE(i)=Tstrsa_KE(i,1)*Tstrsa_KE(i,2)-
     X(Tstrsa_KE(i,6)**2.d0)+Tstrsa_KE(i,2)*Tstrsa_KE(i,3)
     X-(Tstrsa_KE(i,4)**2.d0)+Tstrsa_KE(i,1)*Tstrsa_KE(i,3)
     X-(Tstrsa(i,5)**2.d0)

        tau_at_KE(i)=DSQRT(3.d0*(hyd_pres_KE(i)**2.d0)-sig_at_KE(i))

        END DO

        OPEN(42, file='Pressure.dat', form='formatted')
        DO i=1,natom
        WRITE(42,*) i,hyd_pres(i),tau_at(i)
        END DO
        CLOSE(42)

        OPEN(43, file='Pressure_KE.dat', form='formatted')
        DO i=1,natom
        WRITE(43,*) i,hyd_pres_KE(i),tau_at_KE(i)
        END DO
        CLOSE(43)
       
! Distribution of stress

        p_max=hyd_pres(1)
        p_min=hyd_pres(1)

        DO i=1,natom
        IF(hyd_pres(i).gt.p_max) p_max=hyd_pres(i)
        IF(hyd_pres(i).lt.p_min) p_min=hyd_pres(i)
        END DO

! length of interval

        li=(p_max-p_min)/100.d0

! counting the points for different types of atoms

        DO n=1,100
        DO j=1,7
        countyp(n,j)=0.d0
        END DO
        END DO

        DO n=1,100
        lim_inf(n)=(n-1)*li+p_min
        lim_sup(n)=n*li+p_min
        END DO

        DO i=1,natom
        iti=ityp(i)
        n=INT((hyd_pres(i)-p_min)/li)+1
        if (n.gt.100) n=100
        countyp(n,iti)=countyp(n,iti)+1.d0
        END DO

! Distribution of stress_KE

        p_max_KE=hyd_pres_KE(1)
        p_min_KE=hyd_pres_KE(1)

        DO i=1,natom
        IF(hyd_pres_KE(i).gt.p_max_KE) p_max_KE=hyd_pres_KE(i)
        IF(hyd_pres_KE(i).lt.p_min_KE) p_min_KE=hyd_pres_KE(i)
        END DO

! length of interval_KE

        li_KE=(p_max_KE-p_min_KE)/100.d0

! counting the HYD_PRES points for different types of atoms_KE

        DO n=1,100
        DO j=1,7
        countyp_KE(n,j)=0.d0
        END DO
        END DO

        DO n=1,100
        lim_inf_KE(n)=(n-1)*li_KE+p_min_KE
        lim_sup_KE(n)=n*li_KE+p_min_KE
        END DO

        DO i=1,natom
        iti=ityp(i)
        n=INT((hyd_pres_KE(i)-p_min_KE)/li_KE)+1
        if (n.gt.100) n=100
        countyp_KE(n,iti)=countyp_KE(n,iti)+1.d0
        END DO

! counting the HYD_PRES points for all the atoms

        DO n=1,100
        count(n)=0.d0
        END DO

        DO n=1,100
        lim_inf(n)=(n-1)*li+p_min
        lim_sup(n)=n*li+p_min
        END DO

        DO i=1,natom
        n=INT((hyd_pres(i)-p_min)/li)+1
        if (n.gt.100) n=100
        count(n)=count(n)+1.d0
        END DO


! counting the HYD_PRES points for all the atoms

        DO n=1,100
        count_KE(n)=0.d0
        END DO

        DO n=1,100
        lim_inf_KE(n)=(n-1)*li_KE+p_min_KE
        lim_sup_KE(n)=n*li_KE+p_min_KE
        END DO

        DO i=1,natom
        n=INT((hyd_pres_KE(i)-p_min_KE)/li_KE)+1
        if (n.gt.100) n=100
        count_KE(n)=count_KE(n)+1.d0
        END DO
! Writing out put file for distribution of pressure (without and with KE)

        OPEN(12, file='Distribution_hyd_pres.dat', form='formatted')
        DO n=1,100
        WRITE(12,*) (lim_inf(n)+lim_sup(n))/2.d0,count(n),
     X(lim_inf_KE(n)+lim_sup_KE(n))/2.d0,count_KE(n)
        END DO
        CLOSE(12)

! Writing output file for distribution of pressure for individual atoms

        OPEN(13, file='Distribution_hyd_pres_typ_atom.dat', 
     Xform='formatted')
        DO n=1,100
        WRITE(13,*) (lim_inf(n)+lim_sup(n))/2.d0,countyp(n,1),
     Xcountyp(n,2),countyp(n,3), countyp(n,4), countyp(n,5),
     Xcountyp(n,6),countyp(n,7)    !MODIF_1
        END DO
        CLOSE(13)

! Writing output file distribution of pressure for individual atoms with KE

        OPEN(16,file='Distribution_hyd_pres_typ_atom_KE.dat', 
     Xform='formatted')
        DO n=1,100
        WRITE(16,*) (lim_inf_KE(n)+lim_sup_KE(n))/2.d0,countyp_KE(n,1),
     Xcountyp_KE(n,2),countyp_KE(n,3),countyp_KE(n,4),
     Xcountyp_KE(n,5),countyp_KE(n,6),countyp_KE(n,7)   !MODIF_1
        END DO
        CLOSE(16)
            

! Distribution of shear stress

        t_mx=tau_at(1)
        t_mn=tau_at(1)

        DO i=1,natom
        if(tau_at(i).gt.t_mx) t_mx=tau_at(i)
        if(tau_at(i).lt.t_mn) t_mn=tau_at(i)
        END DO

! length of interval

        lit=(t_mx-t_mn)/100.d0

! counting the points
  
        DO n=1,100
        DO j=1,7    !MODIF_2
        countyp_t(n,j)=0.d0   !MODIF_2
        END DO              !MODIF_2
        END DO

        DO n=1,100
        lim_inft(n)=(n-1)*lit + t_mn
        lim_supt(n)=n*lit + t_mn
        END DO

        DO i=1,natom
        iti=ityp(i)         !MODIF_2
        n=INT((tau_at(i)-t_mn)/lit)+1
        if (n.gt.100) n=100
        countyp_t(n,iti)=countyp_t(n,iti)+1.d0  !MODIF_2
        END DO

! Distribution of shear stress_KE

        t_mx_KE=tau_at_KE(1)
        t_mn_KE=tau_at_KE(1)

        DO i=1,natom
        IF(tau_at_KE(i).gt.t_mx_KE) t_mx_KE=tau_at_KE(i)
        IF(tau_at_KE(i).lt.t_mn_KE) t_mn_KE=tau_at_KE(i)
        END DO

! length of interval

        lit_KE=(t_mx_KE-t_mn_KE)/100.d0

! counting the points

        DO n=1,100
        DO j=1,7            !MODIF_2
        countyp_t_KE(n,j)=0.d0        !MODIF_2
        END DO                      !MODIF_2
        END DO

        DO n=1,100
        lim_inft_KE(n)=(n-1)*lit_KE + t_mn_KE
        lim_supt_KE(n)=n*lit_KE + t_mn_KE
        END DO

        DO i=1,natom
        iti=ityp(i)                          !MODIF_2
        n=INT((tau_at_KE(i)-t_mn_KE)/lit_KE)+1
        if (n.gt.100) n=100
        countyp_t_KE(n,iti)=countyp_t_KE(n,iti)+1.d0     !MODIF_2
        END DO
!!!!!!!!!!!!!!!!!!!!!!!!!MODIF_2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! counting the shear_stress points for all the atoms

        DO n=1,100
        count(n)=0.d0
        END DO

        DO n=1,100
        lim_inft(n)=(n-1)*lit+t_mn
        lim_supt(n)=n*lit+t_mn
        END DO

        DO i=1,natom
        n=INT((tau_at(i)-t_mn)/lit)+1
        if (n.gt.100) n=100
        count(n)=count(n)+1.d0
        END DO

! counting the shear_stress_KE points for all the atoms

        DO n=1,100
        count_KE(n)=0.d0
        END DO

        DO n=1,100
        lim_inft_KE(n)=(n-1)*lit_KE+t_mn_KE
        lim_supt_KE(n)=n*lit_KE+t_mn_KE
        END DO

        DO i=1,natom
        n=INT((tau_at_KE(i)-t_mn_KE)/lit_KE)+1
        if (n.gt.100) n=100
        count_KE(n)=count_KE(n)+1.d0
        END DO
!!!!!!!!!!!!!!!!!!!!!!!!!MODIF_2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!writing the dat for shear distribution
        OPEN(17, file='Distribution_shear.dat', form='formatted')
        DO n=1,100
        WRITE(17,*) (lim_inft(n)+lim_supt(n))/2.d0,counter(n),
     X(lim_inft_KE(n)+lim_supt_KE(n))/2.d0,counter_KE(n)
        END DO
        CLOSE(17)

!!!!!!!!!!!!!!!!!!!!!!!!!MODIF_2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Writing output file for distribution of shear stess for individual atoms

        OPEN(18, file='Distribution_shear_tau_typ_atom.dat', 
     Xform='formatted')
        DO n=1,100
        WRITE(18,*) (lim_inft(n)+lim_supt(n))/2.d0,countyp_t(n,1),
     Xcountyp_t(n,2),countyp_t(n,3), countyp_t(n,4), countyp_t(n,5),
     Xcountyp_t(n,6),countyp_t(n,7)    !MODIF_1
        END DO
        CLOSE(18)

! Writing output file distribution of shear stress for individual atoms with
! KE

        OPEN(19,file='Distribution_shear_tau_typ_atom_KE.dat', 
     Xform='formatted')
        DO n=1,100
        WRITE(19,*) (lim_inft_KE(n)+lim_supt_KE(n))/2.d0,
     Xcountyp_t_KE(n,1),countyp_t_KE(n,2),countyp_t_KE(n,3),
     Xcountyp_t_KE(n,4),countyp_t_KE(n,5),countyp_t_KE(n,6),   !MODIF_1
     Xcountyp_t_KE(n,7)  !MODIF_1
        END DO
        CLOSE(19)
!!!!!!!!!!!!!!!!!!!!!!!!!MODIF_2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Calculating local pressure : p_ltot(k)

        !DO k=1,iconfig
        !p_ltot(k)=0.d0
        !END DO

        !DO k=1,iconfig
        !DO i=1,natom
        !p_ltot(k)=p_ltot(k)+((Tstrsa(i,1)+ &
        !Tstrsa(i,2)+Tstrsa(i,3))/3.d0)
        !END DO
        !END DO

        !DO k=1,iconfig
        !p_ltot(k)=p_ltot(k)/natom
        !p_ltot(k)=p_ltot(k)*1602.177d0        ! ev/A° to katm
        !END DO

! Calculating local pressure : p_ltot_KE(k)

        !DO k=1,iconfig
        !p_ltot_KE(k)=0.d0
        !END DO

        !DO k=1,iconfig
        !DO i=1,natom
        !p_ltot_KE(k)=p_ltot_KE(k)+((Tstrsa_KE(i,1)+ &
        !Tstrsa_KE(i,2)+Tstrsa_KE(i,3))/3.d0)
        !END DO
        !END DO

        !DO k=1,iconfig
        !p_ltot_KE(k)=p_ltot_KE(k)/natom
        !p_ltot_KE(k)=p_ltot_KE(k)*1602.177d0  ! ev/A° to katm
        !END DO


! Calculating total pressure p_tot

        p_tot=0.d0

        DO i=1,natom
        p_tot=p_tot+((Tstrsa(i,1)+ 
     X  Tstrsa(i,2)+Tstrsa(i,3))/3.d0)
        END DO

        p_tot=p_tot/natom

!Calculation total shear T_tot

        T_tot=0.d0
        Do k=1,9
        Tstrsa_tot(k)=0.d0
        End do

        Do n=1,natom
        Do k=1,9
        Tstrsa_tot(k) = Tstrsa_tot(k)+Tstrsa(i,k)
        End do
        End do
        
        Do k=1,9
        Tstrsa_tot(k) = Tstrsa_tot(k)/natom
        End do

        sig_tot=Tstrsa_tot(1)*Tstrsa_tot(2)-(Tstrsa_tot(6)**2.d0)
     X  +Tstrsa_tot(2)*Tstrsa_tot(3)-(Tstrsa_tot(4)**2.d0)+ 
     X  Tstrsa_tot(1)*Tstrsa_tot(3)-(Tstrsa_tot(5)**2.d0)

        T_tot=DSQRT(3.d0*(p_tot**2.d0)-sig_tot)
        p_tot=p_tot*1602.177d0              !ev/A° to katm
        T_tot=T_tot*1602.177d0              !ev/A° to katm

! Calculating total pressure with Ke contribution p_tot_KE

        p_tot_KE=0.d0

        DO i=1,natom
        p_tot_KE=p_tot_KE+((Tstrsa_KE(i,1)+ 
     XTstrsa_KE(i,2)+Tstrsa_KE(i,3))/3.d0)
        END DO

        p_tot_KE=p_tot_KE/natom

!Calculation total shear with KE contribution  T_tot_KE

        T_tot_KE=0.d0
        Do k=1,9
        Tstrsa_tot_KE(k)=0.d0
        End do

        Do n=1,natom
        Do k=1,9
        Tstrsa_tot_KE(k) = Tstrsa_tot_KE(k)+Tstrsa_KE(i,k)
        End do
        End do

        Do k=1,9
        Tstrsa_tot_KE(k) = Tstrsa_tot_KE(k)/natom
        End do

        sig_tot_KE=Tstrsa_tot_KE(1)*Tstrsa_tot_KE(2)-
     X(Tstrsa_tot_KE(6)**2.d0)+Tstrsa_tot_KE(2)*Tstrsa_tot_KE(3)-
     X(Tstrsa_tot_KE(4)**2.d0)+ Tstrsa_tot_KE(1)*Tstrsa_tot_KE(3)-
     X(Tstrsa_tot_KE(5)**2.d0) 
        T_tot_KE=DSQRT(3.d0*(p_tot_KE**2.d0)-sig_tot_KE)
        p_tot_KE=p_tot_KE*1602.177d0           !ev/A° to katm
        T_tot_KE=T_tot_KE*1602.177d0           !ev/A° to katm

! Writing data file for average strees for all atoms in all configs

        !OPEN(19, file='Average_pressure', form='formatted')
        !WRITE(19,*) 'Average pressure :',p_tot,p_tot_KE
        
        !DO k=1,iconfig
        !WRITE(19,*) k, p_ltot(k),p_ltot_KE(k)
        !END DO
        !CLOSE(19)

! Writing data file with total pressure and shear

        OPEN(20, file='Average_total_pressure_shear', form='formatted')
        WRITE(20,*) p_tot, T_tot,p_tot_KE, T_tot_KE
        CLOSE(20)

        if (idepla.eq.1) then !Beginning of the displacement calculation
! Lecture du fichier CONFIG
        
        OPEN (40,file=nomfileconfig,form='formatted',iostat=iend)
        REWIND 40

        DO WHILE(.true.)                    !#19.1
        READ (40,168,end=899) AA4
        READ (40,168,end=899) AA4
        READ(AA4,*)iform,istep,natom 

! Lecture des dimensions
          READ (40,168,end=899) AA4
          READ(AA4,*) LX1,y2,y3
          READ (40,168,end=899) AA4
          READ(AA4,*) y1,LY2,y3
          READ (40,168,end=899) AA4
          READ(AA4,*) y1,y2,LZ2
! Fin de la lecture des dimensions

! Lecture des atomes
          DO i=1,natom                      !#20.1
          READ (40,168,end=899) AA4
          READ(AA4,*) C2,indice

            do j=1,7
            if (C2.eq.CTYP(j)) then
            ityp(indice)=j
            endif
            enddo

          READ (40,168,end=899) AA4
          READ(AA4,*) x4,x5,x6
          XP1(indice,1)=x4
          XP1(indice,2)=x5
          XP1(indice,3)=x6
            IF (iform.ge.1) THEN            !#21.1
            READ (40,168,end=899) AA4
            READ(AA4,*) v4,v5,v6
            VP1(indice,1)=v4
            VP1(indice,2)=v5
            VP1(indice,3)=v6
              IF (iform.ge.2) then          !#22.1
              READ (40,168,end=899) AA4
              READ(AA4,*) f4,f5,f6
              FP1(indice,1)=f4
              FP1(indice,2)=f5
              FP1(indice,3)=f6
              END IF                        !#22.1
              END IF                        !#21.1
              END DO                        !#20.1
              END DO                        !#19.1
 899  CONTINUE
      
! Lecture du fichier REVCON
        OPEN (50,file=nomfilerevcon,form='formatted',iostat=iend)
        REWIND 50
        DO WHILE(.true.)                    !#19.2
        READ (50,168,end=888) AA5
        READ (50,168,end=888) AA5
        READ(AA5,*)iform,istep,natom ! inside the line timestep

! Lecture des dimensions
          READ (50,168,end=888) AA5
          READ(AA5,*) LX2,y2,y3
          READ (50,168,end=888) AA5
          READ(AA5,*) y1,LY2,y3
          READ (50,168,end=888) AA5
          READ(AA5,*) y1,y2,LZ2
! Fin de la lecture des dimensions

! Lecture des atomes
          DO i=1,natom                      !#20.2
          READ (50,168,end=888) AA5
          READ(AA5,*) C2,indice

            DO j=1,7
            IF (C2.eq.CTYP(j)) THEN
            ityp(indice)=j
            END IF
            END DO

          READ (50,168,end=888) AA5
          READ(AA5,*) x7,x8,x9
          XP2(indice,1)=x7
          XP2(indice,2)=x8
          XP2(indice,3)=x9
            IF (iform.ge.1) then            !#21.2
            READ (50,168,end=888) AA5
            READ(AA5,*) v7,v8,v9
            VP2(indice,1)=v7
            VP2(indice,2)=v8
            VP2(indice,3)=v9
              IF (iform.ge.2) then          !#22.2
              READ (50,168,end=888) AA5
              READ(AA5,*) f7,f8,f9
              FP2(indice,1)=f7
              FP2(indice,2)=f8
              FP2(indice,3)=f9
              END IF                        !#22.2
              END IF                        !#21.2
              END DO                        !#20.2
              END DO                        !#19.2

 888  CONTINUE
        CALL DISP (LX1,LY1,LZ1,LX2,LY2,LZ2,xp1,xp2,natom,dis)

        OPEN(41, file='Displacement.dat', form='formatted')
        DO i=1,natom
        WRITE(41,*) i,dis(i),hyd_pres(i)
        END DO

       END IF      
       END

        SUBROUTINE DISP (LX1,LY1,LZ1,LX2,LY2,LZ2,xp1,xp2,natom,dis)

        INTEGER, PARAMETER :: e = 8
        INTEGER, DIMENSION(45000) :: ityp
        REAL(kind=e) :: LX1,LY1,LZ1
        INTEGER :: i,natom
        REAL(kind=e) :: LXS1,LYS1,LZS1
        REAL(kind=e) :: LX2,LY2,LZ2
        REAL(kind=e) :: disx,disy,disz
        REAL(kind=e), DIMENSION(45000)::dis
        !REAL(kind=e), DIMENSION(1500)::hyd_pres
        REAL(kind=e), DIMENSION(45000,3) :: xp1
        REAL(kind=e), DIMENSION(45000,3) :: xp2

        LXS1=LX2/2.d0
        LYS1=LY2/2.d0
        LZS1=LZ2/2.d0

        DO i = 1,natom                        !#26.1

         disx=xp1(i,1)-xp2(i,1)
        IF (disx.ge.LXS1) disx = disx - LX2
        IF (disx.le.-LXS1) disx = disx + LX2

        disy=xp1(i,2)-xp2(i,2)
        IF (disy.ge.LYS1) disy = disy - LY2
        IF (disy.le.-LYS1) disy = disy + LY2

         disz=xp1(i,3)-xp2(i,3)
        IF (disz.ge.LZS1) disz = disz - LZ2
        IF (disz.le.-LZS1) disz = disz + LZ2

          dis(i)=dSQRT(disx**2+disy**2+disz**2)

          END DO                               !#26.1

        RETURN
        END SUBROUTINE DISP

! Calculating the pair, 3 body and kinetic contributions

       SUBROUTINE FORCE (vvdw,gvdw,ke,stresl2,stresl3,LX,LY,LZ,
     XXP,VP,natom,NA,lambda,gama1,gama2,en,ctheta0,rc1tb,rc2tb,
     XRC1,RC2,ityp,A,C,rho,q,m,lcas,delta,rcut)
        INTEGER, PARAMETER :: e = 8
        INTEGER, PARAMETER :: M1=3,N=1
        INTEGER, DIMENSION(45000) :: ityp
        INTEGER, DIMENSION(7) :: na
        INTEGER, DIMENSION(9) :: lcas
        INTEGER :: i, j, k, natom, i1, j1, l
        INTEGER :: ibegin,ifin,jbegin,jfin,kbegin,kfin
        INTEGER :: iti,itj,ipaire,icas
        INTEGER :: iloc,jloc
	REAL(kind=e), PARAMETER :: Pi=3.1415927d0
	REAL(kind=e), DIMENSION(45000,3) :: xp
	REAL(kind=e), DIMENSION(45000,3) :: vp
	REAL(kind=e), DIMENSION(45000,9) :: stresl2, stresl3
        REAL(kind=e), DIMENSION(45000,9) :: Ke
        REAL(kind=e), DIMENSION(15000,28) :: vvdw, gvdw
	REAL(kind=e), DIMENSION(9):: CTHEta0
	REAL(kind=e), DIMENSION(9):: LAMBDA,GAMA1,GAMA2
	REAL(kind=e), DIMENSION(9):: Rc1tb, Rc2tb
        REAL(kind=e), DIMENSION(5) :: q, m
        REAL(kind=e), DIMENSION(28):: A, rho, C
	REAL, DIMENSION (N,M1) :: R3
	REAL, DIMENSION (M1,M1) :: F3
	REAL, DIMENSION (M1) :: R2
	REAL, DIMENSION (M1) :: F2
	REAL(kind=e) :: drijx,drijy,drijz,fijx,fijy,fijz
        REAL(kind=e) :: xij,xik,yij,yik,zij,zik,rij,rik,en,rij2
	REAL(kind=e) :: P_1,P_2,fP_2,P_3
	REAL(kind=e) :: Lx,Ly,Lz,rc,phi2
	REAL(kind=e) :: LXS2,LYS2,LZS2
	REAL(kind=e) :: sigma2,sigma3,vp2,vp3,omega
	REAL(kind=e) :: LAMBDAc,Eu, var_U, cosi, angle, RC1, RC2
	REAL(kind=e) :: Tou_xyL1, Tou_xyL2, Tou_xyL3
	REAL(kind=e) :: Tou_xyi, Tou_xyj, Tou_xyk, Tou_xy
	REAL(kind=e) :: Tou_yzL1, Tou_yzL2, Tou_yzL3
	REAL(kind=e) :: Tou_yzi, Tou_yzj, Tou_yzk, Tou_yz
	REAL(kind=e) :: Tou_zxL1, Tou_zxL2, Tou_zxL3
	REAL(kind=e) :: Tou_zxi, Tou_zxj, Tou_zxk, Tou_zx
	REAL(kind=e) :: TIERS, epsi
        REAL(kind=e) :: Fix,Fjx,Fkx
	REAL(kind=e) :: DU_Drij, DU_Drik,DU_DT
	REAL(kind=e) :: Drij_Dxi,Drik_Dxi,Drij_Dxj,Drik_Dxk
	REAL(kind=e) :: DT_Dxi, DT_Dxj, DT_Dxk
	REAL(kind=e) :: PxxL1, PxxL2, PxxL3, Pxxi, Pxxj, Pxxk, Pxx
	REAL(kind=e) :: Fiy, Fjy, Fky
	REAL(kind=e) :: Drij_Dyi, Drik_Dyi, Drij_Dyj, Drik_Dyk
	REAL(kind=e) :: DT_Dyi, DT_Dyj, DT_Dyk
	REAL(kind=e) :: PyyL1, PyyL2, PyyL3, Pyyi, Pyyj, Pyyk, Pyy
	REAL(kind=e) :: Fiz, Fjz, Fkz 
	REAL(kind=e) :: Drij_Dzi, Drik_Dzi, Drij_Dzj, Drik_Dzk
	REAL(kind=e) :: DT_Dzi, DT_Dzj, DT_Dzk
	REAL(kind=e) :: PzzL1, PzzL2, PzzL3, Pzzi, Pzzj, Pzzk, Pzz
	REAL(kind=e) :: delta,rcut, ppp,gk,gk1,gk2,t1,t2,rgamma

!Initialisation of stresses

        LXS2=LX/2.d0
        LYS2=LY/2.d0
        LZS2=LZ/2.d0

        DO i=1, natom
         DO j=1,9
         stresl2(i,j)=0.d0
         stresl3(i,j)=0.d0
              ke(i,j)=0.d0
         END DO 
        END DO

        TIERS = 1.d0/3.d0                  ! dcos(109°47)
        epsi = 1.602d-9/4.d0/3.14159d0/8.854d-12 !en eV.Ang

        DO i = 1, natom                    !#20
        iti=ityp(i)
        DO j = i + 1, natom                !#21
        itj=ityp(j)
        IF (iti.eq.1.and.itj.eq.1) ipaire=1
        IF (iti.eq.1.and.itj.eq.2) ipaire=2
        IF (iti.eq.2.and.itj.eq.1) ipaire=2
        IF (iti.eq.1.and.itj.eq.3) ipaire=3
        IF (iti.eq.3.and.itj.eq.1) ipaire=3
        IF (iti.eq.1.and.itj.eq.4) ipaire=4
        IF (iti.eq.4.and.itj.eq.1) ipaire=4
        IF (iti.eq.1.and.itj.eq.5) ipaire=5
        IF (iti.eq.5.and.itj.eq.1) ipaire=5
        IF (iti.eq.1.and.itj.eq.6) ipaire=6 
        IF (iti.eq.6.and.itj.eq.1) ipaire=6 
        IF (iti.eq.1.and.itj.eq.7) ipaire=7 
        IF (iti.eq.7.and.itj.eq.1) ipaire=7 
        IF (iti.eq.2.and.itj.eq.2) ipaire=8
        IF (iti.eq.2.and.itj.eq.3) ipaire=9
        IF (iti.eq.3.and.itj.eq.2) ipaire=9
        IF (iti.eq.2.and.itj.eq.4) ipaire=10
        IF (iti.eq.4.and.itj.eq.2) ipaire=10
        IF (iti.eq.2.and.itj.eq.5) ipaire=11
        IF (iti.eq.5.and.itj.eq.2) ipaire=11
        IF (iti.eq.2.and.itj.eq.6) ipaire=12
        IF (iti.eq.6.and.itj.eq.2) ipaire=12
        IF (iti.eq.2.and.itj.eq.7) ipaire=13
        IF (iti.eq.7.and.itj.eq.2) ipaire=13
        IF (iti.eq.3.and.itj.eq.3) ipaire=14
        IF (iti.eq.3.and.itj.eq.4) ipaire=15
        IF (iti.eq.4.and.itj.eq.3) ipaire=15
        IF (iti.eq.3.and.itj.eq.5) ipaire=16
        IF (iti.eq.5.and.itj.eq.3) ipaire=16
        IF (iti.eq.3.and.itj.eq.6) ipaire=17
        IF (iti.eq.6.and.itj.eq.3) ipaire=17
        IF (iti.eq.3.and.itj.eq.7) ipaire=18
        IF (iti.eq.7.and.itj.eq.3) ipaire=18
        IF (iti.eq.4.and.itj.eq.4) ipaire=19
        IF (iti.eq.4.and.itj.eq.5) ipaire=20
        IF (iti.eq.5.and.itj.eq.4) ipaire=20
        IF (iti.eq.4.and.itj.eq.6) ipaire=21
        IF (iti.eq.6.and.itj.eq.4) ipaire=21
        IF (iti.eq.4.and.itj.eq.7) ipaire=22
        IF (iti.eq.7.and.itj.eq.4) ipaire=22
        IF (iti.eq.5.and.itj.eq.5) ipaire=23
        IF (iti.eq.5.and.itj.eq.6) ipaire=24
        IF (iti.eq.6.and.itj.eq.5) ipaire=24
        IF (iti.eq.5.and.itj.eq.7) ipaire=25
        IF (iti.eq.7.and.itj.eq.5) ipaire=25
        IF (iti.eq.6.and.itj.eq.6) ipaire=26
        IF (iti.eq.6.and.itj.eq.7) ipaire=27
        IF (iti.eq.7.and.itj.eq.6) ipaire=27
        IF (iti.eq.7.and.itj.eq.7) ipaire=28

        fijx = 0.d0                      !Initialisation des forces
        fijy = 0.d0
        fijz = 0.d0

! Applying Periodic boundary conditions

        xij = xp(I,1) - xp(j,1)
        IF (xij.ge.LXS2) xij = xij - LX
        IF (xij.le.-LXS2) xij = xij + LX
  
        yij = xp(I,2) - xp(j,2)
        IF (yij.ge.LYS2) yij = yij - LY
        IF (yij.le.-LYS2) yij = yij + LY

        zij = xp(I,3) - xp(j,3)
        IF (zij.ge.LZS2) zij = zij - LZ
        IF (zij.le.-LZS2) zij = zij + LZ

         rij2 =  (xij)**2+(yij)**2+(zij)**2
         rij = DSQRT(rij2)

        drijx = xij / rij
        drijy = yij / rij
        drijz = zij / rij

!       IF (rij < rc1) THEN                !#23

!Part2(P_2) in BMH potential
!     fP_2 = epsi * q(iti) * q(itj) * ((1/rij**3)-(1/rij/(rc1**2))) ! En eV/Ang

!projection of j on sphere and calculate the force of the projected !atom... (Wolf method)
!       fijx = fijx + (fP_2 * rij)* drijx  !No sign ‘-‘
!       fijy = fijy + (fP_2 * rij)* drijy
!       fijz = fijz + (fP_2 * rij)* drijz
!       END IF                             !#23

!Part1(P_1) in BMH potential

!       IF (rij < rc2) THEN                !#24
!       P_1 = A(ipaire)*dexp(-rij/rho(ipaire))/rho(ipaire)
!       P_3 = C(ipaire)/(rij**7.d0)

!       fijx = fijx + (P_1 - 6.d0*P_3)*drijx
!       fijy = fijy + (P_1 - 6.d0*P_3)*drijy
!       fijz = fijz + (P_1 - 6.d0*P_3)*drijz

!       END IF                             !#24

        IF (rij < rcut) THEN               !Utilisation de TABLE
        l   = int(rij/delta)
        ppp = rij/delta - dfloat(l)
        gk  = gvdw(l,ipaire)
        gk1 = gvdw(l+1,ipaire)
        gk2 = gvdw(l+2,ipaire)
        t1 = gk  + (gk1 - gk )*ppp
        t2 = gk1 + (gk2 - gk1)*(ppp - 1.0d0)
        rgamma = (t1 + (t2-t1)*ppp*0.5d0)/rij2
! calculate forces
        fijx = rgamma*xij
        fijy = rgamma*yij
        fijz = rgamma*zij
        END IF                       

!calculate tensor product(T)

        R2(1) = xij 
        R2(2) = yij 
        R2(3) = zij 
        F2(1) = fijx 
        F2(2) = fijy 
        F2(3) = fijz 

        stresl2(i,1) = stresl2(i,1) + F2(1)*R2(1)
        stresl2(j,1) = stresl2(j,1) + F2(1)*R2(1)
        stresl2(i,2) = stresl2(i,2) + F2(2)*R2(2)
        stresl2(j,2) = stresl2(j,2) + F2(2)*R2(2)
        stresl2(i,3) = stresl2(i,3) + F2(3)*R2(3)
        stresl2(j,3) = stresl2(j,3) + F2(3)*R2(3)
        stresl2(i,4) = stresl2(i,4) + F2(2)*R2(3)
        stresl2(j,4) = stresl2(j,4) + F2(2)*R2(3)
        stresl2(i,5) = stresl2(i,5) + F2(1)*R2(3)
        stresl2(j,5) = stresl2(j,5) + F2(1)*R2(3)
        stresl2(i,6) = stresl2(i,6) + F2(1)*R2(2)
        stresl2(j,6) = stresl2(j,6) + F2(1)*R2(2)
        stresl2(i,7) = stresl2(i,7) + F2(3)*R2(2)
        stresl2(j,7) = stresl2(j,7) + F2(3)*R2(2)
        stresl2(i,8) = stresl2(i,8) + F2(3)*R2(1)
        stresl2(j,8) = stresl2(j,8) + F2(3)*R2(1)
        stresl2(i,9) = stresl2(i,9) + F2(2)*R2(1)
        stresl2(j,9) = stresl2(j,9) + F2(2)*R2(1)

        END DO                             !#20
        END DO                             !#21

! End of the treatment of the pair potentials

!Calculation of three body term (phi3); A = teta in degrees, T = teta in radians
! Loop on the 3 possible triplets (OSiO, SiOSi, OBO)

      DO icas=1,9                          !#25
      IF (lcas(icas).ne.0) THEN            !#26

      IF (icas.eq.1) THEN                  !Cas OSiO
      Ibegin = 1
      Ifin = NA(1)
      Jbegin = NA(1) + NA(2) +1
      Jfin = NA(1) + NA(2) + NA(3)
      END IF      

      IF (icas.eq.2) THEN                  !Cas SiOSi
      Ibegin = NA(1) + NA(2) +1
      Ifin = NA(1) + NA(2) + NA(3)
      Jbegin = 1
      Jfin = NA(1)
      END IF      

      IF (icas.eq.3) THEN                  !Cas OBO
      Ibegin = NA(1) + 1
      Ifin = NA(1) + NA(2)
      Jbegin = NA(1) + NA(2) +1
      Jfin = NA(1) + NA(2) + NA(3)
      END IF

      IF (icas.eq.4) THEN                  !Cas OAlO
      Ibegin = NA(1)+NA(2)+NA(3)+NA(4)+1
      Ifin = NA(1)+NA(2)+NA(3)+NA(4)+NA(5)
      Jbegin = NA(1)+NA(2)+1
      Jfin = NA(1)+NA(2)+NA(3)
      END IF

      IF (icas.eq.5) THEN                  !Cas AlOAl
      Ibegin = NA(1)+NA(2)+1
      Ifin = NA(1)+NA(2)+NA(3)
      Jbegin = NA(1)+NA(2)+NA(3)+NA(4)+1
      Jfin = NA(1)+NA(2)+NA(3)+NA(4)+NA(5)
      END IF

      IF (icas.eq.6) THEN                  !Cas SiOAl
      Ibegin = NA(1)+NA(2)+1
      Ifin = NA(1)+NA(2)+NA(3)
      Jbegin = 1
      Jfin = NA(1)
      Kbegin = NA(1)+NA(2)+NA(3)+NA(4)+1
      Kfin = NA(1)+NA(2)+NA(3)+NA(4)+NA(5)
      END IF

      IF (icas.eq.7) THEN                  !Cas SiOH
      Ibegin = NA(1)+NA(2)+1
      Ifin = NA(1)+NA(2)+NA(3)
      Jbegin = 1
      Jfin = NA(1)
      Kbegin = NA(1)+NA(2)+NA(3)+NA(4)+NA(5)+NA(6)+1
      Kfin = NA(1)+NA(2)+NA(3)+NA(4)+NA(5)+NA(6)+NA(7)
      END IF

      IF (icas.eq.8) THEN                  !Cas AlOH
      Ibegin = NA(1)+NA(2)+1
      Ifin = NA(1)+NA(2)+NA(3)
      Jbegin = NA(1)+NA(2)+NA(3)+NA(4)+1
      Jfin = NA(1)+NA(2)+NA(3)+NA(4)+NA(5)
      Kbegin = NA(1)+NA(2)+NA(3)+NA(4)+NA(5)+NA(6)+1
      Kfin = NA(1)+NA(2)+NA(3)+NA(4)+NA(5)+NA(6)+NA(7)
      END IF

      IF (icas.eq.9) THEN                  !Cas HOH
      Ibegin = NA(1)+NA(2)+1
      Ifin = NA(1)+NA(2)+NA(3)
      Jbegin = NA(1)+NA(2)+NA(3)+NA(4)+NA(5)+NA(6)+1
      Jfin = NA(1)+NA(2)+NA(3)+NA(4)+NA(5)+NA(6)+NA(7)
      END IF

! Loop on i, j and k to introduce
        IF (lcas(icas).eq.1) THEN           !#36 Terme a trois corps de type XOX

        DO I = ibegin, ifin                 !#27
        DO J = jbegin, jfin                 !#28

         xij=xp(I,1)-xp(j,1)
        IF (xij.ge.LXS2) xij = xij - LX
        IF (xij.le.-LXS2) xij = xij + LX

         yij=xp(I,2)-xp(j,2)
        IF (yij.ge.LYS2) yij = yij - LY
        IF (yij.le.-LYS2) yij = yij + LY

         zij=xp(I,3)-xp(j,3)
        IF (zij.ge.LZS2) zij = zij - LZ
        IF (zij.le.-LZS2) zij = zij + LZ

         rij=dSQRT(xij**2+yij**2+zij**2)

        IF (rij.lt.(0.95d0*rc1tb(icas))) THEN !#29
        DO K=j+1, jfin                        !#30
       xik=xp(I,1)-xp(k,1)
        IF (xik.ge.LXS2) xik = xik - LX
        IF (xik.le.-LXS2) xik = xik + LX

       yik=xp(I,2)-xp(k,2)
        IF (yik.ge.LYS2) yik = yik - LY
        IF (yik.le.-LYS2) yik = yik + LY

       zik=xp(I,3)-xp(k,3)
        IF (zik.ge.LZS2) zik = zik - LZ
        IF (zik.le.-LZS2) zik = zik + LZ

      rik=dSQRT(xik**2+yik**2+zik**2)

       IF (rik.lt.0.95d0*rc2tb(icas)) THEN      !#31
       cosi=(xij*xik+yij*yik+zij*zik)/(rij*rik) !Cosinus de l angle
       angle = dacos(cosi) !angle is in radian

! Stillinger-Weber potential

      LAMBDAc = LAMBDA(icas) *( (cosi - CTHETA0(icas))**2)
      Eu=dEXP((GAMA1(icas)/(rij-rc1tb(icas)))+(GAMA2(icas)/
     X(rik-rc2tb(icas))))
      var_U = LAMBDAc * Eu

       DU_Drij = var_U * (GAMA1(icas) / (rij - rc1tb(icas))**2)
       DU_Drik = var_U * (GAMA2(icas) / (rik - rc2tb(icas))**2)
       DU_DT = LAMBDA(icas) * Eu * (cosi - CTHETA0(icas)) 

        Drij_Dxi = xij/rij
        Drik_Dxi = xik/rik
        Drij_Dyi = yij/rij
        Drik_Dyi = yik/rik
        Drij_Dzi = zij/rij
        Drik_Dzi = zik/rik

        ! calculating Pxx
        PxxL1=DU_Drij*(xij*Drij_Dxi)+2.d0*DU_DT*cosi*(Drij_Dxi**2)
        PxxL2=DU_Drik*(xik*Drik_Dxi)+2.d0*DU_DT*cosi*(Drik_Dxi**2)
        PxxL3 = -4.d0*DU_DT*Drij_Dxi*Drik_Dxi

        ! calculating Pyy 
        PyyL1=DU_Drij*(yij*Drij_Dyi)+2.d0*DU_DT*cosi*(Drij_Dyi**2)
        PyyL2=DU_Drik*(yik*Drik_Dyi)+2.d0*DU_DT*cosi*(Drik_Dyi**2)
        PyyL3=-4.d0*DU_DT*Drij_Dyi*Drik_Dyi

        ! calculating pzz 
        PzzL1=DU_Drij*(zij*Drij_Dzi)+2.d0*DU_DT*cosi*(Drij_Dzi**2)
        PzzL2=DU_Drik*(zik*Drik_Dzi)+2.d0*DU_DT*cosi*(Drik_Dzi**2)
        PzzL3=-4.d0*DU_DT*Drij_Dzi*Drik_Dzi

        ! Calculation of Non Diagonal terms (Tou-xy) 
        !and its contribution to the i, j, k atoms.
        Tou_xyL1=DU_Drij*Drij_Dxi*yij+2.d0*DU_DT*cosi*Drij_Dxi*Drij_Dyi
        Tou_xyL2=DU_Drik*Drik_Dxi*yik+2.d0*DU_DT*cosi*Drik_Dxi*Drik_Dyi
        Tou_xyL3 =-2.d0*DU_DT*(xij*yik+xik*yij)/(rij*rik) 

        !Calculation of Tou_yz and its contribution to the i, j, k atoms.
        Tou_yzL1=DU_Drij*Drij_Dyi*zij+2.d0*DU_DT*cosi*Drij_Dyi*Drij_Dzi
        Tou_yzL2=DU_Drik*Drik_Dyi*zik+2.d0*DU_DT*cosi*Drik_Dyi*Drik_Dzi
        Tou_yzL3 =-2.d0*DU_DT*(yij*zik+yik*zij)/(rij*rik) 

        !Calculation of Tou_zx and its contribution to the i, j, k atoms.
        Tou_zxL1=DU_Drij*Drij_Dzi*xij+2.d0*DU_DT*cosi*Drij_Dzi*Drij_Dxi
        Tou_zxL2=DU_Drik*Drik_Dzi*xik+2.d0*DU_DT*cosi*Drik_Dzi*Drik_Dxi
        Tou_zxL3 =-2.d0*DU_DT*(zij*xik+zik*xij)/(rij*rik) 

!Calculate the stress from 3body terms (stresl3)
        stresl3(i,1)=stresl3(i,1)+(PxxL1+PxxL2+PxxL3)/2.d0
        stresl3(j,1)=stresl3(j,1)+PxxL1/2.d0+PxxL3/4.d0
        stresl3(k,1)=stresl3(k,1)+PxxL2/2.d0+PxxL3/4.d0
        stresl3(i,2)=stresl3(i,2)+(PyyL1+PyyL2+PyyL3)/2.d0
        stresl3(j,2)=stresl3(j,2)+PyyL1/2.d0+PyyL3/4.d0
        stresl3(k,2)=stresl3(k,2)+PyyL2/2.d0+PyyL3/4.d0
        stresl3(i,3)=stresl3(i,3)+(PzzL1+PzzL2+PzzL3)/2.d0
        stresl3(j,3)=stresl3(j,3)+PzzL1/2.d0+PzzL3/4.d0
        stresl3(k,3)=stresl3(k,3)+PzzL2/2.d0+PzzL3/4.d0
        stresl3(i,4)=stresl3(i,4)+(Tou_yzL1+Tou_yzL2+Tou_yzL3)/2.d0
        stresl3(j,4)=stresl3(j,4)+Tou_yzL1/2.d0+Tou_yzL3/4.d0
        stresl3(k,4)=stresl3(k,4)+Tou_yzL2/2.d0+Tou_yzL3/4.d0
        stresl3(i,5)=stresl3(i,5)+(Tou_zxL1+Tou_zxL2+Tou_zxL3)/2.d0
        stresl3(j,5)=stresl3(j,5)+Tou_zxL1/2.d0+Tou_zxL3/4.d0
        stresl3(k,5)=stresl3(k,5)+Tou_zxL2/2.d0+Tou_zxL3/4.d0
        stresl3(i,6)=stresl3(i,6)+(Tou_xyL1+Tou_xyL2+Tou_xyL3)/2.d0
        stresl3(j,6)=stresl3(j,6)+Tou_xyL1/2.d0+Tou_xyL3/4.d0
        stresl3(k,6)=stresl3(k,6)+Tou_xyL2/2.d0+Tou_xyL3/4.d0
        stresl3(i,7)=stresl3(i,7)+(Tou_yzL1+Tou_yzL2+Tou_yzL3)/2.d0
        stresl3(j,7)=stresl3(j,7)+Tou_yzL1/2.d0+Tou_yzL3/4.d0
        stresl3(k,7)=stresl3(k,7)+Tou_yzL2/2.d0+Tou_yzL3/4.d0
        stresl3(i,8)=stresl3(i,8)+(Tou_zxL1+Tou_zxL2+Tou_zxL3)/2.d0
        stresl3(j,8)=stresl3(j,8)+Tou_zxL1/2.d0+Tou_zxL3/4.d0
        stresl3(k,8)=stresl3(k,8)+Tou_zxL2/2.d0+Tou_zxL3/4.d0
        stresl3(i,9)=stresl3(i,9)+(Tou_xyL1+Tou_xyL2+Tou_xyL3)/2.d0
        stresl3(j,9)=stresl3(j,9)+Tou_xyL1/2.d0+Tou_xyL3/4.d0
        stresl3(k,9)=stresl3(k,9)+Tou_xyL2/2.d0+Tou_xyL3/4.d0

        END IF          !#31
        END DO          !#30
        END IF          !#29
        END DO          !#28
        END DO          !#27
        END IF          !#36


        IF (lcas(icas).eq.2) THEN           !#36 Terme a trois corps de type XOY

        DO I = ibegin, ifin                 !#27
        DO J = jbegin, jfin                 !#28

         xij=xp(I,1)-xp(j,1)
        IF (xij.ge.LXS2) xij = xij - LX
        IF (xij.le.-LXS2) xij = xij + LX

         yij=xp(I,2)-xp(j,2)
        IF (yij.ge.LYS2) yij = yij - LY
        IF (yij.le.-LYS2) yij = yij + LY

         zij=xp(I,3)-xp(j,3)
        IF (zij.ge.LZS2) zij = zij - LZ
        IF (zij.le.-LZS2) zij = zij + LZ

         rij=dSQRT(xij**2+yij**2+zij**2)

        IF (rij.lt.(0.95d0*rc1tb(icas))) THEN !#29
        DO K=kbegin, kfin                     !#30
       xik=xp(I,1)-xp(k,1)
        IF (xik.ge.LXS2) xik = xik - LX
        IF (xik.le.-LXS2) xik = xik + LX

       yik=xp(I,2)-xp(k,2)
        IF (yik.ge.LYS2) yik = yik - LY
        IF (yik.le.-LYS2) yik = yik + LY

       zik=xp(I,3)-xp(k,3)
        IF (zik.ge.LZS2) zik = zik - LZ
        IF (zik.le.-LZS2) zik = zik + LZ

      rik=dSQRT(xik**2+yik**2+zik**2)

       IF (rik.lt.0.95d0*rc2tb(icas)) THEN      !#31
       cosi=(xij*xik+yij*yik+zij*zik)/(rij*rik) !Cosinus de l angle
       angle = dacos(cosi) !angle is in radian

! Stillinger-Weber potential

      LAMBDAc = LAMBDA(icas) *( (cosi - CTHETA0(icas))**2)
      Eu=dEXP((GAMA1(icas)/(rij-rc1tb(icas)))+(GAMA2(icas)/
     X(rik-rc2tb(icas))))
      var_U = LAMBDAc * Eu

       DU_Drij = var_U * (GAMA1(icas) / (rij - rc1tb(icas))**2)
       DU_Drik = var_U * (GAMA2(icas) / (rik - rc2tb(icas))**2)
       DU_DT = LAMBDA(icas) * Eu * (cosi - CTHETA0(icas)) 

        Drij_Dxi = xij/rij
        Drik_Dxi = xik/rik
        Drij_Dyi = yij/rij
        Drik_Dyi = yik/rik
        Drij_Dzi = zij/rij
        Drik_Dzi = zik/rik

        ! calculating Pxx
        PxxL1=DU_Drij*(xij*Drij_Dxi)+2.d0*DU_DT*cosi*(Drij_Dxi**2)
        PxxL2=DU_Drik*(xik*Drik_Dxi)+2.d0*DU_DT*cosi*(Drik_Dxi**2)
        PxxL3 = -4.d0*DU_DT*Drij_Dxi*Drik_Dxi

        ! calculating Pyy 
        PyyL1=DU_Drij*(yij*Drij_Dyi)+2.d0*DU_DT*cosi*(Drij_Dyi**2)
        PyyL2=DU_Drik*(yik*Drik_Dyi)+2.d0*DU_DT*cosi*(Drik_Dyi**2)
        PyyL3=-4.d0*DU_DT*Drij_Dyi*Drik_Dyi

        ! calculating pzz 
        PzzL1=DU_Drij*(zij*Drij_Dzi)+2.d0*DU_DT*cosi*(Drij_Dzi**2)
        PzzL2=DU_Drik*(zik*Drik_Dzi)+2.d0*DU_DT*cosi*(Drik_Dzi**2)
        PzzL3=-4.d0*DU_DT*Drij_Dzi*Drik_Dzi

        ! Calculation of Non Diagonal terms (Tou-xy) 
        !and its contribution to the i, j, k atoms.
        Tou_xyL1=DU_Drij*Drij_Dxi*yij+2.d0*DU_DT*cosi*Drij_Dxi*Drij_Dyi
        Tou_xyL2=DU_Drik*Drik_Dxi*yik+2.d0*DU_DT*cosi*Drik_Dxi*Drik_Dyi
        Tou_xyL3 =-2.d0*DU_DT*(xij*yik+xik*yij)/(rij*rik) 

        !Calculation of Tou_yz and its contribution to the i, j, k atoms.
        Tou_yzL1=DU_Drij*Drij_Dyi*zij+2.d0*DU_DT*cosi*Drij_Dyi*Drij_Dzi
        Tou_yzL2=DU_Drik*Drik_Dyi*zik+2.d0*DU_DT*cosi*Drik_Dyi*Drik_Dzi
        Tou_yzL3 =-2.d0*DU_DT*(yij*zik+yik*zij)/(rij*rik) 

        !Calculation of Tou_zx and its contribution to the i, j, k atoms.
        Tou_zxL1=DU_Drij*Drij_Dzi*xij+2.d0*DU_DT*cosi*Drij_Dzi*Drij_Dxi
        Tou_zxL2=DU_Drik*Drik_Dzi*xik+2.d0*DU_DT*cosi*Drik_Dzi*Drik_Dxi
        Tou_zxL3 =-2.d0*DU_DT*(zij*xik+zik*xij)/(rij*rik) 

!Calculate the stress from 3body terms (stresl3)
        stresl3(i,1)=stresl3(i,1)+(PxxL1+PxxL2+PxxL3)/2.d0
        stresl3(j,1)=stresl3(j,1)+PxxL1/2.d0+PxxL3/4.d0
        stresl3(k,1)=stresl3(k,1)+PxxL2/2.d0+PxxL3/4.d0
        stresl3(i,2)=stresl3(i,2)+(PyyL1+PyyL2+PyyL3)/2.d0
        stresl3(j,2)=stresl3(j,2)+PyyL1/2.d0+PyyL3/4.d0
        stresl3(k,2)=stresl3(k,2)+PyyL2/2.d0+PyyL3/4.d0
        stresl3(i,3)=stresl3(i,3)+(PzzL1+PzzL2+PzzL3)/2.d0
        stresl3(j,3)=stresl3(j,3)+PzzL1/2.d0+PzzL3/4.d0
        stresl3(k,3)=stresl3(k,3)+PzzL2/2.d0+PzzL3/4.d0
        stresl3(i,4)=stresl3(i,4)+(Tou_yzL1+Tou_yzL2+Tou_yzL3)/2.d0
        stresl3(j,4)=stresl3(j,4)+Tou_yzL1/2.d0+Tou_yzL3/4.d0
        stresl3(k,4)=stresl3(k,4)+Tou_yzL2/2.d0+Tou_yzL3/4.d0
        stresl3(i,5)=stresl3(i,5)+(Tou_zxL1+Tou_zxL2+Tou_zxL3)/2.d0
        stresl3(j,5)=stresl3(j,5)+Tou_zxL1/2.d0+Tou_zxL3/4.d0
        stresl3(k,5)=stresl3(k,5)+Tou_zxL2/2.d0+Tou_zxL3/4.d0
        stresl3(i,6)=stresl3(i,6)+(Tou_xyL1+Tou_xyL2+Tou_xyL3)/2.d0
        stresl3(j,6)=stresl3(j,6)+Tou_xyL1/2.d0+Tou_xyL3/4.d0
        stresl3(k,6)=stresl3(k,6)+Tou_xyL2/2.d0+Tou_xyL3/4.d0
        stresl3(i,7)=stresl3(i,7)+(Tou_yzL1+Tou_yzL2+Tou_yzL3)/2.d0
        stresl3(j,7)=stresl3(j,7)+Tou_yzL1/2.d0+Tou_yzL3/4.d0
        stresl3(k,7)=stresl3(k,7)+Tou_yzL2/2.d0+Tou_yzL3/4.d0
        stresl3(i,8)=stresl3(i,8)+(Tou_zxL1+Tou_zxL2+Tou_zxL3)/2.d0
        stresl3(j,8)=stresl3(j,8)+Tou_zxL1/2.d0+Tou_zxL3/4.d0
        stresl3(k,8)=stresl3(k,8)+Tou_zxL2/2.d0+Tou_zxL3/4.d0
        stresl3(i,9)=stresl3(i,9)+(Tou_xyL1+Tou_xyL2+Tou_xyL3)/2.d0
        stresl3(j,9)=stresl3(j,9)+Tou_xyL1/2.d0+Tou_xyL3/4.d0
        stresl3(k,9)=stresl3(k,9)+Tou_xyL2/2.d0+Tou_xyL3/4.d0

        END IF          !#31
        END DO          !#30
        END IF          !#29
        END DO          !#28
        END DO          !#27
        END IF          !#36

        END IF          !#26
        END DO          !#25

! Kinetic Energy term calculation

         DO i=1,natom
         DO j=1,9
         Ke(i,j)=0.d0
         END DO
         END DO

        DO i=1, natom
        iti=ITYP(i)
        mi=m(iti)
        Ke(i,1)= mi*vp(i,1)*vp(i,1)
        Ke(i,2)= mi*vp(i,2)*vp(i,2)
        Ke(i,3)= mi*vp(i,3)*vp(i,3)
        Ke(i,4)= mi*vp(i,2)*vp(i,3)
        Ke(i,5)= mi*vp(i,1)*vp(i,3)
        Ke(i,6)= mi*vp(i,1)*vp(i,2)
        Ke(i,7)= mi*vp(i,3)*vp(i,2)
        Ke(i,8)= mi*vp(i,3)*vp(i,1)
        Ke(i,9)= mi*vp(i,2)*vp(i,1)
        END DO
        
        RETURN
        END SUBROUTINE FORCE
