C ***************************************
C *       PROGRAMME  analyse_output.f   *
C ***************************************

C*
C* Usage : ./analyse_output
C*

      module routine
        integer, parameter              :: e      = 8
      end module

      PROGRAM ANALYSE
      use routine
      character(len=116) :: aa1, aa2, aa3, aa4, aa5
      real(kind=e) :: cpu, volume, temp_shl, eng_shl
      real(kind=e) :: vir_shl, alpha, beta, gamm, vir_pmf, press
      real(kind=e) :: time, eng_pv, temp_rot, vir_cfg
      real(kind=e) :: vir_src, vir_cou, vir_bnd, vir_ang, vir_con
      real(kind=e) :: vir_tet, eng_tot, temp_tot
      integer*4 :: istep, ilign

      OPEN (70,file='OUTPUT',form='formatted',iostat=iend) 
      rewind 70
      OPEN (72,file='OUTPUT_ANA',form='formatted',iostat=iend) 
      rewind 72
      read(70,'(A116)',end=999) aa1 
      read(70,'(A116)',end=999) aa2 
      read(70,'(A116)',end=999) aa3 
      read(70,'(A116)',end=999) aa4 
      read(70,'(A116)',end=999) aa5 
      ilign=0
         do while (.true.) !Lecture d un fichier OUTPUT
         aa1=aa2
         aa2=aa3
         aa3=aa4
         aa4=aa5
         read(70,'(A116)',end=999) aa5 
         if (aa5(3:9)=='rolling') then
         read(aa3,*) cpu, volume, temp_shl, eng_shl,
     X vir_shl, alpha, beta, gamm, vir_pmf, press
         read(aa2,*) time, eng_pv, temp_rot, vir_cfg,
     X vir_src, vir_cou, vir_bnd, vir_ang, vir_con, vir_tet 
         read(aa1,*) istep, eng_tot, temp_tot
         if (ilign.eq.0) then
         volume0=volume
         press0=press
         endif
!        write(72,*) istep, volume/volume0, press/press0
         write(6,*) ilign, istep
         write(72,*) istep, volume, press, eng_tot, vir_cfg, vir_src,
     X vir_ang, temp_tot
         ilign=ilign+1
         endif
         enddo !Fin du do while

  999 CONTINUE

      STOP
      END
