!     ================================================= !
!             ____  _       _   ____  _____   _         !
!            |  _ \| |     |_| |  _ \|  ___| |_|        !
!            | |_) | |___   _  | |_) | |___   _         !
!            |  _ /|  _  | | | |  _ /|___  | | |        !
!            | |   | | | | | | | |    ___| | | |        !
!            |_|   |_| |_| |_| |_|   |_____| |_|        !
!     ================================================= !
!     PhiPsi:     a general-purpose computational       !
!                 mechanics program written in Fortran. !
!     Website:    http://phipsi.top                     !
!     Author:     Shi Fang, Huaiyin Institute of        !
!                 Technology, Huaian, JiangSu, China    !
!     Email:      shifang@hyit.edu.cn                   !
!     ------------------------------------------------- !
!     Please cite the following papers:                 !
!     (1)Shi F., Lin C. Modeling fluid-driven           !
!        propagation of 3D complex crossing fractures   !
!        with the extended finite element method.       !
!        Computers and Geotechnics, 2024, 172, 106482.  !
!     (2)Shi F., Wang D., Li H. An XFEM-based approach  !
!        for 3D hydraulic fracturing simulation         !
!        considering crack front segmentation. Journal  !
!        of Petroleum Science and Engineering, 2022,    !
!        214, 110518.                                   !
!     (3)Shi F., Wang D., Yang Q. An XFEM-based         !
!        numerical strategy to model three-dimensional  !
!        fracture propagation regarding crack front     !
!        segmentation. Theoretical and Applied Fracture !
!        Mechanics, 2022, 118, 103250.                  !
!     (4)Shi F., Liu J. A fully coupled hydromechanical !
!        XFEM model for the simulation of 3D non-planar !
!        fluid-driven fracture propagation. Computers   !
!        and Geotechnics, 2021, 132: 103971.            !
!     (5)Shi F., Wang X.L., Liu C., Liu H., Wu H.A. An  !
!        XFEM-based method with reduction technique     !
!        for modeling hydraulic fracture propagation    !
!        in formations containing frictional natural    !
!        fractures. Engineering Fracture Mechanics,     !
!        2017, 173: 64-90.                              !
!     ------------------------------------------------- !
 
      SUBROUTINE Write_Lfrac_Vars_to_Memory_HF(ifra,Counter_Iter,
     &                        Old_total_Time,total_Time,Yes_Growth)
      ! Store the relevant data at the end moment of the previous fracture step into memory
      
c     ****************************
c     Read public variable module
c     ****************************
      use Global_Float_Type
      use Global_Common
      use Global_Filename
      use Global_Crack
      use Global_HF
      
c     *********************
c     Variable Declaration
c     *********************
      implicit none     
      integer,intent(in)::ifra,Counter_Iter
      real(kind=FT),intent(in)::Old_total_Time(Num_Frac),total_Time
      logical,intent(in)::Yes_Growth(Max_Num_Cr,2)
      logical c_Yes_Last_Growth
      
      print *, '    Writting variables of last fracture step '
     &       // 'to memory....' 
     
c     ************************************************************************
c     If it is the end of the first rupture step, then allocate memory space.
c     ************************************************************************
      if(ifra==1)then
         ALLOCATE(L_Cracks_CalP_Num(Max_Num_Cr))
         ALLOCATE(L_Cracks_CalP_Pres(Max_Num_Cr,Max_Num_Cr_CalP))
         ALLOCATE(L_Cracks_CalP_Aper(Max_Num_Cr,Max_Num_Cr_CalP))
         ALLOCATE(L_Cracks_CalP_Coors(Max_Num_Cr,Max_Num_Cr_CalP,2))
         ALLOCATE(L_Cracks_CalP_Orient(Max_Num_Cr,Max_Num_Cr_CalP))
         ALLOCATE(L_Cracks_CalP_Seg(Max_Num_Cr,Max_Num_Cr_CalP))
         ALLOCATE(L_Cracks_CalP_Elem(Max_Num_Cr,Max_Num_Cr_CalP))
         ALLOCATE(L_Cracks_CalP_Pgra(Max_Num_Cr,Max_Num_Cr_CalP))
         ALLOCATE(L_Cracks_CalP_Velo(Max_Num_Cr,Max_Num_Cr_CalP))
         ALLOCATE(L_Cracks_CalP_Quan(Max_Num_Cr,Max_Num_Cr_CalP))  
         ALLOCATE(L_Cracks_CalP_Conc(Max_Num_Cr,Max_Num_Cr_CalP))
         ALLOCATE(Map_L_Cracks_CalP_Conc(Max_Num_Cr,Max_Num_Cr_CalP))
         ! Initialized to 0
         L_Cracks_CalP_Num(1:Max_Num_Cr) = 0
         L_Cracks_CalP_Pres(1:Max_Num_Cr,1:Max_Num_Cr_CalP) = ZR
         L_Cracks_CalP_Aper(1:Max_Num_Cr,1:Max_Num_Cr_CalP) = ZR
         L_Cracks_CalP_Coors(1:Max_Num_Cr,1:Max_Num_Cr_CalP,1:2)=.0D0
         L_Cracks_CalP_Orient(1:Max_Num_Cr,1:Max_Num_Cr_CalP)=ZR
         L_Cracks_CalP_Seg(1:Max_Num_Cr,1:Max_Num_Cr_CalP)   = 0
         L_Cracks_CalP_Elem(1:Max_Num_Cr,1:Max_Num_Cr_CalP)  = 0
         L_Cracks_CalP_Pgra(1:Max_Num_Cr,1:Max_Num_Cr_CalP)  =ZR
         L_Cracks_CalP_Velo(1:Max_Num_Cr,1:Max_Num_Cr_CalP)  =ZR
         L_Cracks_CalP_Quan(1:Max_Num_Cr,1:Max_Num_Cr_CalP)  =ZR
         L_Cracks_CalP_Conc(1:Max_Num_Cr,1:Max_Num_Cr_CalP)  =ZR
         Map_L_Cracks_CalP_Conc(1:Max_Num_Cr,1:Max_Num_Cr_CalP)=ZR
      endif
c     **********
c     Save Data
c     **********
      if(ifra < Num_Frac)then
          L_Cracks_CalP_Num    = Cracks_CalP_Num
          L_Cracks_CalP_Pres   = Cracks_CalP_Pres
          L_Cracks_CalP_Aper   = Cracks_CalP_Aper
          L_Cracks_CalP_Coors  = Cracks_CalP_Coors
          L_Cracks_CalP_Orient = Cracks_CalP_Orient
          L_Cracks_CalP_Seg    = Cracks_CalP_Seg
          L_Cracks_CalP_Elem   = Cracks_CalP_Elem
          L_Cracks_CalP_Pgra   = Cracks_CalP_Pgra
          L_Cracks_CalP_Velo   = Cracks_CalP_Velo
          L_Cracks_CalP_Quan   = Cracks_CalP_Quan
          L_Cracks_CalP_Conc   = Cracks_CalP_Conc
      endif
c     ***************************************************************
c     If it is the last breaking step, clear the relevant variables.
c     ***************************************************************
      if (any(Yes_Growth).eqv..True.)then
          c_Yes_Last_Growth = .True.
      ! If no cracks have propagated, exit the program.
      else
          c_Yes_Last_Growth = .False.
      endif
      
      if(ifra==Num_Frac .or.  
     &  (c_Yes_Last_Growth.eqv..False..and.ifra<Num_Frac)) then
          DEALLOCATE(L_Cracks_CalP_Num)
          DEALLOCATE(L_Cracks_CalP_Pres)  
          DEALLOCATE(L_Cracks_CalP_Aper)
          DEALLOCATE(L_Cracks_CalP_Coors)
          DEALLOCATE(L_Cracks_CalP_Orient)
          DEALLOCATE(L_Cracks_CalP_Seg)  
          DEALLOCATE(L_Cracks_CalP_Elem) 
          DEALLOCATE(L_Cracks_CalP_Pgra) 
          DEALLOCATE(L_Cracks_CalP_Velo)
          DEALLOCATE(L_Cracks_CalP_Quan)  
          DEALLOCATE(L_Cracks_CalP_Conc)
          DEALLOCATE(Map_L_Cracks_CalP_Conc)
      endif
      
      RETURN
      END SUBROUTINE Write_Lfrac_Vars_to_Memory_HF
