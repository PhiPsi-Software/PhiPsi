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
 
      SUBROUTINE HF_Logic_Checking
      ! Some logical checks and data corrections for hydraulic fracturing analysis, the currently
      ! implemented inspection items include:
      ! (1) At the initial moment, it is obvious that friction cracks cannot contain water (only for the
      ! third type of natural cracks);
      ! (2) The current version of the program can only have one initial water-driven fracture.
      
      !*****************************
      ! Read public variable module
      !*****************************
      use Global_Float_Type
      use Global_Common
      use Global_Filename
      use Global_Crack
      use Global_Crack_Common
      use Global_HF
      use Global_Model
      
      implicit none
      integer i_C,c_Na_Cr_num,Count_Fluid_Driven_Cr
      
      print *, "    Logic checking of the HF analysis...." 
      
      !*****************************************************
      ! Check item 1: Friction cracks cannot contain water.
      !*****************************************************
 1001 FORMAT
     & (5X,'-- Frictioanl crack ',I4,' should not be fluid-driven!')  
      if(num_Na_Crack>=1 .and. Key_Na_Crack_Type==3)then
          do i_C =1,num_Crack
              ! Corresponding natural fracture number
              c_Na_Cr_num = Cracks_fric_NF_num(i_C)
              ! If the current calculated fracture number corresponds to a natural friction fracture, check
              ! whether it contains water. If it contains water, it needs to be adjusted to be water-free.
              if (c_Na_Cr_num /= 0) then
                  if(Cracks_HF_State(i_C)==1)then
                      Cracks_HF_State(i_C)=0
                      write(*,1001) i_C
                  endif
              endif
          enddo
      endif

      !***************************************************************************************
      ! Check iTerm 2: The current version of the program can only have one initial hydraulic
      ! fracture.
      !***************************************************************************************
      Count_Fluid_Driven_Cr = 0
      do i_C =1,num_Crack
          if(Cracks_HF_State(i_C)==1)then
              Count_Fluid_Driven_Cr = Count_Fluid_Driven_Cr +1
          endif
      enddo
      if(Count_Fluid_Driven_Cr >1)then
          print *,'    -- Error:: more than one fluid-driven'
     &                    // ' fractures found, which is illegal!'
          call Warning_Message('S',Keywords_Blank)
      endif
      
      !***************************************************************************************
      ! Check iTerm 2: The current version of the program can only have one initial hydraulic
      ! fracture.
      !***************************************************************************************
      Count_Fluid_Driven_Cr = 0
      do i_C =1,num_Crack
          if(Cracks_HF_State(i_C)==1)then
              Count_Fluid_Driven_Cr = Count_Fluid_Driven_Cr +1
          endif
      enddo
      if(Count_Fluid_Driven_Cr >1)then
          print *,'    -- Error:: more than one fluid-driven'
     &                    // ' fractures found, which is illegal!'
          call Warning_Message('S',Keywords_Blank)
      endif
      RETURN
      END SUBROUTINE HF_Logic_Checking
