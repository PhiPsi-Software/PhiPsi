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
 
      subroutine Tool_Check_Oscillation_by_6_Variables(Input_Var,
     &                                                 Yes_Oscill)
C     This subroutine check if the input six variables vary in an oscillating form.
c     In other wrods, check if the six variables arrange like this:
c     A,b,A,b,A,b  or  A,B,C,A,B,C

c     The rule for determining whether there is volatility is:
c     First case: for: A, b, A, b, A, b
c     (1) 1, 3, and 5 are basically the same (with an error not exceeding Yita1, e.g., 5%), and 2, 4, and 6 are basically the same (with an error not exceeding Yita1, e.g., 5%)
c     (2) On the premise of satisfying rule (1), 1, 3, and 5 should have significant differences from 2, 4, and 6, with a difference greater than Yita2 (e.g., 10%)
c      Specifically, compare the sum of 1, 3, and 5 with the sum of 2, 4, and 6.
c     Second case: for: A, B, C, A, B, C
c     (1) 1 and 4 are basically the same (with an error not exceeding Yita1, e.g., 5%), and 2 and 5 are basically the same (with an error not exceeding Yita1, e.g., 5%).
c     , 3, and 6 are basically the same (with an error not exceeding Yita1 (e.g., 5%))
c     (2) On the premise of meeting rule (1), there must be a significant difference between at least two of A, B, and C.
      use Global_Float_Type      
      implicit none
      real(kind=FT),intent(in):: Input_Var(6)
      logical,intent(out):: Yes_Oscill
      
      real(kind=FT) Error_1_3_5,Error_2_4_6
      real(kind=FT) Yita1,Yita2
      real(kind=FT) sum_1_3_5,sum_2_4_6,Error_sum
      real(kind=FT) Error_1_4,Error_2_5,Error_3_6
      real(kind=FT) Error_1_2,Error_2_3,Error_1_3
      
      Yes_Oscill = .False.
      Yita1 = 5.0D-2;
      Yita2 = 10.0D-2;
      
      call Tool_Error_Between_3_Variables(Input_Var(1),
     &                                        Input_Var(3),
     &                                        Input_Var(5),
     &                                        Error_1_3_5)
      call Tool_Error_Between_3_Variables(Input_Var(2),
     &                                        Input_Var(4),
     &                                        Input_Var(6),
     &                                        Error_2_4_6)
      if(Error_1_3_5 <= Yita1  .and. Error_2_4_6 <= Yita1 )then
          sum_1_3_5 = Input_Var(1) + Input_Var(3) + Input_Var(5)
          sum_2_4_6 = Input_Var(2) + Input_Var(4) + Input_Var(6)
          Error_sum = abs(sum_1_3_5- sum_2_4_6 ) /
     &                max(sum_1_3_5, sum_2_4_6 )
          if(Error_sum >=Yita2 )then
              Yes_Oscill = .True.
              return
          endif
      endif
      
      Error_1_4 = abs(Input_Var(1)- Input_Var(4) ) /
     &            max(Input_Var(1), Input_Var(4) )
      Error_2_5 = abs(Input_Var(2)- Input_Var(5) ) /
     &            max(Input_Var(2), Input_Var(5) )
      Error_3_6 = abs(Input_Var(3)- Input_Var(6) ) /
     &            max(Input_Var(3), Input_Var(6) )
      if(Error_1_4 <= Yita1  .and. 
     &   Error_2_5 <= Yita1  .and.
     &   Error_3_6 <= Yita1 )    then
          Error_1_2 = abs(Input_Var(1)- Input_Var(2) ) /
     &                max(Input_Var(1), Input_Var(2) )
          Error_2_3 = abs(Input_Var(2)- Input_Var(3) ) /
     &                max(Input_Var(2), Input_Var(3) )
          Error_1_3 = abs(Input_Var(1)- Input_Var(3) ) /
     &                max(Input_Var(1), Input_Var(3) )
          if(Error_1_2 >=Yita2  .or.
     &       Error_2_3 >=Yita2  .or.
     &       Error_1_3 >=Yita2 )then
              Yes_Oscill = .True.
              return
          endif
      endif
      
      return 
      end SUBROUTINE Tool_Check_Oscillation_by_6_Variables                         
