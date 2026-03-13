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
 
      SUBROUTINE Force_Factor(Lambda,isub,Yes_Last_Growth)
c     Calculate the load factor.

c     ----------------------------
c     Read public variable module
c     ----------------------------
      use Global_Float_Type
      use Global_Common
      use Global_Dynamic
      
      implicit none
      logical Yes_Last_Growth
      integer isub
      real(kind=FT) Lambda
      
      select case(Key_Force_Control)
      case(1)
            if (Key_Analysis_Type==1 .or.
     &        Key_Analysis_Type==3 .or. Key_Analysis_Type==4)then
              Lambda = ONE
          elseif (Key_Analysis_Type.eq.2) then
              if (isub .le. IDy_Num_force_Itr) then
                  Lambda = ONE
            else 
            Lambda = ZR
            end if
          elseif (Key_Analysis_Type.eq.6) then
              if (isub .le. EDy_Num_force_Itr) then
                  Lambda = ONE
            else 
            Lambda = ZR
            end if
          end if
      case(2)
            if (Key_Analysis_Type==1 .or.
     &        Key_Analysis_Type==3 .or. Key_Analysis_Type==4)then
              Lambda = Lambda + ONE/Num_Substeps
          elseif (Key_Analysis_Type.eq.2) then
              if (isub .le. IDy_Num_force_Itr) then
                  Lambda = ONE
            else 
            Lambda = ZR
            end if
          elseif (Key_Analysis_Type.eq.6) then
              if (isub .le. EDy_Num_force_Itr) then
                  Lambda = ONE
            else 
            Lambda = ZR
            end if
          end if
      case(3)    
            if (Key_Analysis_Type==1 .or.
     &        Key_Analysis_Type==3 .or. Key_Analysis_Type==4)then
              if (Yes_Last_Growth .eqv. .False.) then
                  Lambda = Lambda + ONE/Num_Force_Divs
              else
            Lambda = 1/Num_Force_Divs
               end if
          elseif (Key_Analysis_Type.eq.2) then
              if (isub .le. IDy_Num_force_Itr) then
                  Lambda = ONE
            else 
            Lambda = ZR
            end if
          elseif (Key_Analysis_Type.eq.6) then
              if (isub .le. EDy_Num_force_Itr) then
                  Lambda = ONE
            else 
            Lambda = ZR
            end if
          end if
      case(4)
          Lambda = ONE
      end select
      
      RETURN
      END SUBROUTINE Force_Factor
