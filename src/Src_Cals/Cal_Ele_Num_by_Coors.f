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
 
      subroutine Cal_Ele_Num_by_Coors(x,y,OUT_Elem)
C     Calculate element number based on coordinates

      use Global_Float_Type
      use Global_Model
      
      implicit none
      
      integer i
      real(kind=FT) x,y
      integer,intent(out)::OUT_Elem
      real(kind=FT) c_x_max,c_x_min,c_y_max,c_y_min
      
      integer c_count,Potent_Elem(20)
      logical Yes_In, Yes_On
      real(kind=FT) xpol(5),ypol(5)
      
      c_count = 0
      OUT_Elem = 0
      do i=1,Num_Elem
          c_x_max = x_max_Elements(i)
          c_x_min = x_min_Elements(i)
          c_y_max = y_max_Elements(i)
          c_y_min = y_min_Elements(i)
          if   ((x.le.c_x_max). and.(x.ge.c_x_min). and.(y.le.c_y_max).
     &      and.(y.ge.c_y_min))then
              c_count = c_count +1
              Potent_Elem(c_count) = i
          end if
      end do 
      
      do i =1,c_count
          xpol(1:4) = G_X_NODES(1:4,Potent_Elem(i))
          ypol(1:4) = G_Y_NODES(1:4,Potent_Elem(i))
          xpol(5)   = G_X_NODES(1,Potent_Elem(i))
          ypol(5)   = G_Y_NODES(1,Potent_Elem(i))
          call Tool_Yes_In_Poly(x,y,xpol,ypol,5,Yes_In)
          call Tool_Yes_On_Poly(x,y,xpol,ypol,5,Yes_On)
          if ((Yes_In.eqv..True.).or.(Yes_On.eqv..True.))then
              OUT_Elem = Potent_Elem(i)
              exit
          end if
      end do
      
      return 
      end SUBROUTINE Cal_Ele_Num_by_Coors                          
