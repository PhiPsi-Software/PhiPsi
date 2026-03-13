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
 
      subroutine Cal_Equal_Division_Points_3D(A,B,
     &                            Num_Div,Yes_Include_Endpoint,
     &                            Div_Points,Num_Div_Points)
      !Get the equal diversion points of line AB in three-dimension.
      !Diversion point are arranged from A to B.   
      use Global_Float_Type
      use Global_Crack
      
      implicit none
      
      real(kind=FT),intent(in)::A(3),B(3)
      integer,intent(in)::Num_Div
      logical,intent(in)::Yes_Include_Endpoint
      real(kind=FT),intent(out)::Div_Points(Max_Num_Seg_CalP,3)
      integer,intent(out)::Num_Div_Points
      
      integer i
      real(kind=FT) a_x,a_y,a_z,b_x,b_y,b_z
      
      a_x = A(1)
      a_y = A(2)
      a_z = A(3)
      b_x = B(1)
      b_y = B(2)
      b_z = B(3)
      
      Num_Div_Points = 0
      if (Yes_Include_Endpoint.eqv..False.) then
          do i = 1,Num_Div-1
              Div_Points(i,1) = (i*b_x+(Num_Div-i)*a_x)/Num_Div
              Div_Points(i,2) = (i*b_y+(Num_Div-i)*a_y)/Num_Div
              Div_Points(i,3) = (i*b_z+(Num_Div-i)*a_z)/Num_Div
          end do
          Num_Div_Points = Num_Div-1
      else
          Div_Points(1,:)= A
          Num_Div_Points = 1
          do i = 1,Num_Div-1
              Num_Div_Points  = Num_Div_Points + 1
              Div_Points(Num_Div_Points,1)=
     &                            (i*b_x+(Num_Div-i)*a_x)/Num_Div
              Div_Points(Num_Div_Points,2)=
     &                            (i*b_y+(Num_Div-i)*a_y)/Num_Div
              Div_Points(Num_Div_Points,3)=
     &                            (i*b_z+(Num_Div-i)*a_z)/Num_Div     
          end do 
          Num_Div_Points  = Num_Div_Points + 1
          Div_Points(Num_Div_Points,:) = B
      end if
      
      return 
      end SUBROUTINE Cal_Equal_Division_Points_3D              
