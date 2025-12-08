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
 
      subroutine Check_Hole_Crack_Generate(ifra,iter,Yes_Generated_out)
c     Check for crack initiation at the edge of the hole according to the weighted average maximum principal tensile stress criterion.
      
      ! Description of some variables:

c     **************
c     Public Module
c     **************
      use Global_Float_Type
      use Global_Crack
      use Global_Crack_Common
      use Global_DISP
      use Global_Elem_Area_Vol
      use Global_Model
      use Global_Common
      use Global_Material
      
c     *********************
c     Variable Declaration
c     *********************
      implicit none
      integer,intent(in)::ifra,iter
      logical,intent(out)::Yes_Generated_out
      integer i_C,j_C,i_Tip,Tip_Elem,mat_num,Old_num_Poi
      
      real(kind=FT) L_New_Crack
      integer i_Hole
      real(kind=FT) c_Hole_x,c_Hole_y,c_Hole_r
      integer num_fineness,i_fine
      real(kind=FT) alpha
      real(kind=FT) c_x,c_y
      integer c_Elem_Num
      real(kind=FT) Search_R,a_Weight,l_Ordina
      real(kind=FT) WA_S_xx,WA_S_yy,WA_S_xy
      real(kind=FT) WA_S_1,WA_S_3,WA_angle
      real(kind=FT) Max_WA_S_1,WA_angle_of_Max,i_fine_of_Max
      real(kind=FT) St_tem,St_Critical
      logical Yes_Generated
      real(kind=FT) c_fine_alpha,c_fine_x,c_fine_y
      real(kind=FT) Point_A(2),Point_B(2)
      real(kind=FT) Crack_P1(2),Crack_P2(2),
     &              Line_AB(2,2),new_Line_AB(2,2)
      
c     ********************
c     Interactive display
c     ********************
      print *,'    Checking crack generation of each hole......'
 1001 FORMAT(7X,'-- Crack ',I3,' generated form hole ',I3,'.')  
      
      L_New_Crack        = 3.5D0*Ave_Elem_L_Enrich
      Yes_Generated_out = .False.
c     *****************************
c     Circulation between cavities
c     *****************************
      do i_Hole =1,num_Hole
          ! The premise is that the number of cracks already generated in the current hole is less than the
          ! allowed number.
          if(num_Hole_Crack_Generated(i_Hole)<=
     &                      Num_Crack_Hole_Generated-1) then
              c_Hole_x = Hole_Coor(i_Hole,1) 
              c_Hole_y = Hole_Coor(i_Hole,2) 
              c_Hole_r = Hole_Coor(i_Hole,3) 
              
              num_fineness = 200
              
              Max_WA_S_1 = -Con_Big_20
              Yes_Generated = .False.
              !....................................................
              ! Looping between each subdivision point on the hole
              !....................................................
              do i_fine = 1,num_fineness+1
                  alpha = TWO*pi/num_fineness*(i_fine-1)
                  c_x = c_Hole_x + c_Hole_r*cos(alpha)
                  c_y = c_Hole_y + c_Hole_r*sin(alpha)
                  call Cal_Ele_Num_by_Coors(c_x,c_y,c_Elem_Num)
                  !Setting controlling parameter
                  !-------------------------------------------------------
                  ! Calculate the principal stresses and their directions
                  !-------------------------------------------------------
                  !For option 1, if Key_Ave_Stress==1, SHI's formulation:
                  Search_R = Factor_Ave_R*Ave_Elem_L_Enrich
                  a_Weight = a_Ave_Shi
                  !For option 2, if Key_Ave_Stress==2, Ordinary formulation:
                  l_Ordina = THR*Ave_Elem_L_Enrich
                                                     !                          Theory and Applications
                  call Cal_Weighted_Ave_Stress_of_Point(
     &                Key_Ave_Stress,
     &                [c_x,c_y],Search_R,a_Weight,l_Ordina,
     &                WA_S_xx,WA_S_yy,WA_S_xy)
                  !Calculate the angle of the principal stress.
                  call Tool_Principal_Stresses_2D(
     &                           WA_S_xx,WA_S_yy,WA_S_xy,
     &                           WA_S_1,WA_S_3,WA_angle)
                  !GUI output.
                  St_tem = WA_S_1
                  mat_num = Elem_Mat(c_Elem_Num)
                  St_Critical  = Material_Para(mat_num,5)
                  !-------------------------------------------------------------
                  ! Find the maximum principal stress that meets the conditions
                  !-------------------------------------------------------------
                  if (St_tem >= St_Critical)then
                      Yes_Generated = .True.
                      if (St_tem > Max_WA_S_1)then
                          Max_WA_S_1 = St_tem
                          WA_angle_of_Max = WA_angle
                          i_fine_of_Max = i_fine
                      endif
                  endif             
              enddo
              !......................
              ! Generate a new crack
              !......................
              if (Yes_Generated) then
                  Yes_Generated_out = .True.
                  num_Crack = num_Crack +1
                  ! Confirm crack coordinates
                  c_fine_alpha = TWO*pi/num_fineness*(i_fine_of_Max-1)
                  c_fine_x = c_Hole_x + c_Hole_r*cos(c_fine_alpha)
                  c_fine_y = c_Hole_y + c_Hole_r*sin(c_fine_alpha)
                  Line_AB(1,1:2) = [c_Hole_x,c_Hole_y]
                  Line_AB(2,1:2) = [c_fine_x,c_fine_y]
                  call Tool_Shorten_or_Extend_Line(Line_AB,L_New_Crack,
     &                                          'B',
     &                                          new_Line_AB,Crack_P2)     
                  Crack_Coor(num_Crack,1,1:2) = [c_fine_x,c_fine_y]
                  Crack_Coor(num_Crack,2,1:2) = Crack_P2
                  Each_Cr_Poi_Num(num_Crack)  = 2
                  Cracks_Allow_Propa(num_Crack) = 1
                  ! Information Update
                  num_Hole_Crack_Generated(i_Hole) =
     &                 num_Hole_Crack_Generated(i_Hole)  +1   
                  Hole_Crack_Generated_num(i_Hole,
     &                 num_Hole_Crack_Generated(i_Hole)) = num_Crack     
                  write(*,1001) num_Crack,i_Hole
              endif
          endif
      enddo

      
      RETURN
      END SUBROUTINE Check_Hole_Crack_Generate