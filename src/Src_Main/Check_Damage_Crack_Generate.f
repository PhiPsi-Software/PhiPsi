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
 
      subroutine Check_Damage_Crack_Generate(ifra,iter,Total_Num_G_P,
     &                                       Yes_Generated_out)
c     Determine whether a crack has formed based on the damage factor, and determine the crack 
c     direction according to the direction of the perpendicular of the maximum principal tensile
c     stress at the location of maximum damage.
c     Applicable to, Analysis Type 1, with elastic damage material (Material Type 3)
      
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
      integer,intent(in)::ifra,iter,Total_Num_G_P
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
     
     
      integer Max_Damage_Location
      real(kind=FT) Coor_X_Max_Dam,Coor_Y_Max_Dam,Max_Damage
      real(kind=FT) angle_crack,P1(2),P2(2)
      
c     ********************
c     Interactive display
c     ********************
      print *,'    Checking crack generation for damage material......'
 1001 FORMAT(7X,'-- Crack ',I3,' generated form damage material.')  
      
      ! L_New_Crack = 3.5D0 * Ave_Elem_L_Enrich ! Length of the newly formed crack, coefficient at least
      ! 5.5 times the scale of the enriched element
      L_New_Crack        = 3.5D0*Ave_Elem_L
      Yes_Generated_out = .False.
      
      Max_Damage_Location = maxloc(Damage_Gauss(1:Total_Num_G_P),1)
      Max_Damage          = Damage_Gauss(Max_Damage_Location)
      Coor_X_Max_Dam      = Gauss_CoorX(Max_Damage_Location)
      Coor_Y_Max_Dam      = Gauss_CoorY(Max_Damage_Location)
      
      ! Temporarily only generating one damage crack
      if (Crack_Gen_from_Damage == 0) then
          if(Max_Damage >= Material_Dam2Frac_Value)then
          
              Yes_Generated_out = .True.
              
              Crack_Gen_from_Damage  =  Crack_Gen_from_Damage  +1
              
              
              call Cal_Ele_Num_by_Coors(Coor_X_Max_Dam,Coor_Y_Max_Dam,
     &                                  c_Elem_Num)
              !-------------------------------------------------------
              ! Calculate the principal stresses and their directions
              !-------------------------------------------------------
              !For option 1, if Key_Ave_Stress==1, SHI's formulation:
              Search_R = Factor_Ave_R*Ave_Elem_L_Enrich
              a_Weight = a_Ave_Shi
              !For option 2, if Key_Ave_Stress==2, Ordinary formulation:
              l_Ordina = THR*Ave_Elem_L_Enrich
                                                         !                          Theory and Applications
              call Cal_Weighted_Ave_Stress_of_Point(Key_Ave_Stress,
     &                    [Coor_X_Max_Dam,Coor_Y_Max_Dam],
     &                    Search_R,a_Weight,l_Ordina,
     &                    WA_S_xx,WA_S_yy,WA_S_xy)
              !Calculate the angle of the principal stress.
              call Tool_Principal_Stresses_2D(
     &                    WA_S_xx,WA_S_yy,WA_S_xy,
     &                    WA_S_1,WA_S_3,WA_angle) 
              
              angle_crack = WA_angle+pi/TWO
              ! Update Crack Number
              num_Crack = num_Crack +1
              
              ! Confirm crack coordinates
              P1(1) = Coor_X_Max_Dam + L_New_Crack/TWO*cos(angle_crack)
              P1(2) = Coor_Y_Max_Dam + L_New_Crack/TWO*sin(angle_crack)
              P2(1) = Coor_X_Max_Dam - L_New_Crack/TWO*cos(angle_crack)
              P2(2) = Coor_Y_Max_Dam - L_New_Crack/TWO*sin(angle_crack)
              
              Crack_Coor(num_Crack,1,1:2) = P1
              Crack_Coor(num_Crack,2,1:2) = P2
              Each_Cr_Poi_Num(num_Crack)  = 2
              Cracks_Allow_Propa(num_Crack) = 1
              write(*,1001) num_Crack
          endif
      endif
      
      
      RETURN
      END SUBROUTINE Check_Damage_Crack_Generate