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
 
subroutine Check_Crack_Grows(ifra,iter,Yes_Grow)
! Check whether the crack has propagated according to the maximum circumferential tensile stress
! criterion or the weighted average maximum principal tensile stress criterion, and determine the
! crack propagation direction and the new crack tip.
      
!***************
! Public Module
!***************
use Global_Float_Type
use Global_Crack
use Global_Crack_Common
use Global_DISP
use Global_Elem_Area_Vol
use Global_Model
use Global_HF
use Global_Common
use Global_Material
use Global_Dynamic 
use Global_Inclusion
      
!**********************
! Variable Declaration
!**********************
implicit none
integer,intent(in)::ifra,iter
logical,intent(out)::Yes_Grow(Max_Num_Cr,2)
integer i_C,j_C,i_Tip,Tip_Elem,mat_num,Old_num_Poi
real(kind=FT) cr_Tip(2,2),omega_Tip(2),delta_angle(2),  &
              theta(2),tem1,K_tem,KI_Critical,new_tip(2,2), &
              delta_L_Growth,tem_Crack_Coor(Max_Num_Cr_P,2)
integer EndByEnd_Outline(size(Outline,1)+1) 
logical Yes_In_Model
real(kind=FT) New_Seg_A(2),New_Seg_B(2),c_Seg_C(2),c_Seg_D(2), &
              Intersec_x,Intersec_y,Signed_Dis,c_Seg_CD(2,2), &
              Check_Distance,tried_New_Tip(2,2), &
              Tri_New_Seg_A(2),Tri_New_Seg_B(2), &
              Tri_Intersec_x,Tri_Intersec_y
integer nPt,jPt,num_pt
logical Yes_Cross,Tri_Yes_Cross
! ---------Related to Natural Fracture Intersections--------
integer j_NC,nPt_Na                                                      
logical Yes_HF_met_NF
real(kind=FT) New_Crack(10,2,2)
integer Na_Cr_num(10)
integer Positive_Cr_num(10)
integer num_New_Crack
real(kind=FT) tem_Intersec(2)
integer tem_NF_num,tem_Positive_num
real(kind=FT) NF_Tip1(2),NF_Tip2(2),NF_Tip(2)
real(kind=FT) L_Tip1,L_Tip2
real(kind=FT) L_New_Crack
integer Statu_Tip1,Statu_Tip2,Statu_Tip
integer New_Crack_Statu(30,2)  
integer i_New_C,c_Na_Cr
real(kind=FT) L_C_Tip1,L_C_Tip2,L_C_Tip
! ---------Related to the legality check of new crack tips (crack fragments)--------
integer c_Ele,i_Node
logical Yes_In,Yes_on_Ele_Edge
real(kind=FT) Modifi_Factor,c_Node_x,c_Node_y
logical Yes_Node_OnSeg,Yes_Inj_On_Cr
!---------Parameters about propagation criterion--------
real(kind=FT) Search_R,a_Weight
real(kind=FT) WA_S_xx,WA_S_yy,WA_S_xy,WA_S_1,WA_S_3,WA_angle
real(kind=FT) St_Critical,St_tem,l_Ordina,pi_180,pi_2
real(kind=FT) Possible_Angle,Possible_Delta_Angle
real(kind=FT) New_Na_Tip(2),Vector_J_Tip1(2),Vector_J_Tip2(2)
real(kind=FT) Point_A(2),Vector_J_A(2),dot_J_Tip1,dot_J_Tip2
real(kind=FT) k_I,k_II,tem_K
logical Yes_Tip_Growth
integer c_Na_Tip
logical Yes_In_Frac_zone
!-------------
logical Yes_In_OneEl
integer Old_OUT_Elem,New_OUT_Elem
real(kind=FT) Old_Tip(2),Old_Old_Tip(2),Ori_delta_L_Growth,c_L,c_L2
!-------Holes related--------
real(kind=FT) c_Hole_x,c_Hole_y,c_Hole_r,c_Dis,c_N_x,c_N_y
real(kind=FT) Tool_Function_2Point_Dis

integer c_Na_Crack,Real_Cr_NF, c_Na_Cr_Num 
real(kind=FT) New_Crack_Tip_1(2),New_Crack_Tip_2(2)
real(kind=FT) delta_angle_second,tem3,K_tem1,tem_delta_L
real(kind=FT) Lambda,mui,Velocity_long,Velocity_tran
real(kind=FT) c_E,c_v,c_density,Velocity_Rayl,c_KIc
integer i_Hole
real(kind=FT) x0,y0,R0
integer c_num_Inter,c_State
real(kind=FT) c_Inter(2,2)
integer iii_try
! -------Variables Related to the Formation of Arc Cracks from the Intersection of Cracks and
! Inclusions
logical Yes_In_Incl
integer i_Incl
real(kind=FT) c_Incl_x,c_Incl_y,c_Incl_r
real(kind=FT) c_Inter_P(2,2)
logical Yes_New_Arc_Cr,Yes_Update_Arc_Cr
real(kind=FT) c_Arc_Cr(1:11)
real(kind=FT) c_P1_P2(2,2)
integer c_Point_Status
real(kind=FT) Want_Point(2)
logical c_Yes_feasible
real(kind=FT) c_Arc_r,c_Arc_Radian_S,c_Arc_Radian_E,c_Arc_Radian
real(kind=FT) L_AP1,L_AP2
real(kind=FT) c_Arc_Direction
logical c_Yes_in_Incl
integer c_Incl_Num
real(kind=FT) delta_L_Growth_Arc
logical Extend_Status
real(kind=FT) c_Seed_Number,Factor_L
! Crack and polygon inclusion interaction
integer n_Vertex
logical Yes_in_Poly_Inc
integer i_S 
real(kind=FT) Side_A(2),Side_B(2),c_inter_X,c_inter_Y
logical c_Yes_Cross
integer c_Incl,c_Side
real(kind=FT) cc_Line_AB(2,2),c_new_Line_AB(2,2)
real(kind=FT) Possoble_P_1(2),Possoble_P_2(2)
real(kind=FT) c_O(2),L_P_1_O,L_P_2_O
logical c_Yes_ON
real(kind=FT) c_I(2),c_Alpha
real(kind=FT) WA_S_normal,WA_S_shear
real(kind=FT)  c_dis_1,c_dis_2
real(kind=FT) x0_Hole,y0_Hole,R0_Hole,theta_Hole,a_Hole,b_Hole
integer j_S
real(kind=FT) S1(2),S2(2)
real(kind=FT) S3(2),S4(2)
real(kind=FT) i_C_Tip(2)
real(kind=FT) point_seg_1(2),point_seg_2(2)                 
real(kind=FT) Shear_Strength,Shear_Stress
real(kind=FT) Twist_Angle,Direction_Parameter
real(kind=FT) random_num_for_angle
real(kind=FT) tem_value_k_I
                       
!*********************
! Interactive display
!*********************
print *,'    Checking propagation of each crack...'   
1001 FORMAT(5X,'-- Tip ',I1,' of crack',I4,' propagated!')  
1601 FORMAT(5X,'-- Tip ',I1,' of crack',I4,' propagated and the')  
1602 FORMAT(5X,'   originsl tip turned into an arc tip!')  
1701 FORMAT(5X,'-- Arc tip ',I1,' of crack',I4,' propagated!') 
1002 FORMAT(5X,'-- R_Del_L of ',I1,' of crack',I4, ' is ',E15.7, ' m.')  
1012 FORMAT(5X,'-- Delta_L of ',I1,' of crack',I4, ' is ',E15.7, ' m.')  
1003 FORMAT(5X,'-- Warning :: Tip ',I1,' of crack',I4,' was stopped!')
1004 FORMAT(5X,'               Because it was going beyong the model!')
1005 FORMAT(5X,'               Because it was going beyong fracture zone!')
1013 FORMAT(5X,'-- Caution :: Tip ',I1,' of crack',I4,' met crack',I4,'!')
1023 FORMAT(5X,'-- Caution :: Tip ',I1,' of crack',I4,' met hole',I4,'!')
1024 FORMAT(5X,'              Tip ',I1,' of crack',I4,' stopped growing!')
1014 FORMAT(5X,'              Crack',I4,' stopped growing!')     
2013 FORMAT(5X,'-- Caution :: Tip ',I1,' of crack',I4,' met natural crack',I4,'!')
2014 FORMAT(5X,'              Crack',I4,' stopped growing!')  
2015 FORMAT(5X,'-- Caution :: Fracture tip on element edge modified!')
2016 FORMAT(5X,'-- Caution :: Fracture segment passes through node!')
3001 FORMAT(5X,'-- Max principal tensile stress of tip ',I1,' of crack',I4,' is',F15.7,' MPa.')
3002 FORMAT(5X,'-- Max shear stress of tip ',I1,' of crack',I4,' is',F15.7,' MPa.')
2161 FORMAT(5X,'-- Error :: Circle inclusion',I3,' is too small!') 
2162 FORMAT(5X,'            Growth length of arc crack > 2/3*2*π*r!') 
2163 FORMAT(5X,'            Error in subroutine Check_Crack_Grows.f!')  
2261 FORMAT(5X,'-- Error :: Return value Error, Point_Status = -1') 
2262 FORMAT(5X,'            in Tool_2_Points_of_Cir_give_Point_L.f!')  
2263 FORMAT(5X,'            Called by subroutine Check_Crack_Grows.f!')  
2271 FORMAT(5X,'-- Error :: Failed to extend arc crack in ') 
2272 FORMAT(5X,'            Tool_Arc_Extend.f called by subroutine')  
2273 FORMAT(5X,'            Check_Crack_Grows.f!') 

!************
! Initialize
!************
Yes_Grow(1:Max_Num_Cr,1:2) = .False.
New_Crack(1:10,1:2,1:2)    = ZR
num_New_Crack              = 0
Na_Cr_num(1:10)            = 0

!*********************
! Some variables used
!*********************
pi_180                     = pi/Con_180
pi_2                       = TWO*pi

!**************************************************
! Critical control variables of the program kernel
!**************************************************
Check_Distance     = Factor_Check_Dis*Ave_Elem_L_Enrich
delta_L_Growth     = Factor_Propagation*Ave_Elem_L_Enrich
L_New_Crack        = Factor_L_Jun_NewCr*Ave_Elem_L_Enrich

! If local encryption is used, the cell feature length before local encryption is adopted, because
! crack propagation still occurs according to the scale of the larger grid before encryption
! (2021-08-22)
if (Key_Local_Mesh_Refine>=1) then
  Check_Distance = Factor_Check_Dis*Ave_Elem_L_Enrich_Unlocalrefined
  delta_L_Growth = Factor_Propagation*Ave_Elem_L_Enrich_Unlocalrefined
  L_New_Crack    = Factor_L_Jun_NewCr*Ave_Elem_L_Enrich_Unlocalrefined
endif


! Given a fixed extension step size (2021-07-28)
if (Propagation_Length> Tol_10)then   
  delta_L_Growth = Propagation_Length
  !Check Propagation_Length, 2021-08-16
  if (delta_L_Growth < 1.0D0*Ave_Elem_L_Enrich)then
      print *, '    ERROR :: *Propagation_Length is too small!'
      print *, '             Please check *Propagation_Length!'     
      print *, '             Error in Check_Crack_Grows.f'    
      call Warning_Message('S',Keywords_Blank)
  endif
  if (delta_L_Growth > 10.0D0*Ave_Elem_L_Enrich)then
      print *, '    ERROR :: *Propagation_Length is too large!'
      print *, '             Please check *Propagation_Length!'     
      print *, '             Error in Check_Crack_Grows.f'                 
      call Warning_Message('S',Keywords_Blank)
  endif          
endif




Ori_delta_L_Growth = delta_L_Growth
delta_L_Growth_Arc = ONEP5*delta_L_Growth

!*******************
! Some preparations
!*******************
! Obtain a model boundary that is connected end to end
EndByEnd_Outline(1:size(Outline,1)) = Outline(:,1)
EndByEnd_Outline(size(Outline,1)+1) = Outline(1,1)   
      
!**********************************
!                               *
!                               *
! Cycle Between Cracks           *
!                               *
!                               *
!**********************************
loop_i_C:do i_C =1,num_Crack    
      ! If the current crack is not allowed to extend, skip it directly.
      if(Cracks_Allow_Propa(i_C) == 0)then  
          goto 100
      endif
      
      ! Read the current crack tip coordinates from the global variable
      cr_Tip(1,1:2) = Cr_First_Tip(i_C,1:2)
      cr_Tip(2,1:2) = Cr_Second_Tip(i_C,1:2)
      
      ! Read the inclination of the crack segment at the crack tip from the global variable
      omega_Tip(1) = Cr_First_Tip_Ori(i_C)
      omega_Tip(2) = Cr_Second_Tip_Ori(i_C)
      
      !----------------------
      ! Two crack tip cycles
      !----------------------
      loop_i_Tip:do i_Tip = 1,2
          ! Initialization: Mark encountered natural fractures
          Yes_HF_met_NF = .False.
          ! If the current crack tip is neither a boundary crack nor a junction point, then calculate its
          ! propagation direction.
          if ((Crack_Tip_Type(i_C,i_Tip) .ne. 1) .and.  &
              (Crack_Tip_Type(i_C,i_Tip) .ne.-1) .and.   &
              (Crack_Tip_Type(i_C,i_Tip) .ne.-2)) then
              !~~~~~~~~~~~~~~~~~~~~~~~~~
              !Angle given by the MCSC.
              !~~~~~~~~~~~~~~~~~~~~~~~~~
              if(CFCP==1)then
                  delta_angle(i_Tip) = TWO*atan(-2*KII(i_C,i_Tip)/KI(i_C,i_Tip)/  &
                            (ONE+sqrt(ONE+EIG*(KII(i_C,i_Tip)/KI(i_C,i_Tip))**2)))
                            
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              !Angle given by the WAMPTSC(weighted average of maximum 
              !principal tensile stress criterion).
              !-------------------
              !In this case, stresses of all Gauss points should be
              !avaliable!
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              elseif(CFCP==2)then
                  !Setting controlling parameter
                  !For option 1, if Key_Ave_Stress==1, SHI's formulation:
                  Search_R = Factor_Ave_R*Ave_Elem_L_Enrich
                  a_Weight = a_Ave_Shi
                  !For option 2, if Key_Ave_Stress==2, Ordinary formulation:
                  l_Ordina = THR*Ave_Elem_L_Enrich
                                                   !                          Theory and Applications
                  call Cal_Weighted_Ave_Stress_of_Point( Key_Ave_Stress,cr_Tip(i_Tip,1:2),Search_R,a_Weight,l_Ordina, &
                     WA_S_xx,WA_S_yy,WA_S_xy)
                  !Calculate the angle of the principal stress.
                  call Tool_Principal_Stresses_2D(WA_S_xx,WA_S_yy,WA_S_xy,WA_S_1,WA_S_3,WA_angle)
                  !GUI output.
                  St_tem = WA_S_1
                  write(*,3001) i_Tip,i_C,St_tem/1.0D6
                  !Possible angle.
                  if(abs(WA_angle-omega_Tip(i_Tip))<=pi/TWO) then
                      Possible_Angle = WA_angle
                  else
                      Possible_Angle = WA_angle + pi
                  endif
                  !Get delta_angle.
                  delta_angle(i_Tip)=Possible_Angle-omega_Tip(i_Tip) 
                  !Checking range (-2*pi to 2*p).
                  if(delta_angle(i_Tip) > pi_2)then
                      delta_angle(i_Tip) = delta_angle(i_Tip)-pi_2
                  elseif(delta_angle(i_Tip) < -pi_2)then
                      delta_angle(i_Tip) = delta_angle(i_Tip)+pi_2
                  endif
                  !Keep abs(delta_angle) < pi.
                  if(delta_angle(i_Tip) > pi)then
                      delta_angle(i_Tip) = delta_angle(i_Tip)-pi_2
                  endif
                  if(delta_angle(i_Tip) < -pi)then
                      delta_angle(i_Tip) = delta_angle(i_Tip)+pi_2
                  endif
                  
              !~~~~~~~~~~~~~~~~~~~~~~
              !Schollmann criterion.
              !~~~~~~~~~~~~~~~~~~~~~~
              elseif(CFCP==3)then
                  
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              !Maximum shear stress criterion.
              ! Maximum shear stress criterion.
              !2024-06-24. NEWFTU2024062401.
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              elseif(CFCP==4)then
                  !Setting controlling parameter
                  !For option 1, if Key_Ave_Stress==1, SHI's formulation:
                  Search_R = Factor_Ave_R*Ave_Elem_L_Enrich
                  a_Weight = a_Ave_Shi
                  !For option 2, if Key_Ave_Stress==2, Ordinary formulation:
                  l_Ordina = THR*Ave_Elem_L_Enrich
                                                   !                          Theory and Applications
                  call Cal_Weighted_Ave_Stress_of_Point( Key_Ave_Stress,cr_Tip(i_Tip,1:2),Search_R,a_Weight,l_Ordina, &
                     WA_S_xx,WA_S_yy,WA_S_xy)
                  !Calculate the angle of the principal stress.
                  call Tool_Principal_Stresses_2D(WA_S_xx,WA_S_yy,WA_S_xy,WA_S_1,WA_S_3,WA_angle)
                  !GUI output.
                  !St_tem = WA_S_1
                  Shear_Stress = (WA_S_1 - WA_S_3)/TWO
                  write(*,3002) i_Tip,i_C,Shear_Stress/1.0D6
                  
                  ! The principal stress at 45 degrees is the direction of the shear stress.
                  WA_angle = WA_angle + pi/FOU 
                  
                  ! Zigzag cracks. For demonstration purposes only. Add a jagged correction angle.
                  if(Key_Sawtooth_Crack_Path==1)then
                      
                      if(mod(iter,2)==0) then
                          Direction_Parameter =  ONE
                      else
                          Direction_Parameter = -ONE
                      endif
                      
                      
                      call Tool_Generate_Random_Number(random_num_for_angle)
                      !random_num_for_angle - 0.5D0
                      
                      ! Random angle range: (20-10, 20+10), that is (10-30)
                      Twist_Angle = (20.0D0+10.0D0*(random_num_for_angle - 0.5D0))*pi/180.0D0*Direction_Parameter
                      
                      
                      WA_angle = WA_angle + Twist_Angle
                  endif
                  
                  !Possible angle.
                  if(abs(WA_angle-omega_Tip(i_Tip))<=pi/TWO) then
                      Possible_Angle = WA_angle
                  else
                      Possible_Angle = WA_angle + pi
                  endif
                  
                  !Get delta_angle.
                  delta_angle(i_Tip)=Possible_Angle-omega_Tip(i_Tip) 
                  
                  !Checking range (-2*pi to 2*p).
                  if(delta_angle(i_Tip) > pi_2)then
                      delta_angle(i_Tip) = delta_angle(i_Tip)-pi_2
                  elseif(delta_angle(i_Tip) < -pi_2)then
                      delta_angle(i_Tip) = delta_angle(i_Tip)+pi_2
                  endif
                  
                  !Keep abs(delta_angle) < pi.
                  if(delta_angle(i_Tip) > pi)then
                      delta_angle(i_Tip) = delta_angle(i_Tip)-pi_2
                  endif
                  if(delta_angle(i_Tip) < -pi)then
                      delta_angle(i_Tip) = delta_angle(i_Tip)+pi_2
                  endif
                  
              endif
              
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! Adjust the angle to prevent it from being too large (<=Prop_Angle_Allowed)
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              if(delta_angle(i_Tip)>Prop_Angle_Alowed*pi_180)then  
                  delta_angle(i_Tip) = Prop_Angle_Alowed*pi_180
              elseif(delta_angle(i_Tip) < -Prop_Angle_Alowed*pi_180)   then
                  delta_angle(i_Tip)=-Prop_Angle_Alowed*pi_180
              end if

              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              !Angle of future crack growth (local coordinates) 
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              theta(i_Tip)=delta_angle(i_Tip) + omega_Tip(i_Tip)   
              
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! Obtain the element number where the crack tip is located
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              call Cal_Ele_Num_by_Coors(cr_Tip(i_Tip,1),cr_Tip(i_Tip,2),Tip_Elem)
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! Material number of the crack tip element
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              mat_num = Elem_Mat(Tip_Elem)
              if(CFCP==1)then
                  if (Material_Type(mat_num) .eq. 1)then
                      KI_Critical  = Material_Para(mat_num,6)
                      ! If the current element is a broken unit, its fracture toughness is 0, 2020-02-28
                      if(Elem_Break(Tip_Elem))then
                          KI_Critical  =ZR
                      endif
                  else
                      KI_Critical  = Material_Para(mat_num,6)                 
                  end if
                  tem1  = delta_angle(i_Tip)/TWO
                  K_tem = cos(tem1)*(KI(i_C,i_Tip)*cos(tem1)**2 - 1.5D0*sin(delta_angle(i_Tip)))
              elseif(CFCP==2)then
                  if (Material_Type(mat_num) .eq. 1  )then
                      St_Critical  = Material_Para(mat_num,5)
                  else
                      St_Critical  = Material_Para(mat_num,5)                 
                  end if
              elseif(CFCP==3)then
                  !For 3D problems only.
              !Maximum shear stress criterion. 2024-06-24. NEWFTU2024062401.
              elseif(CFCP==4)then
                  Shear_Strength = Material_Para_Added(mat_num,12)
              endif
              
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              !checking if the tip grows or not
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              Yes_Tip_Growth = .False.
              if(CFCP==1)then
                  if (K_tem >= KI_Critical) then
                      Yes_Tip_Growth = .True.
                  endif
              elseif(CFCP==2)then
                  if (St_tem >= St_Critical) then
                      Yes_Tip_Growth = .True. 
                  endif
              elseif(CFCP==3)then
                  !For 3D problems only.
              !Maximum shear stress criterion. 2024-06-24. NEWFTU2024062401.
              elseif(CFCP==4)then
                  if (Shear_Stress >= Shear_Strength) then
                      Yes_Tip_Growth = .True. 
                  endif
              endif
              ! If the current crack tip is an arc-shaped crack tip, the crack tip is allowed to extend for
              ! subsequent further calculations and assessments.
              if(Crack_Tip_Type(i_C,i_Tip)==7)then
                  Yes_Tip_Growth = .True. 
              endif
              
              
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! Determination of Load Step for Fatigue Analysis
              ! The step size in the first step should be fixed!!!!!! Otherwise, due to fatigue, the
              ! step size is generally very small, making it difficult to form the first correct
              ! rupture step.
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              if(Key_Static_Fatigue==1 .and. iter >=2 )then
                  k_I = KI(i_C,i_Tip)
                  k_II = KII(i_C,i_Tip)
                  tem_K = (k_I**4 +EIG*k_II**4)**ZP25
                  delta_L_Growth = Fatigue_Paris_C*iter*tem_K**Fatigue_Paris_m
                  ! Information Output
                  write(*,1002) i_Tip,i_C,delta_L_Growth
              endif
              
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! Whether implicit dynamic analysis is extended and how to determine the extension
              ! step size
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              if(Key_Analysis_Type==2)then
                  ! Mark as Do Not Expand First
                  !Yes_Tip_Growth =.False.
                  c_KIc = KIc(mat_num,1)
                  ! Calculate the Rayleigh wave velocity for dynamic analysis
                  c_E  = Material_Para(mat_num,1) 
                  c_v  = Material_Para(mat_num,2)  
                  c_density = density(mat_num)
                  Lambda = c_E*c_v/(ONE+c_v)/(ONE-TWO*c_v)
                  mui    = c_E/TWO/(ONE+c_v)
                  Velocity_long = sqrt((Lambda+TWO*mui)/c_density)
                  Velocity_tran = sqrt(mui/c_density)
                  Velocity_Rayl = Velocity_tran*(0.875-0.2*c_v-0.05*(c_v+0.25)**3)
                  ! Crack Tip Stress Intensity Factor
                  k_I = KI(i_C,i_Tip)
                  k_II = KII(i_C,i_Tip)
                  tem_K = (k_I**4 +EIG*k_II**4)**ZP25
                  
                  !delta_angle_second = TWO*atan(-TWO*k_II/k_I/(ONE+sqrt(ONE+EIG*(k_II/k_I)**2)))
                  
                  ! The denominator cannot be 0. BUGFIX2024091451.
                  tem_value_k_I = k_I
                  if (k_I == ZR)then
                      tem_value_k_I = Tol_20
                  endif
                  delta_angle_second = TWO*atan(-TWO*k_II/tem_value_k_I/(ONE+sqrt(ONE+EIG*(k_II/tem_value_k_I)**2)))
                  
                  tem3  = delta_angle_second/TWO
                  K_tem1 = cos(tem3)*(k_I*cos(tem3)**2 -1.5D0*sin(delta_angle_second))
                  
                  ! The denominator cannot be 0. BUGFIX2024091451.
                  if (K_tem1==ZR) K_tem1=Tol_20
                  tem_delta_L = delt_time_NewMark*Velocity_Rayl*(ONE-(c_KIc/K_tem1)**2)
                  
                  delta_L_Growth = maxval([tem_delta_L,Factor_Propagation*Ave_Elem_L_Enrich])
                  ! However, it must not exceed the maximum step length for crack propagation in dynamic analysis
                  if(delta_L_Growth>=Factor_Prop_Dy*Ave_Elem_L_Enrich)then
                    delta_L_Growth=Factor_Prop_Dy*Ave_Elem_L_Enrich
                  endif
                  ! Information Output
                  !write(*,1002) i_Tip,i_C,tem_delta_L
                  write(*,1012) i_Tip,i_C,delta_L_Growth
              endif
              
              !~~~~~~~~~~~~~~~~~~~~~~~~
              !If the tip grows, then:
              !~~~~~~~~~~~~~~~~~~~~~~~~
              if (Yes_Tip_Growth.eqv..True.) then
                !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                ! If the analysis involves hydraulic fracturing with natural fractures, it is necessary to determine
                ! whether it has already been completed.
                ! It becomes a parasite of natural cracks and is treated differently.
                !-------------------------------
                ! Case 1: If it is a natural bidirectional expansion cemented fracture.
                !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                if((Key_Analysis_Type==3 .or.Key_Analysis_Type==4) .and. &
                    num_Na_Crack>=1 .and. Key_Na_Crack_Type==1 )then
                  ! If the current crack tip of the crack has already completed parasitism, the crack extension
                  ! direction is no longer restricted.
                  if(Cracks_NF_JiS_Stat(i_C,i_Tip)/=-1)then
                    new_tip(i_Tip,1) = cr_Tip(i_Tip,1) + delta_L_Growth*cos(theta(i_Tip))
                    new_tip(i_Tip,2) = cr_Tip(i_Tip,2) + delta_L_Growth*sin(theta(i_Tip))
                  ! If the current crack's current crack tip has not completed parasitism, it will continue to
                  ! propagate along the natural crack.
                  elseif(Cracks_NF_JiS_Stat(i_C,i_Tip)==-1)then
                    c_Na_Cr = Cracks_NF_JiS_Cr_Num(i_C)
                    NF_Tip = Na_Crack_Coor(c_Na_Cr,i_Tip,1:2)
                    ! Tip1 point of the new fracture generated from the HF and NF intersection along the Tip1 direction
                    call Tool_Point_Giv_2P_and_L(  &
                        cr_Tip(i_Tip,1:2),NF_Tip,delta_L_Growth,new_tip(i_Tip,1:2),Statu_Tip,L_C_Tip)
                    ! Check if the parasitism is completed
                    if(Statu_Tip==1 .or. Statu_Tip==0)then            
                        Cracks_NF_JiS_Stat(i_C,i_Tip) = 1
                    endif
                  endif
                !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                ! case2: If it is a natural one-way extension cemented fracture.
                !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                elseif((Key_Analysis_Type==3 .or. Key_Analysis_Type==4) .and. &
                    num_Na_Crack>=1 .and. Key_Na_Crack_Type==2)then
                  ! If the current crack tip of the crack has already completed parasitism, the crack extension
                  ! direction is no longer restricted.
                  if(Cracks_NF_JiS_Stat(i_C,i_Tip)/=-1)then
                    new_tip(i_Tip,1) = cr_Tip(i_Tip,1) + delta_L_Growth*cos(theta(i_Tip))
                    new_tip(i_Tip,2) = cr_Tip(i_Tip,2) + delta_L_Growth*sin(theta(i_Tip))
                  ! If the current crack's current crack tip has not completed parasitism, it will continue to
                  ! propagate along the natural crack.
                  elseif(Cracks_NF_JiS_Stat(i_C,i_Tip)==-1)then
                    c_Na_Cr = Cracks_NF_JiS_Cr_Num(i_C)
                    ! After learning that crack i_C intersects with the cemented crack, does it propagate along the tip
                    ! of crack 1 or the tip of crack 2 of the natural fracture?
                    c_Na_Tip = Cracks_NF_Cement_Tip_Num(i_C,i_Tip)    
                    NF_Tip = Na_Crack_Coor(c_Na_Cr,c_Na_Tip,1:2)
                    ! Tip1 point of the new fracture generated from the HF and NF intersection along the Tip1 direction
                    call Tool_Point_Giv_2P_and_L(      &
                        cr_Tip(i_Tip,1:2),NF_Tip,delta_L_Growth,&
                        new_tip(i_Tip,1:2),Statu_Tip,L_C_Tip)
                    ! Check if the parasitism is completed
                    if(Statu_Tip==1 .or. Statu_Tip==0)then            
                        Cracks_NF_JiS_Stat(i_C,i_Tip) = 1
                    endif
                  endif
                !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                ! case 3: If it is a natural friction crack (the handling code below is exactly the
                ! same as case 1
                ! Same):
                !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                elseif((Key_Analysis_Type==3 .or. Key_Analysis_Type==4) .and. &
                    num_Na_Crack>=1 .and. Key_Na_Crack_Type==3 )then
                  ! If the current crack tip of the crack has already completed parasitism, the crack extension
                  ! direction is no longer restricted.
                  if(Cracks_NF_JiS_Stat(i_C,i_Tip)/=-1)then
                    new_tip(i_Tip,1) = cr_Tip(i_Tip,1) + delta_L_Growth*cos(theta(i_Tip))
                    new_tip(i_Tip,2) = cr_Tip(i_Tip,2) + delta_L_Growth*sin(theta(i_Tip))
                  ! If the current crack's current crack tip has not completed parasitism, it will continue to
                  ! propagate along the natural crack.
                  elseif(Cracks_NF_JiS_Stat(i_C,i_Tip)==-1)then
                    c_Na_Cr = Cracks_NF_JiS_Cr_Num(i_C)
                    
                    NF_Tip = Na_Crack_Coor(c_Na_Cr,i_Tip,1:2)
                    ! Tip1 point of the new fracture generated from the HF and NF intersection along the Tip1 direction
                    !--------------------------------------------------------------------------------
                    ! TDL -  TDL  
                    ! Note, here at the intersection of HF and friction-type NF, there is a major
                    ! bug, because friction-type natural fractures already exist, so it cannot
                    ! follow the original path
                    ! The natural cracks gradually expanded, and this treatment method is completely
                    ! incorrect; cr_Tip coincides with the endpoints of the natural cracks.
                    ! Added description time: 2018-10-31
                    !--------------------------------------------------------------------------------
                    call Tool_Point_Giv_2P_and_L(     &
                        cr_Tip(i_Tip,1:2),NF_Tip,delta_L_Growth,&
                        new_tip(i_Tip,1:2),Statu_Tip,L_C_Tip)
                    ! Check if the parasitism is completed
                    if(Statu_Tip==1 .or. Statu_Tip==0)then    
                        Cracks_NF_JiS_Stat(i_C,i_Tip) = 1
                    endif
                  endif
                !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                ! Generally, there are no issues related to natural cracks.
                !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                else
                  new_tip(i_Tip,1) = cr_Tip(i_Tip,1) + delta_L_Growth*cos(theta(i_Tip))
                  new_tip(i_Tip,2) = cr_Tip(i_Tip,2) + delta_L_Growth*sin(theta(i_Tip)) 
                endif
              else
                  ! If the load control mode is 5, CFCP = 2, and there are polygonal inclusion cracks, then: although
                  ! The maximum tensile stress criterion is not met, but the inclusion interface may still propagate,
                  ! so it is necessary
                  ! Continue running downward -- goto 88
                  !2018-01-11
                  if(num_Poly_Incl/=0 .and. CFCP==2 .and. Key_Force_Control == 5)then
                      goto 88
                  ! Otherwise, proceed to handle the next crack tip
                  else
                      goto 99
                  endif
              end if

              
              88  continue                
              
              
              
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! Handling the Intersection and Separation of Cracks and Polygonal Inclusions
            !    -----------------------------------------------------------------------------------
            ! (1) Crack_Tip_Type = -8 indicates a crack tip that has just extended to the polygonal
            ! inclusion interface.
            ! Crack_Tip_Type = -9 indicates a crack tip that has extended along the polygonal
            ! inclusion interface
            ! Note: Once an 8-shaped crack tip has formed, it should extend at least one step along
            ! the inclusion interface before
            ! Mixed interface
            !   ------------------------------------------------------------------------------------
            !   Added on 2018-01-07
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              if (num_Poly_Incl/=0)then
                  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                  ! If it is a standard notch tip
                  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                  if (Crack_Tip_Type(i_C,i_Tip)==0)then
                      ! Check whether the current crack tip coordinates are within an inclusion (a closed-form polygonal
                      ! inclusion)
                      do i_Incl = 1,num_Poly_Incl
                        n_Vertex = Poly_Inclu_Edges_Num(i_Incl)
                        call Tool_Yes_In_Poly(new_tip(i_Tip,1),new_tip(i_Tip,2),&
                         Poly_Incl_Coor_x_Cl(i_Incl,1:n_Vertex+1),&
                         Poly_Incl_Coor_y_Cl(i_Incl,1:n_Vertex+1),n_Vertex+1,Yes_in_Poly_Inc)
                        ! If present, check whether the maximum principal stress exceeds the tensile strength of the
                        ! inclusion.
                        ! If it exceeds, the requirement is met; if it does not, adjust the crack tip to fit.
                        ! On the surface, mark the crack tip as the crack tip on the polygonal inclusion interface:
                        !Crack_Tip_Type=8
                        if(Yes_in_Poly_Inc)then
                          if (St_tem < Material_Para(Poly_Inclu_Mat_Num(i_Incl),5))then
                            !.....................................................................
                            ! Correct the crack tip back onto the interface
                            ! Note: In subsequent expansions, priority will be given to extending
                            ! along the interface.
                            !.....................................................................
                            ! (1) Find the intersecting edges and the corresponding intersection points
                            do i_S=1,n_Vertex
                                Side_A(1) = Poly_Incl_Coor_x_Cl(i_Incl,i_S)
                                Side_A(2) = Poly_Incl_Coor_y_Cl(i_Incl,i_S)     
                                Side_B(1) = Poly_Incl_Coor_x_Cl(i_Incl,i_S+1)
                                Side_B(2) = Poly_Incl_Coor_y_Cl(i_Incl,i_S+1)          
                                call Tool_Intersection(cr_Tip(i_Tip,1:2),new_tip(i_Tip,1:2), &
                                  Side_A,Side_B,c_inter_X,c_inter_Y,c_Yes_Cross)
                                if (c_Yes_Cross)then
                                    exit
                                endif
                            enddo
                            ! Correct the crack tip back onto the interface
                            new_tip(i_Tip,1) = c_inter_X
                            new_tip(i_Tip,2) = c_inter_Y
                            ! Fracture Tip Type Label
                            Crack_Tip_Type(i_C,i_Tip)=-8
                            ! Label the corresponding inclusion number and the edge number of the inclusion
                            Crack_Tip_Ploy_Inc_Info(i_C,i_Tip,1)= i_Incl
                            Crack_Tip_Ploy_Inc_Info(i_C,i_Tip,2)= i_S   
                          endif
                          exit
                        endif
                      enddo
                  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                  ! If it is a type 8 crack tip, the crack needs to be verified to extend at least
                  ! one step along the inclusion boundary.
                  ! After extending one step, it becomes a type 9 split tip
                  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                  elseif (Crack_Tip_Type(i_C,i_Tip)==-8)then
                      c_Incl = Crack_Tip_Ploy_Inc_Info(i_C,i_Tip,1)
                      c_Side = Crack_Tip_Ploy_Inc_Info(i_C,i_Tip,2)
                      Side_A(1)=Poly_Incl_Coor_x_Cl(c_Incl,c_Side)
                      Side_A(2)=Poly_Incl_Coor_y_Cl(c_Incl,c_Side)     
                      Side_B(1)=Poly_Incl_Coor_x_Cl(c_Incl,c_Side+1)
                      Side_B(2)=Poly_Incl_Coor_y_Cl(c_Incl,c_Side+1)   
                      ! There are two possible crack tips, and you need to choose the one with the obtuse angle. The
                      ! principle is in my notes V5-p143.
                      cc_Line_AB(1,1:2) = cr_Tip(i_Tip,1:2)
                      cc_Line_AB(2,1:2) = Side_A
                      call Tool_Shorten_or_Extend_Line(cc_Line_AB,-delta_L_Growth,'A',c_new_Line_AB,Possoble_P_1)
                      cc_Line_AB(2,1:2) = Side_B
                      call Tool_Shorten_or_Extend_Line(cc_Line_AB,-delta_L_Growth,'A',c_new_Line_AB,Possoble_P_2)   
                      ! Obtain O-point coordinates
                      if(i_Tip==1) then
                          c_O(1:2) = Crack_Coor(i_C,2,1:2)
                      elseif(i_Tip==2)then
                          num_pt = Each_Cr_Poi_Num(i_C)
                          c_O(1:2) =  Crack_Coor(i_C,num_pt-1,1:2)
                      endif
                      L_P_1_O=Tool_Function_2Point_Dis(c_O,Possoble_P_1)     
                      L_P_2_O=Tool_Function_2Point_Dis(c_O,Possoble_P_2)  
                      ! The one that determines the obtuse angle based on distance; the principle can be seen in my notes
                      ! V5-p143.
                      if(L_P_1_O .GE. L_P_2_O) then
                          new_tip(i_Tip,1:2) = Possoble_P_1
                      else
                          new_tip(i_Tip,1:2) = Possoble_P_2
                      endif
                      ! Determine whether the new crack tip is still on the original edge. If it is, it becomes a type -9
                      ! crack tip; if not, it is corrected to a regular crack tip.
                      call Tool_Yes_On_Line(new_tip(i_Tip,1),new_tip(i_Tip,2),Side_A,Side_B,c_Yes_ON)
                      if(c_Yes_ON)then
                          Crack_Tip_Type(i_C,i_Tip)=-9
                      else
                          Crack_Tip_Type(i_C,i_Tip)=0
                      endif
                      ! Label the corresponding inclusion number and the edge number of the inclusion
                      ! No changes needed
                      !new_tip(i_Tip,1) = cr_Tip(i_Tip,1:2)
                  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                  ! If it is a type 9 crack tip, it is necessary to determine whether it will
                  ! continue to propagate along the interface or detach from the interface.
                  ! Extending towards the matrix? See the idea on V5-P143
                  ! PhiPsi's disposal strategy is as follows: first, determine whether the interface
                  ! meets the extension conditions, that is, case (1) if the boundary
                  ! If the tensile strength on the surface is less than the tensile stress at the
                  ! interface, it will preferentially propagate along the interface.
                  ! If not satisfied, then determine in case (2) whether the crack tip calculated
                  ! according to the maximum principal stress criterion is in the direction.
                  ! Within the substrate, if so, it will extend into the substrate. If neither of
                  ! the above two conditions is met, the expansion stops.
                  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                  elseif (Crack_Tip_Type(i_C,i_Tip)==-9)then
                      c_Incl = Crack_Tip_Ploy_Inc_Info(i_C,i_Tip,1)
                      c_Side = Crack_Tip_Ploy_Inc_Info(i_C,i_Tip,2)
                      Side_A(1)=Poly_Incl_Coor_x_Cl(c_Incl,c_Side)
                      Side_A(2)=Poly_Incl_Coor_y_Cl(c_Incl,c_Side)     
                      Side_B(1)=Poly_Incl_Coor_x_Cl(c_Incl,c_Side+1)
                      Side_B(2)=Poly_Incl_Coor_y_Cl(c_Incl,c_Side+1)   
                      !................................................
                      ! First, determine whether case (1) is satisfied
                      !................................................
                      ! There is one possible crack tip, which is Possoble_P_1
                      ! Obtain the coordinates of point I
                      if(i_Tip==1) then
                          c_I(1:2) = Crack_Coor(i_C,2,1:2)
                      elseif(i_Tip==2)then
                          num_pt = Each_Cr_Poi_Num(i_C)
                          c_I(1:2) =  Crack_Coor(i_C,num_pt-1,1:2)
                      endif                          
                      cc_Line_AB(1,1:2) = c_I
                      cc_Line_AB(2,1:2) = cr_Tip(i_Tip,1:2)
                      call Tool_Shorten_or_Extend_Line(cc_Line_AB,delta_L_Growth,'B',c_new_Line_AB,Possoble_P_1)
                      ! Calculate the weighted stress at point P1
                      call Cal_Weighted_Ave_Stress_of_Point(Key_Ave_Stress, &
                          Possoble_P_1,Search_R,a_Weight,l_Ordina, &
                          WA_S_xx,WA_S_yy,WA_S_xy)
                      ! Calculate the inclination of I_P1, denoted as c_Alpha
                      c_Alpha= atan2(Possoble_P_1(2)-c_I(2),Possoble_P_1(1)-c_I(1))
                      ! Calculate the normal stress and shear stress.
                      call Tool_Normal_and_Shear_Stresses_2D(WA_S_xx,WA_S_yy,WA_S_xy,c_Alpha, &
                                 WA_S_normal,WA_S_shear)   
                      ! If the normal stress on the interface is greater than the tensile strength of the interface, then
                      ! case 1 is satisfied.
                      if(WA_S_normal>=Material_Interface(1)*0.99D0)then
                          new_tip(i_Tip,1:2) = Possoble_P_1
                          ! Determine whether P1 is still on the inclusion boundary. If it is, it remains a -9 type crack tip;
                          ! if not, correct it to a regular crack tip.
                          call Tool_Yes_On_Line(new_tip(i_Tip,1),new_tip(i_Tip,2),Side_A,Side_B,c_Yes_ON)
                          if(c_Yes_ON)then
                              Crack_Tip_Type(i_C,i_Tip)=-9
                          else
                              Crack_Tip_Type(i_C,i_Tip)=0
                          endif
                          goto 1542
                      endif 
                      !.....................................
                      ! Determine whether it meets case (2)
                      !.....................................
                      ! Check whether the current new crack tip coordinates are within the matrix, that is, not inside any
                      ! inclusion (a closed-form polygonal inclusion).
                      do i_Incl = 1,num_Poly_Incl
                        n_Vertex = Poly_Inclu_Edges_Num(i_Incl)
                        call Tool_Yes_In_Poly( new_tip(i_Tip,1),new_tip(i_Tip,2), &
                          Poly_Incl_Coor_x_Cl(i_Incl,1:n_Vertex+1),Poly_Incl_Coor_y_Cl(i_Incl,1:n_Vertex+1), &
                          n_Vertex+1,Yes_in_Poly_Inc)
                        if (Yes_in_Poly_Inc) then
                            ! If it is in a certain inclusion, then it means that case (2) is impossible, jump to 1542
                            goto 1542
                        endif
                      enddo
                      ! Reaching this point means there is no diversion, indicating that the potential crack tip is within
                      ! the matrix. At this point, the crack tip can be used for the calculation.
                      goto 1542
                      !.......................................................................
                      ! If neither case (1) nor case (2) is satisfied, the crack tip does not
                      ! propagate and remains a Type-9 crack.
                      !.......................................................................
                      new_tip(i_Tip,1:2) = c_I(1:2)
                      Crack_Tip_Type(i_C,i_Tip)=-9
                      
                      1542 continue    

                  endif
              endif
              
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! Check whether the new crack tip coordinates exceed the model range, and if they do
              ! Then the crack tip does not propagate and will not propagate afterward!
              ! If not exceeded, then Yes_Grow(i_C, i_Tip) = .True.
              !-----
              ! In addition, if a rupture zone is defined, check whether it exceeds the rupture zone.
              ! The scope, added on 2016-03-28
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              if (Key_Fracture_Zone==0) then
                  Call Tool_Yes_In_Poly(new_tip(i_Tip,1),new_tip(i_Tip,2), &
                      Coor(EndByEnd_Outline,1),Coor(EndByEnd_Outline,2), &
                      size(EndByEnd_Outline),Yes_In_Model)
                  if (Yes_In_Model) then              
                      Yes_Grow(i_C,i_Tip) = .True.
                  else
                      Crack_Tip_Type(i_C,i_Tip) = -1
                      write(*,1003) i_Tip,i_C
                      write(*,1004) 
                      !Yes_Grow(i_C,i_Tip) = .True.
                  end if 
              elseif (Key_Fracture_Zone==1) then
                  Yes_In_Frac_zone = .True.
                  if (new_tip(i_Tip,1) < Frac_Zone_MinX .or. new_tip(i_Tip,1) > Frac_Zone_MaxX .or. &
                      new_tip(i_Tip,2) < Frac_Zone_MinY .or.new_tip(i_Tip,2) > Frac_Zone_MaxY) then
                      Yes_In_Frac_zone = .False.
                  endif
                  ! The new crack tip coordinates are within the rupture zone
                  if (Yes_In_Frac_zone) then              
                      Yes_Grow(i_C,i_Tip) = .True.
                  ! The new crack tip coordinates are not within the fracture zone
                  else
                      Crack_Tip_Type(i_C,i_Tip) = -1
                      write(*,1003) i_Tip,i_C
                      write(*,1005) 
                  end if  
              endif
              
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! Check whether the new crack tip is on the element boundary; if it is, it must be
              ! corrected (move slightly toward the lower-left corner).
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              call Cal_Ele_Num_by_Coors_YesInOn(new_tip(i_Tip,1),new_tip(i_Tip,2),c_Ele,Yes_In,Yes_on_Ele_Edge)
              Modifi_Factor = 1.0D-5
              if (Yes_on_Ele_Edge .eqv. .True.) then
                  new_tip(i_Tip,1) = new_tip(i_Tip,1) -Modifi_Factor*Ave_Elem_L_Enrich
                  new_tip(i_Tip,2) = new_tip(i_Tip,2) -Modifi_Factor*Ave_Elem_L_Enrich
                  write(*,2015)
              end if
              
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! Check whether the new crack segment passes exactly through a certain node, and if it
              ! does, it must be corrected.
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              Yes_Node_OnSeg = .False.
              Modifi_Factor = 1.0D-4
              ! Loop through all nodes (Note: this scheme has room for improvement, only the nodes near the cracks
              ! need to be checked, TDL)
              do i_Node = 1,Num_Node
                  ! Current node coordinates
                  c_Node_x = Coor(i_Node,1)
                  c_Node_y = Coor(i_Node,2)
                  ! Check if the current node is on the current crack segment
                  call Tool_Yes_On_Line(c_Node_x,c_Node_y,cr_Tip(i_Tip,1:2),new_tip(i_Tip,1:2),Yes_Node_OnSeg)
                  ! If on the current crack segment, then correct the left endpoint.
                  if (Yes_Node_OnSeg.eqv..True.) then
                      new_tip(i_Tip,1) = new_tip(i_Tip,1) -Modifi_Factor*Ave_Elem_L_Enrich
                      new_tip(i_Tip,2) = new_tip(i_Tip,2) -TWO*Modifi_Factor*Ave_Elem_L_Enrich
                      write(*,2016)
                      exit
                  end if
              end do
              
              


              
              
              
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! Handling the Issue of Crack Intersections with Circular Inclusions and the Continued
            ! Propagation of Arc-shaped Crack Tips
            !    ------------------------------------------------------------------------------------
            ! (1) If there are spherical inclusions and the crack tip is a standard crack tip
            ! (Crack_Tip_Type=0),
            ! Then check whether the new crack tip coordinates are inside the circular inclusion; if
            ! so, a redo is required.
            ! Judgment: Does it generate arc-shaped cracks along the inclusion interface, or does it
            ! penetrate the inclusion directly?
            ! If an arc crack is generated, add 2 new crack coordinate points.
            ! (2) If an arc-shaped crack tip has already formed (or was originally present)
            ! (Crack_Tip_Type=7)
            !    ----------------------------------------------------------------------------------
            ! Note: If the new crack tip is close to the inclusion (less than one extension step),
            ! and the old and new crack tips
            ! If the extension of the connecting line intersects with this inclusion, it is also
            ! necessary to determine whether
            ! Generate arched cracks; the purpose of handling it this way is to prevent them from
            ! being too close together.
            ! A unit has a crack segment
            !   -----------------------------------------------------------------------------------
            ! added on 2017-07-20, corresponding note description: V5-58_P117
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              !:::::::::::::::::::::::
              ! First, initialize the relevant variables
              !:::::::::::::::::::::::
              ! Generate new arc-shaped crack marks?
              Yes_New_Arc_Cr    = .False.   
              ! Whether to update the coordinates (crack tip) of the curved crack
              Yes_Update_Arc_Cr = .False.
              if (num_Circ_Incl/=0) then 
                  ! Temporarily store data related to the arcs of curved cracks (10 data points)
                  c_Arc_Cr(1:11)    = ZR
                  !::::::::::::::::::::::::::::::::::::::::::
                  ! If it is a standard crack tip (an arc-shaped crack tip has not yet formed)
                  !::::::::::::::::::::::::::::::::::::::::::
                  if (Crack_Tip_Type(i_C,i_Tip)==0)then
                      Yes_In_Incl = .False.
                      ! Cycle among various interspersed intervals
                      do i_Incl =1,num_Circ_Incl
                          c_Incl_x =  Circ_Inclu_Coor(i_Incl,1)
                          c_Incl_y =  Circ_Inclu_Coor(i_Incl,2)
                          c_Incl_r =  Circ_Inclu_Coor(i_Incl,3) 
                          c_Dis  =Tool_Function_2Point_Dis([new_tip(i_Tip,1),new_tip(i_Tip,2)],[c_Incl_x,c_Incl_y])
                          if(c_Dis<=c_Incl_r)then
                            Yes_In_Incl = .True.
                            !//////////////////////////////////////////////////////////////////////
                            ! Calculate the coordinates of the intersection points between the old
                            ! and new crack tip line segments and inclusions
                            !//////////////////////////////////////////////////////////////////////
                            call Tool_Intersection_Line_and_Circle(c_Incl_x,c_Incl_y,c_Incl_r, &
                               cr_Tip(i_Tip,1:2),new_tip(i_Tip,1:2),c_num_Inter,c_State,c_Inter_P)
                            !////////////////////////////////////////////////////////////////////////
                            ! Determine whether the arc crack propagation step (1.5 times the linear
                            ! propagation step) is greater than the inclusion diameter
                            ! If the inclusion is too small, then the subsequent potential crack
                            ! point search algorithm (finding the intersection of two circles)
                            ! About to expire
                            !////////////////////////////////////////////////////////////////////////
                            if(delta_L_Growth_Arc >= TWO*c_Incl_r) then
                                write (*,2161) i_Incl
                                write (*,2162)
                                write (*,2163)
                                call Warning_Message('S',Keywords_Blank)
                            end if
                            !/////////////////////////////////////////////////////////////////////
                            ! Calculating potential crack points: P1 and P2 in the notebook image
                            ! (V5-P117)
                            !/////////////////////////////////////////////////////////////////////
                            call Tool_2_Points_of_Cir_give_Point_L(c_Incl_x,c_Incl_y,c_Incl_r, &
                                   c_Inter_P(1,1:2),delta_L_Growth_Arc,c_Point_Status,c_P1_P2) 
                            ! If two points are successfully found
                            if(c_Point_Status ==1 .and. (sum(c_P1_P2(2,1:2))>Tol_11)) then                                  
                              !----------------------------
                              !         OPTION 1
                              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                              ! Determine the points to choose based on the distance.
                              ! Algorithm see Notes V5-58_P117
                              !----------------------------
                              L_AP1 = Tool_Function_2Point_Dis(cr_Tip(i_Tip,1:2),c_P1_P2(1,1:2))
                              L_AP2 = Tool_Function_2Point_Dis(cr_Tip(i_Tip,1:2),c_P1_P2(2,1:2))
                              if (L_AP1 >= L_AP2)then
                                  Want_Point(1:2) = c_P1_P2(1,1:2)
                              else
                                  Want_Point(1:2) = c_P1_P2(2,1:2)
                              endif
                            ! An unexpected unknown error
                            elseif(c_Point_Status ==-1)then
                                write (*,2261)
                                write (*,2262)
                                write (*,2263)
                                call Warning_Message('S',Keywords_Blank)
                            endif
                            
                            !/////////////////////////////////////////////////////////////////////
                            ! Calculate and temporarily store the related data of the arcs of the
                            ! arc-shaped cracks (10 pieces of data):
                            !x,y,r,Radian_Start,Radian_End,Radian,Point_Start_x,
                            !Point_Start_y,Point_End_x,Point_End_y
                            !/////////////////////////////////////////////////////////////////////
                            ! Update the local arc crack marking in the subroutine
                            Yes_New_Arc_Cr    = .True.
                            ! Updated global PhiPsi arc fracture markers, Debugging on 2017-07-22
                            Yes_Arc_Crack     = .True.
                            c_Arc_Cr(1:2)  = [c_Incl_x,c_Incl_y]
                            c_Arc_Cr(4)    =  c_Incl_r
                            ! Determine whether it is a clockwise arc or a counterclockwise arc. Counterclockwise is 1,
                            ! clockwise is -1.
                            ! Algorithm: First, calculate the radian angle corresponding to the arc counterclockwise. If the
                            ! radian angle is <= 180 degrees,
                            ! Then accept it; otherwise, take it as clockwise.
                            c_Arc_Direction = ONE
                            ! Calculate other data of the arc
                            call Tool_Arc_r_and_Radian_Given_Coors(c_Arc_Direction, &
                                c_Inter_P(1,1:2),Want_Point,[c_Incl_x,c_Incl_y], &
                                c_Yes_feasible,c_Arc_r,c_Arc_Radian_S,c_Arc_Radian_E,c_Arc_Radian)
                            if (c_Arc_Radian > Con_180)then
                              c_Arc_Direction = -ONE
                              call Tool_Arc_r_and_Radian_Given_Coors(c_Arc_Direction,c_Inter_P(1,1:2),Want_Point, &
                                  [c_Incl_x,c_Incl_y],c_Yes_feasible,c_Arc_r,c_Arc_Radian_S, &
                                  c_Arc_Radian_E,c_Arc_Radian) 
                            endif
                            ! Save the other data of the arc to the temporary variable c_Arc_Cr
                            c_Arc_Cr(3)    = c_Arc_Direction
                            c_Arc_Cr(5)    = c_Arc_Radian_S
                            c_Arc_Cr(6)    = c_Arc_Radian_E
                            c_Arc_Cr(7)    = c_Arc_Radian
                            c_Arc_Cr(8)    = c_Inter_P(1,1)
                            c_Arc_Cr(9)    = c_Inter_P(1,2)
                            c_Arc_Cr(10)   = Want_Point(1) 
                            c_Arc_Cr(11)   = Want_Point(2) 
                            exit
                          endif
                      enddo 
                  !::::::::::::::::::::::::::::::::::::::::::::::::
                  ! If an arc-shaped crack tip has already formed (or was initially present)
                  !  -------------------------------------------
                  ! There are three possible situations:
                  ! (1) Stop expanding;
                  ! (2) Continue to extend along the inclusion interface;
                  ! (3) Extend into the matrix;
                  ! (4) Expand into interleaving (temporarily unavailable).
                  !::::::::::::::::::::::::::::::::::::::::::::::::
                  elseif (Crack_Tip_Type(i_C,i_Tip)==7)then
                      Yes_Grow(i_C,i_Tip) = .True.
                      !////////////////////////////////////////////////////////
                      !                       CASE 1
                      !         -----------------------------------
                      ! If the mixed-mode interfacial crack propagation criterion is not met: the stress intensity factor
                      ! is less than the interface KIc.
                      ! (or if the tensile strength is less than the interface St), the crack tip stops propagating
                      !////////////////////////////////////////////////////////
                      if(CFCP==1)then
                          if(K_tem < Material_Interface(2))then
                              Yes_Grow(i_C,i_Tip) = .False.
                          endif
                      elseif(CFCP==2)then
                          if(St_tem < Material_Interface(1))then
                              Yes_Grow(i_C,i_Tip) = .False.
                          endif
                      endif
                      
                      !/////////////////////////////////////////////////////
                      !                       CASE 2
                      !         -----------------------------------
                      ! If the current arc's angle has already exceeded 240 degrees, stop expanding it as well.
                      !/////////////////////////////////////////////////////
                      if(i_Tip ==1)then
                          c_Arc_Radian = Arc_Crack_Coor(i_C,1,7)
                      elseif(i_Tip ==2)then
                          c_Arc_Radian = Arc_Crack_Coor(i_C,Each_Cr_Poi_Num(i_C)-1,7)
                      endif
                      if(c_Arc_Radian >=Con_240)then
                          Yes_Grow(i_C,i_Tip) = .False.
                      endif
                      
                      !/////////////////////////////////////////////////////////
                      !                       CASE 3
                      !         -----------------------------------
                      ! If the interfacial crack propagation criterion is met, and the crack deflection angle is greater
                      ! than 30 degrees,
                      ! And if the new crack tip is inside any inclusion, the propagation is stopped.
                      !/////////////////////////////////////////////////////////
                      if(CFCP==1)then
                          if(K_tem >= Material_Interface(2))then
                            if (abs(delta_angle(i_Tip))>Con_30)then
                                ! Check whether the potential new crack tip coordinates are within any inclusion (in fact, this can
                                ! be further optimized to correspond to the inclusion associated with the previously formed curved
                                ! crack).
                                call Tool_Yes_Point_in_Inclusions(new_tip(i_Tip,1),&
                                         new_tip(i_Tip,2),c_Yes_in_Incl,c_Incl_Num)
                                if(c_Yes_in_Incl)then
                                   Yes_Grow(i_C,i_Tip) = .False.
                                endif
                                !new_tip(i_Tip,1:2)
                            endif
                          endif
                      elseif(CFCP==2)then
                          if(St_tem >= Material_Interface(1))then
                              Yes_Grow(i_C,i_Tip) = .False.
                          endif
                      endif
                      
                      !//////////////////////////////////////////////////////////
                      !                       CASE 4
                      !         -----------------------------------
                      ! If the stress intensity factor is greater than the interface K_Ic, and the crack propagation
                      ! deflection angle is greater than 30 degrees,
                      ! If the new crack tip is outside the inclusion, it will extend toward the matrix.
                      !//////////////////////////////////////////////////////////
                      if(CFCP==1)then
                          if(K_tem >= Material_Interface(2))then
                            if (abs(delta_angle(i_Tip))>Con_30)then
                                ! Check whether the potential new crack tip coordinates are within any inclusion (in fact, this can
                                ! be further optimized to correspond to the inclusion associated with the previously formed curved
                                ! crack).
                                call Tool_Yes_Point_in_Inclusions(new_tip(i_Tip,1),new_tip(i_Tip,2),&
                                                                  c_Yes_in_Incl,c_Incl_Num)
                                ! If the new crack tip is not within any inclusion
                                if(c_Yes_in_Incl .eqv. .False.)then
                                   Yes_Grow(i_C,i_Tip) = .True.
                                endif
                                ! The crack tip type has been changed from a curved tip to a standard tip.
                                Crack_Tip_Type(i_C,i_Tip) = 0
                                !new_tip(i_Tip,1:2)
                            endif
                          endif
                      elseif(CFCP==2)then
                          if(St_tem >= Material_Interface(1))then
                              Yes_Grow(i_C,i_Tip) = .False.
                          endif
                      endif

                      !///////////////////////////////////////////////////////////
                      !                       CASE 5
                      !         -----------------------------------
                      ! If the interfacial crack propagation criterion is met, and the crack deflection angle is less than
                      ! 30 degrees,
                      ! Then continue to expand along the inclusion interface: Yes_Update_Arc_Cr = .False. and update
                      ! Arc-shaped crack tip (number of coordinate points does not increase)
                      !///////////////////////////////////////////////////////////
                      if(CFCP==1)then
                          if(K_tem >= Material_Interface(2))then
                            ! Obtain the length factor of the newly developed arc-shaped crack (subtract one ten-thousandth at
                            ! each step)
                            Factor_L = ONE-iter*ZPZZZ1
                            if (abs(delta_angle(i_Tip))<Con_30)then
                                Yes_Grow(i_C,i_Tip) = .True.
                                ! Activate update arc marker
                                Yes_Update_Arc_Cr = .True.
                                ! New arc (arc direction remains the same, length is increased (length is random, with a variation
                                ! range of 1%))
                                if(i_Tip==1)then
                                  call Tool_Arc_Extend_or_Shorten(delta_L_Growth_Arc*Factor_L,  &
                                       Arc_Crack_Coor(i_C,1,1:11),c_Arc_Cr,Extend_Status) 
                                  ! If there is no problem with updating the arc, then:
                                  if(Extend_Status)then
                                  !nothing
                                  print *,c_Arc_Cr
                                  ! If there is a problem updating the arc, terminate the program:
                                  else
                                      write (*,2271)
                                      write (*,2272)
                                      write (*,2273)
                                      call Warning_Message('S',Keywords_Blank)
                                  endif
                                elseif(i_Tip==2)then
                                  call Tool_Arc_Extend_or_Shorten(&
                                      delta_L_Growth_Arc*Factor_L, &
                                      Arc_Crack_Coor(i_C,Each_Cr_Poi_Num(i_C)-1,1:11),&
                                      c_Arc_Cr,Extend_Status)
                      
                                  ! If there is no problem with updating the arc, then:
                                  if(Extend_Status)then
                                  !nothing
                                  ! If there is a problem updating the arc, terminate the program:
                                  else
                                      write (*,2271)
                                      write (*,2272)
                                      write (*,2273)
                                      call Warning_Message('S',Keywords_Blank)
                                  endif
                                endif
                            endif
                          endif
                      elseif(CFCP==2)then
                          if(St_tem >= Material_Interface(1))then
                            ! Obtain the length factor of the newly developed arc-shaped crack (subtract one ten-thousandth at
                            ! each step)
                            Factor_L = ONE-iter*ZPZZZ1
                          endif
                      endif
                      
                  endif
              endif

              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! If it is the final rupture step, no further expansion is allowed.
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              if(ifra==Num_Frac)then
                   Yes_Grow(i_C,i_Tip) = .False.
              endif
              
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! Update crack coordinates and crack segments (divided into fixed step size and
              ! automatically calculated step size cases)
              ! (1) Fixed step size: No need to check whether the element contains independent crack
              ! segments; cracks are added directly.
              ! Coordinate point;
              ! (2) Calculating step size: It is necessary to determine whether an element contains
              ! independent crack segments, and if it does, then
              ! Directly adjust the original crack tip position to the new position; if not
              ! included, add a new crack coordinate point.
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              
              ! Then all are forcibly changed to (2), that is, it is necessary to determine whether a unit
              ! contains independent crack segments, and if it does,
              ! Directly adjust the original crack tip position to the new position; if not included, add a new
              ! crack coordinate point.
              ! Modified Date: 2017-04-16
              Key_Propa_Type = 2
              ! Directly adjust the original crack tip position to the new position; if not included, add a new
              ! crack coordinate point.
              if (Yes_Grow(i_C,i_Tip) .eqv. .True.) then
                  !/////////////////////////////
                  !Fixed propagation increment.
                  !/////////////////////////////
                  if (Key_Propa_Type==1)then
                    ! Adjust first
                    if (i_Tip==1) then
                      ! Number of previous crack points
                      Old_num_Poi = Each_Cr_Poi_Num(i_C)
                      ! Temporary variable: stores the previous crack coordinates
                      tem_Crack_Coor(1:Old_num_Poi,1:2) = Crack_Coor(i_C,1:Old_num_Poi,1:2)
                      ! Update the number of crack coordinate points
                      Each_Cr_Poi_Num(i_C) = Each_Cr_Poi_Num(i_C) +1
                      ! Add the newly created crack tip
                      Crack_Coor(i_C,1,1:2) = new_tip(i_Tip,1:2)   
                      ! Move the coordinates of the crack at the back one position backward
                      Crack_Coor(i_C,2:Old_num_Poi+1,1:2) = tem_Crack_Coor(1:Old_num_Poi,1:2)
                    elseif(i_Tip==2) then
                      Old_num_Poi = Each_Cr_Poi_Num(i_C)
                      ! Update the number of crack coordinate points
                      Each_Cr_Poi_Num(i_C) = Each_Cr_Poi_Num(i_C) +1
                      ! Add the newly created crack tip
                      Crack_Coor(i_C,Old_num_Poi+1,1:2) = new_tip(i_Tip,1:2)
                    end if
                    
                  !/////////////////////////////////////////////////////////////////////////
                  !Caculated propagation increment, the delta_L may be relatively small.
                  ! Make sure the newly generated crack segments are long enough! >= 0.99 *
                  ! Ori_delta_L_Growth
                  !/////////////////////////////////////////////////////////////////////////
                  elseif(Key_Propa_Type==2)then
                    Yes_In_OneEl = .False.
                    call Cal_Ele_Num_by_Coors(new_tip(i_Tip,1),new_tip(i_Tip,2),New_OUT_Elem)
                    !IIIIIIIIIIIIIIIIIIIIIIIIIIIIII
                    ! If it is the number 1 tip
                    !IIIIIIIIIIIIIIIIIIIIIIIIIIIIII
                    if (i_Tip==1) then
                      !::::::::::::::::::::::
                      ! Number of previous crack points
                      !::::::::::::::::::::::
                      Old_num_Poi = Each_Cr_Poi_Num(i_C)
                      Old_Tip(1:2) = Crack_Coor(i_C,1,1:2)
                      Old_Old_Tip(1:2) = Crack_Coor(i_C,2,1:2)
                      c_L = sqrt((Old_Old_Tip(2)-new_tip(1,2))**2+(Old_Old_Tip(1)-new_tip(1,1))**2)
                      c_L2 = sqrt((Old_Tip(2)-new_tip(1,2))**2+(Old_Tip(1)-new_tip(1,1))**2)
                      !:::::::::::::::::::::::::::::
                      ! Temporary variable: stores the previous crack coordinates
                      !:::::::::::::::::::::::::::::
                      tem_Crack_Coor(1:Old_num_Poi,1:2) = Crack_Coor(i_C,1:Old_num_Poi,1:2)
                      !::::::::::::::::::::::::::
                      ! If it is a newly formed arc-shaped crack
                      !::::::::::::::::::::::::::
                      if (Yes_New_Arc_Cr.eqv. .True.)then
                          ! Updated the number of crack coordinate points: increased by 2
                          Each_Cr_Poi_Num(i_C)=Each_Cr_Poi_Num(i_C)+2
                          ! Add the newly created crack tip
                          Crack_Coor(i_C,1,1) = c_Arc_Cr(10)   
                          Crack_Coor(i_C,1,2) = c_Arc_Cr(11) 
                          Crack_Coor(i_C,2,1) = c_Arc_Cr(8) 
                          Crack_Coor(i_C,2,2) = c_Arc_Cr(9) 
                          ! Move the coordinates of the rear crack 2 positions backward
                          Crack_Coor(i_C,3:Old_num_Poi+2,1:2) = tem_Crack_Coor(1:Old_num_Poi,1:2)
                          ! Mark the crack tip as an arc-shaped crack tip (crack tip type 7)
                          Crack_Tip_Type(i_C,i_Tip)=7
                          ! Add the data for the corresponding arc
                          Arc_Crack_Coor(i_C,1,1:11) =  c_Arc_Cr(1:11)   
                      !:::::::::::::::::::::::::::::::::::::::    
                      ! If updating the coordinates of the arc crack points (crack tips)
                      !:::::::::::::::::::::::::::::::::::::::   
                      elseif(Yes_Update_Arc_Cr.eqv. .True.) then
                          ! The number of crack coordinate points remains unchanged
                          !nothing
                          ! Add the newly created crack tip
                          Crack_Coor(i_C,1,1) = c_Arc_Cr(10)   
                          Crack_Coor(i_C,1,2) = c_Arc_Cr(11)
                          ! The crack tip continues to be marked as an arc-shaped crack tip (crack tip type 7)
                          Crack_Tip_Type(i_C,i_Tip)=7
                          ! Update the data for the corresponding arc
                          Arc_Crack_Coor(i_C,1,1:11)=c_Arc_Cr(1:11)                    
                      !:::::::::::::::::::   
                      ! If it is a standard split tip
                      !::::::::::::::::::: 
                      else
                          !Check if the old and new tips are in a common element.
                          call Cal_Ele_Num_by_Coors(Old_Tip(1),Old_Tip(2),Old_OUT_Elem)
                          if(Old_OUT_Elem == New_OUT_Elem)then
                              Yes_In_OneEl = .True.
                          endif
                          !if the old and new tips are in a common element,then:
                          if(Yes_In_OneEl .eqv..True.)then
                              !Modify the original tip location.
                              Crack_Coor(i_C,1,1:2) = new_tip(i_Tip,1:2)   
                          !if the old and new tips are not in a common element,but c_L is not large enough then:
                          elseif((Yes_In_OneEl .eqv..False.) .and. (c_L < 0.99*Ori_delta_L_Growth))then
                              !Modify the original tip location.
                              Crack_Coor(i_C,1,1:2) = new_tip(i_Tip,1:2)   
                          !if the old and new tips are not in a common element,and c_L is large enough then:
                          elseif((Yes_In_OneEl .eqv..False.) .and. (c_L >=0.99*Ori_delta_L_Growth))then
                             ! Update the number of crack coordinate points
                             Each_Cr_Poi_Num(i_C)=Each_Cr_Poi_Num(i_C)+1
                             ! Add the newly formed crack tip
                             Crack_Coor(i_C,1,1:2) = new_tip(i_Tip,1:2)   
                             ! Move the coordinates of the crack at the back one position backward
                             Crack_Coor(i_C,2:Old_num_Poi+1,1:2) = tem_Crack_Coor(1:Old_num_Poi,1:2)
                          endif
                      endif
                      
                      
                    !IIIIIIIIIIIIIIIIIIIIIIIIIIIIII
                    ! If it is tip number 2
                    !IIIIIIIIIIIIIIIIIIIIIIIIIIIIII
                    elseif(i_Tip==2) then
                      Old_num_Poi = Each_Cr_Poi_Num(i_C)
                      Old_Tip(1:2) = Crack_Coor(i_C,Old_num_Poi,1:2)
                      Old_Old_Tip(1:2) =  Crack_Coor(i_C,Old_num_Poi-1,1:2)
                      c_L = sqrt((Old_Old_Tip(2)-new_tip(2,2))**2+(Old_Old_Tip(1)-new_tip(2,1))**2)
                      c_L2 = sqrt((Old_Tip(2)-new_tip(2,2))**2+(Old_Tip(1)-new_tip(2,1))**2)

                      !::::::::::::::::::::::::::
                      ! If it is a newly formed arc-shaped crack
                      !::::::::::::::::::::::::::
                      if (Yes_New_Arc_Cr .eqv. .True.)then
                          ! Updated the number of crack coordinate points: 2 added
                          Each_Cr_Poi_Num(i_C)=Each_Cr_Poi_Num(i_C)+2
                          ! Add the newly created crack tip
                          Crack_Coor(i_C,Old_num_Poi+1,1) = c_Arc_Cr(8)   
                          Crack_Coor(i_C,Old_num_Poi+1,2) = c_Arc_Cr(9) 
                          Crack_Coor(i_C,Old_num_Poi+2,1) = c_Arc_Cr(10) 
                          Crack_Coor(i_C,Old_num_Poi+2,2) = c_Arc_Cr(11) 
                          ! Mark the crack tip as an arc-shaped crack tip (crack tip type 7)
                          Crack_Tip_Type(i_C,i_Tip)=7
                          ! Add the data for the corresponding arc
                          Arc_Crack_Coor(i_C,Each_Cr_Poi_Num(i_C)-1,1:11) = c_Arc_Cr(1:11)   
                      !:::::::::::::::::::::::::::::::::::::::    
                      ! If updating the coordinates of the arc crack points (crack tips)
                      !:::::::::::::::::::::::::::::::::::::::   
                      elseif(Yes_Update_Arc_Cr .eqv. .True.) then
                          ! The number of crack coordinate points remains unchanged
                          !nothing
                          ! Update crack tip
                          Crack_Coor(i_C,Old_num_Poi,1) = c_Arc_Cr(10) 
                          Crack_Coor(i_C,Old_num_Poi,2) = c_Arc_Cr(11) 
                          ! The crack tip continues to be labeled as an arc-shaped crack tip (crack tip type 7)
                          Crack_Tip_Type(i_C,i_Tip)=7
                          ! Update the data for the corresponding arc
                          Arc_Crack_Coor(i_C,Old_num_Poi-1,1:11) = c_Arc_Cr(1:11)       
                      !:::::::::::::::::::   
                      ! If it is a standard split tip
                      !::::::::::::::::::: 
                      else     
                          !Check if the old and new tips are in a common element.
                          call Cal_Ele_Num_by_Coors(Old_Tip(1),Old_Tip(2),Old_OUT_Elem)
                          if(Old_OUT_Elem == New_OUT_Elem)then
                              Yes_In_OneEl = .True.
                          endif
                          !if the old and new tips are in a common element,then:
                          if(Yes_In_OneEl .eqv..True.)then
                              !Modify the original tip location.
                              Crack_Coor(i_C,Old_num_Poi,1:2) = new_tip(i_Tip,1:2)
                          !if the old and new tips are not in a common element,but c_L is not large enough then:
                          elseif((Yes_In_OneEl .eqv..False.) .and.  (c_L < 0.99*Ori_delta_L_Growth))then
                              !Modify the original tip location.
                              Crack_Coor(i_C,Old_num_Poi,1:2) =  new_tip(i_Tip,1:2)
                          !if the old and new tips are not in a common element,then:
                          elseif((Yes_In_OneEl .eqv..False.) .and. (c_L >= 0.99*Ori_delta_L_Growth))then
                              ! Update the number of crack coordinate points
                              Each_Cr_Poi_Num(i_C)=Each_Cr_Poi_Num(i_C)+1
                              ! Add the newly created crack tip
                              Crack_Coor(i_C,Old_num_Poi+1,1:2) = new_tip(i_Tip,1:2)
                          endif
                      endif
                      
                    end if
                  endif
                  ! Information Output (Different Scenarios)
                  if (Yes_New_Arc_Cr.eqv. .True.)then
                      write(*,1601) i_Tip,i_C    
                      write(*,1602)     
                  elseif(Yes_Update_Arc_Cr.eqv. .True.) then
                      write(*,1701) i_Tip,i_C                    
                  else
                      write(*,1001) i_Tip,i_C
                  endif
              end if   
              
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! Junction Detection Related
              !      ---------------------------------------------------------
              ! Check the intersection of cracks to see if there is a junction point (divided into HF cases with
              ! naturally bidirectionally extended bonded cracks, case 1,
              ! HF with natural cemented fractures case 2, and friction natural fractures case 3, and ordinary
              ! situation case 4)
              ! In both the first and third cases, it is necessary to consider the problem of the active crack tip
              ! being too close to the passive crack.
              ! The second case generates an L-shaped crack, without a junction intersection.
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
              ! case1: Hydraulic fracturing involving natural fractures, specifically bilateral propagation
              ! cemented natural fractures
              ! Current functional limitations: (1) Assume that HF only intersects with NF, and HF itself does not
              ! intersect with HF.
              ! (2) A natural fracture can have only one intersection point
              !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
              if((Key_Analysis_Type==3 .or.Key_Analysis_Type==4)  &
&                .and. num_Na_Crack>=1 .and. Key_Na_Crack_Type==1)then
                  !The new segment.
                  New_Seg_A(1:2)= new_tip(i_Tip,1:2)
                  New_Seg_B(1:2)= cr_Tip(i_Tip,1:2) 
              !Loop through each natural crack
              do j_NC =1,num_Na_Crack    
                  !/////////////////////////////////////////////////////////////////////////////////
                  ! If the current natural fracture has already been parasitized by i_C, no further
                  ! detection will be performed, and the next loop variable will be executed.
                  !/////////////////////////////////////////////////////////////////////////////////
                  if(j_NC ==Cracks_NF_JiS_Cr_Num(i_C)) then
                      goto 789
                  endif
                  ! Number of coordinate points of j_C cracks
                  nPt_Na = Each_Na_Cr_Poi_Num(j_NC)
                  do jPt = 2,nPt_Na
                    !The current segment.
                    c_Seg_C = Na_Crack_Coor(j_NC,jPt-1,1:2)
                    c_Seg_D = Na_Crack_Coor(j_NC,jPt  ,1:2) 
                    !Check if the new segment intersects with the current segment.
                    call Tool_Intersection(New_Seg_A,New_Seg_B,c_Seg_C,c_Seg_D,Intersec_x,Intersec_y,Yes_Cross)
                    !If intersect, then:
                    if(Yes_Cross .eqv. .True.) then
                      ! A crucial correction process to ensure that each intersection only requires one Junction
                      ! enhancement element.
                      if(Key_Junction_Check==1)then
                          call Cal_and_Check_Junction_Point &
                         (ifra,iter,New_Seg_B,c_Seg_C,c_Seg_D,Intersec_x,Intersec_y)
                      endif
                      Yes_HF_met_NF = .True.
                      tem_Intersec = [Intersec_x,Intersec_y]
                      tem_Positive_num = i_C
                      tem_NF_num = j_NC
                      ! Added detection: if the distance from the intersection point to the tip of a natural crack is less
                      ! than twice the length of the reinforcement unit, the intersection point is ignored. 2018-10-31
                      ! Explanation: The reason for adding this check is that if the intersection point is very close to
                      ! the tip of a natural crack, it often triggers the following error.
                      !Error :: Junction enrichment and tip enrichment cannot coexist for a node!
                      c_dis_1 = Tool_Function_2Point_Dis(tem_Intersec,c_Seg_C) 
                      c_dis_2 = Tool_Function_2Point_Dis(tem_Intersec,c_Seg_D)      
                      if(c_dis_1<= 2.0D0*Ave_Elem_L_Enrich .or. c_dis_2<= 2.0D0*Ave_Elem_L_Enrich) then
                         Yes_HF_met_NF = .False. 
                         goto 953
                      endif                          
                      if (i_Tip==1) then
                        ! Correct crack tip coordinates
                        Crack_Coor(i_C,1,1:2) = [Intersec_x,Intersec_y]  
                        write(*,2013) i_Tip,i_C,j_NC
                        write(*,2014) i_C
                      elseif(i_Tip==2) then
                        ! Correct crack tip coordinates
                        num_pt = Each_Cr_Poi_Num(i_C)
                        Crack_Coor(i_C,num_pt,1:2) = [Intersec_x,Intersec_y]
                        write(*,2013) i_Tip,i_C,j_NC
                        write(*,2014) i_C
                      end if
                      953 continue
                    !If not intersect, then calculate the distance from 
                    !the new tip to the currentsegment.                                                                          
                    elseif (Yes_Cross .eqv..False.) then 
                      c_Seg_CD(1,1:2) = c_Seg_C
                      c_Seg_CD(2,1:2) = c_Seg_D
                      call Cal_Signed_Distance(c_Seg_CD,new_tip(i_Tip,1:2),Signed_Dis)
                      !If the new tip is too close to the current segment, then try a new tip by 2 times △l.
                      if (abs(Signed_Dis) <= Check_Distance) then
                          tried_New_Tip(i_Tip,1) =  cr_Tip(i_Tip,1)  + TWO*delta_L_Growth*cos(theta(i_Tip))
                          tried_New_Tip(i_Tip,2) =  cr_Tip(i_Tip,2)  + TWO*delta_L_Growth*sin(theta(i_Tip)) 
                          !Recheck if the tried new tip segment intersects with the current segment.
                          Tri_New_Seg_A(1:2)=tried_New_Tip(i_Tip,1:2)
                          Tri_New_Seg_B(1:2)=cr_Tip(i_Tip,1:2)  
                          call Tool_Intersection(Tri_New_Seg_A,Tri_New_Seg_B,c_Seg_C,c_Seg_D, &
                                   Tri_Intersec_x,Tri_Intersec_y,Tri_Yes_Cross)  
                          !if intersect, then: 
                          if(Tri_Yes_Cross.eqv..True.)then  
                             ! A crucial correction process to ensure that each intersection only requires one Junction
                             ! enhancement element.
                             if(Key_Junction_Check==1)then
                                 call Cal_and_Check_Junction_Point(ifra,iter,Tri_New_Seg_B, &
                                 c_Seg_C,c_Seg_D,Tri_Intersec_x,Tri_Intersec_y)
                             endif
                             Yes_HF_met_NF = .True.
                             tem_Intersec = [Tri_Intersec_x,Tri_Intersec_y]
                             tem_NF_num = j_NC
                             tem_Positive_num = i_C
                             ! Added detection: if the distance from the intersection point to the tip of a natural crack is less
                             ! than twice the length of the reinforcement unit, the intersection point is ignored. 2018-10-31
                             ! Explanation: The reason for adding this check is that if the intersection point is very close to
                             ! the tip of a natural crack, the following error often appears.
                             !Error :: Junction enrichment and tip enrichment cannot coexist for a node!
                             c_dis_1 = Tool_Function_2Point_Dis(tem_Intersec,c_Seg_C) 
                             c_dis_2 = Tool_Function_2Point_Dis(tem_Intersec,c_Seg_D)      
                             if(c_dis_1<= 2.0D0*Ave_Elem_L_Enrich .or. c_dis_2<= 2.0D0*Ave_Elem_L_Enrich)then
                                 Yes_HF_met_NF = .False. 
                                 goto 954
                             endif
                             
                             if (i_Tip==1) then
                               ! Correct crack tip coordinates
                              Crack_Coor(i_C,1,1:2)=[Tri_Intersec_x,Tri_Intersec_y]
                               write(*,2013) i_Tip,i_C,j_NC
                               write(*,2014) i_C
                             elseif(i_Tip==2) then
                               ! Correct crack tip coordinates
                               num_pt = Each_Cr_Poi_Num(i_C)
                               Crack_Coor(i_C,num_pt,1:2)=[Tri_Intersec_x,Tri_Intersec_y]
                               write(*,2013) i_Tip,i_C,j_NC
                               write(*,2014) i_C
                             end if   
                             954 continue                              
                           end if    
                         end if   
                       end if
                  enddo
                  789 continue
                enddo
                ! If a convergence occurs:
                ! This forms new preliminary cracks (New_Crack) and updates the number of cracks.
                ! During the program's iterative calculations, num_crack and the crack coordinates cannot be updated
                ! and can only be changed here.
                ! Calculate first, then update at the end of the program. The length of the new crack is five times
                ! the element length.
                ! Each side is 2.5 times long
                if(Yes_HF_met_NF .eqv. .True.)then    
                    num_New_Crack = num_New_Crack + 1
                    Na_Cr_num(num_New_Crack) = tem_NF_num
                    Positive_Cr_num(num_New_Crack)=tem_Positive_num
                    NF_Tip1 = Na_Crack_Coor(tem_NF_num,1,1:2)
                    NF_Tip2 = Na_Crack_Coor(tem_NF_num,nPt_Na,1:2)
                    ! Tip1 point of the new fracture generated from the HF and NF intersection along the Tip1 direction
                    call Tool_Point_Giv_2P_and_L(      &
                        tem_Intersec,NF_Tip1,HLF*L_New_Crack,&
                        New_Crack(num_New_Crack,1,1:2),&
                        Statu_Tip1,L_C_Tip1)
                    ! Tip2 point of a new fracture generated from the HF and NF intersection along the Tip2 direction
                    call Tool_Point_Giv_2P_and_L(       &
                        tem_Intersec,NF_Tip2,HLF*L_New_Crack,&
                        New_Crack(num_New_Crack,2,1:2),  &
                        Statu_Tip2,L_C_Tip2)
                    ! Determine the status of the newly generated crack (for each crack tip): =0, non-parasitic crack;
                    ! =-1, incomplete parasitic; =1, parasitic completed
                    if(Statu_Tip1==1 .or. Statu_Tip1==0)then
                        New_Crack_Statu(num_New_Crack,1) = 1
                    elseif(Statu_Tip1==-1)then
                        New_Crack_Statu(num_New_Crack,1) = -1
                    endif
                    ! Determine the status of the newly generated cracks (for each crack tip): =0, non-parasitic crack;
                    ! =-1, incomplete parasitic; =1, parasitic completed
                    if(Statu_Tip2==1 .or. Statu_Tip2==0)then
                        New_Crack_Statu(num_New_Crack,2) = 1
                    elseif(Statu_Tip2==-1)then
                        New_Crack_Statu(num_New_Crack,2) = -1
                    endif
                endif
                !Cracks_NF_JiS_Stat
              !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
              ! case2: In the situation of hydraulic fracturing containing natural fractures, if the natural
              ! fractures are cemented, then no junction occurs.
              ! And directly generate L-shaped cracks
              !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
              elseif((Key_Analysis_Type==3 .or.Key_Analysis_Type==4) &
                      .and. num_Na_Crack >= 1.and. Key_Na_Crack_Type ==2)then
                !The new segment.
                New_Seg_A(1:2)= new_tip(i_Tip,1:2)
                New_Seg_B(1:2)= cr_Tip(i_Tip,1:2) 
              !Loop through each natural crack
              do j_NC =1,num_Na_Crack    
                  !/////////////////////////////////////////////////////////////////////////////////
                  ! If the current natural fracture has already been parasitized by i_C, no further
                  ! detection will be performed, and the next loop variable will be executed.
                  !/////////////////////////////////////////////////////////////////////////////////
                  if(j_NC ==Cracks_NF_JiS_Cr_Num(i_C)) then
                      goto 790
                  endif
                  ! Number of coordinate points of j_C cracks
                  nPt_Na = Each_Na_Cr_Poi_Num(j_NC)
                  do jPt = 2,nPt_Na
                    !The current segment.
                    c_Seg_C = Na_Crack_Coor(j_NC,jPt-1,1:2)
                    c_Seg_D = Na_Crack_Coor(j_NC,jPt  ,1:2) 
                    !Check if the new segment intersects with the current segment.
                    call Tool_Intersection(New_Seg_A,New_Seg_B,c_Seg_C,c_Seg_D,Intersec_x,Intersec_y,Yes_Cross)
                    !If intersect, then:
                    if(Yes_Cross .eqv. .True.) then
                      Yes_HF_met_NF = .True.
                      tem_Intersec = [Intersec_x,Intersec_y]
                      tem_Positive_num = i_C
                      tem_NF_num = j_NC
                      if (i_Tip==1) then
                        ! Correct crack tip coordinates
                        Crack_Coor(i_C,1,1:2) = [Intersec_x,Intersec_y]  
                        write(*,2013) i_Tip,i_C,j_NC
                        write(*,2014) i_C
                      elseif(i_Tip==2) then
                        ! Correct crack tip coordinates
                        num_pt = Each_Cr_Poi_Num(i_C)
                        Crack_Coor(i_C,num_pt,1:2) = [Intersec_x,Intersec_y]
                        write(*,2013) i_Tip,i_C,j_NC
                        write(*,2014) i_C
                      end if
                    !If not intersect, then calculate the distance 
                    !from the new tip to the current segment.                                                                      &
                    elseif (Yes_Cross .eqv..False.) then 
                      c_Seg_CD(1,1:2) = c_Seg_C
                      c_Seg_CD(2,1:2) = c_Seg_D
                      call Cal_Signed_Distance(c_Seg_CD,new_tip(i_Tip,1:2),Signed_Dis)
                      !If the new tip is too close to the current segment, then try a new tip by 2 times △l.
                      if (abs(Signed_Dis) <= Check_Distance) then
                          tried_New_Tip(i_Tip,1) =  cr_Tip(i_Tip,1)  + TWO*delta_L_Growth*cos(theta(i_Tip))
                          tried_New_Tip(i_Tip,2) =  cr_Tip(i_Tip,2)  + TWO*delta_L_Growth*sin(theta(i_Tip)) 
                          !Recheck if the tried new tip segment intersects with the current segment.
                          Tri_New_Seg_A(1:2)=tried_New_Tip(i_Tip,1:2)
                          Tri_New_Seg_B(1:2)=cr_Tip(i_Tip,1:2)  
                          call Tool_Intersection(Tri_New_Seg_A,Tri_New_Seg_B, &
                                   c_Seg_C,c_Seg_D,Tri_Intersec_x,Tri_Intersec_y, &
                                   Tri_Yes_Cross)  
                          !if intersect, then: 
                          if(Tri_Yes_Cross.eqv..True.)then  
                             Yes_HF_met_NF = .True.
                             tem_Intersec = [Tri_Intersec_x,Tri_Intersec_y]
                             tem_NF_num = j_NC
                             tem_Positive_num = i_C
                             if (i_Tip==1) then
                               ! Correct crack tip coordinates
                              Crack_Coor(i_C,1,1:2)=[Tri_Intersec_x,Tri_Intersec_y]  
                               write(*,2013) i_Tip,i_C,j_NC
                               write(*,2014) i_C
                             elseif(i_Tip==2) then
                               ! Correct crack tip coordinates
                               num_pt = Each_Cr_Poi_Num(i_C)
                               Crack_Coor(i_C,num_pt,1:2)=[Tri_Intersec_x,Tri_Intersec_y]
                               write(*,2013) i_Tip,i_C,j_NC
                               write(*,2014) i_C
                             end if                                
                           end if    
                         end if   
                       end if
                  enddo
                  790 continue
                enddo
                ! If an intersection occurs, a new preliminary crack (New_Crack) will be formed, and the number of
                ! cracks will be updated.
                ! During the program's iterative calculations, num_crack and the crack coordinates cannot be updated
                ! and can only be changed here.
                ! Calculate first, and then update at the end of the program. The length of the new crack is five
                ! times the element length.
                ! Each side is 2.5 times long
                if(Yes_HF_met_NF .eqv. .True.)then    
                    NF_Tip1 = Na_Crack_Coor(tem_NF_num,1,1:2)
                    NF_Tip2 = Na_Crack_Coor(tem_NF_num,nPt_Na,1:2)
                    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                    ! Determine c_Na_Tip=1 or 2, that is, decide whether the extension is along
                    ! crack tip 1 or crack tip 2 of the bonded crack.
                    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                    ! Option 1, determine based on the size of the angle; see my notes V5 for details.
                    Vector_J_Tip1 = NF_Tip1 - tem_Intersec
                    Vector_J_Tip2 = NF_Tip2 - tem_Intersec
                    Point_A = Crack_Coor(i_C,Each_Cr_Poi_Num(i_C)-1,1:2)
                    Vector_J_A = Point_A - tem_Intersec 
                    dot_J_Tip1=dot_product(Vector_J_A,Vector_J_Tip1)
                    dot_J_Tip2=dot_product(Vector_J_A,Vector_J_Tip2)
                    ! Take the angle with a negative dot product (obtuse angle)
                    if(dot_J_Tip1 <= ZR) c_Na_Tip=1
                    if(dot_J_Tip2 <= ZR) c_Na_Tip=2
                    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                    ! If it extends along the tip of crack No. 1 in the cemented joint
                    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                    if(c_Na_Tip==1)then
                        ! Generate a new crack tip of i_C from the intersection of HF and NF along the Tip1 direction
                        call Tool_Point_Giv_2P_and_L(     &
                            tem_Intersec,NF_Tip1,HLF*L_New_Crack,New_Na_Tip,Statu_Tip1,L_C_Tip1)
                        ! Determine the status of the newly generated crack (for each crack tip): =0, non-parasitic crack;
                        ! =-1, incomplete parasitic; =1, parasitic completed
                        if(Statu_Tip1==1 .or. Statu_Tip1==0)then
                            Cracks_NF_JiS_Stat(i_C,i_Tip) = 1
                        elseif(Statu_Tip1==-1)then
                            Cracks_NF_JiS_Stat(i_C,i_Tip) = -1
                        endif
                    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                    ! If it extends along the No. 2 tip of the bonded crack
                    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                    elseif(c_Na_Tip==2)then
                        ! Generate a new crack tip of i_C from the intersection of HF and NF along the Tip2 direction
                        call Tool_Point_Giv_2P_and_L(  &
                            tem_Intersec,NF_Tip2,HLF*L_New_Crack,New_Na_Tip,Statu_Tip2,L_C_Tip2)
                        ! Determine the status of the newly generated crack (for each crack tip): =0, non-parasitic crack;
                        ! =-1, incomplete parasitic; =1, parasitic completed
                        if(Statu_Tip2==1 .or. Statu_Tip2==0)then
                            Cracks_NF_JiS_Stat(i_C,i_Tip) = 1
                        elseif(Statu_Tip2==-1)then
                            Cracks_NF_JiS_Stat(i_C,i_Tip) = -1
                        endif
                    endif
                    ! After crack i_C intersects with the cemented crack, does it propagate along the tip 1 of the
                    ! natural crack or along the tip 2?
                    Cracks_NF_Cement_Tip_Num(i_C,i_Tip) = c_Na_Tip                
                    ! Current natural fracture number corresponding to the crack
                    Cracks_NF_JiS_Cr_Num(i_C) =tem_NF_num
                    ! Add new crack coordinate point
                    Each_Cr_Poi_Num(i_C) = Each_Cr_Poi_Num(i_C) + 1
                    Crack_Coor(i_C,Each_Cr_Poi_Num(i_C),1:2) = New_Na_Tip
                endif
              !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
              ! case 3: Hydraulic fracturing involving natural fractures, specifically friction-type natural
              ! fractures
              ! Current feature limitations: None
              !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
              elseif((Key_Analysis_Type==3 .or.Key_Analysis_Type==4) &
               .and. num_Na_Crack>=1 .and. Key_Na_Crack_Type==3) then
                !The new segment.
                New_Seg_A(1:2)= new_tip(i_Tip,1:2)
                New_Seg_B(1:2)= cr_Tip(i_Tip,1:2) 
                ! Loop through each crack (except i_C)
                do j_C =1,num_Crack   
                 if(j_C .ne. i_C) then
                  c_Na_Crack  = Cracks_fric_NF_num(j_C)
                  ! Number of coordinate points of j_C cracks
                  nPt = Each_Cr_Poi_Num(j_C)  
                  do jPt = 2,nPt
                    !The current segment.
                    c_Seg_C = Crack_Coor(j_C,jPt-1,1:2)
                    c_Seg_D = Crack_Coor(j_C,jPt  ,1:2) 
                    !Check if the new segment intersects with the current segment.
                    call Tool_Intersection(New_Seg_A,New_Seg_B, &
                                c_Seg_C,c_Seg_D,Intersec_x,Intersec_y,Yes_Cross)
                    !If intersect, then:
                    if(Yes_Cross .eqv. .True.) then
                      ! A crucial correction process to ensure that each intersection only requires one Junction
                      ! enhancement element.
                      if(Key_Junction_Check==1)then
                          call Cal_and_Check_Junction_Point &
                         (ifra,iter,New_Seg_B,c_Seg_C,c_Seg_D,Intersec_x,Intersec_y)
                      endif
                      ! If the current passive fracture happens to be a naturally occurring frictional crack, then
                      if(c_Na_Crack /= 0)then
                          Yes_HF_met_NF = .True.
                      endif
                      tem_Intersec = [Intersec_x,Intersec_y]
                      tem_Positive_num = i_C
                      Real_Cr_NF = j_C
                      if (i_Tip==1) then
                        ! Correct crack tip coordinates
                        Crack_Coor(i_C,1,1:2) = [Intersec_x,Intersec_y]  
                        write(*,2013) i_Tip,i_C,j_C
                        write(*,2014) i_C
                      elseif(i_Tip==2) then
                        ! Correct crack tip coordinates
                        num_pt = Each_Cr_Poi_Num(i_C)
                        Crack_Coor(i_C,num_pt,1:2) = [Intersec_x,Intersec_y]
                        write(*,2013) i_Tip,i_C,j_C
                        write(*,2014) i_C
                      end if
                    !If not intersect, then calculate the distance 
                    !from the new tip to the current segment. 
                    elseif (Yes_Cross .eqv..False.) then 
                      c_Seg_CD(1,1:2) = c_Seg_C
                      c_Seg_CD(2,1:2) = c_Seg_D
                      call Cal_Signed_Distance(c_Seg_CD,new_tip(i_Tip,1:2),Signed_Dis)
                      !If the new tip is too close to the current segment, then try a new tip by 2 times △l.
                      if (abs(Signed_Dis) <= Check_Distance) then
                          tried_New_Tip(i_Tip,1) =  cr_Tip(i_Tip,1)  + TWO*delta_L_Growth*cos(theta(i_Tip))
                          tried_New_Tip(i_Tip,2) =  cr_Tip(i_Tip,2)  + TWO*delta_L_Growth*sin(theta(i_Tip)) 
                          !Recheck if the tried new tip segment intersects with the current segment.
                          Tri_New_Seg_A(1:2)=tried_New_Tip(i_Tip,1:2)
                          Tri_New_Seg_B(1:2)=cr_Tip(i_Tip,1:2)  
                          call Tool_Intersection(Tri_New_Seg_A,Tri_New_Seg_B, &
                                   c_Seg_C,c_Seg_D,Tri_Intersec_x,Tri_Intersec_y,Tri_Yes_Cross)  
                          !if intersect, then: 
                          if(Tri_Yes_Cross.eqv..True.)then  
                             ! If the current passive crack happens to be a naturally occurring frictional crack, then
                             if(c_Na_Crack /= 0)then
                                 Yes_HF_met_NF = .True.
                             endif
                             ! A crucial correction process to ensure that each intersection only requires one Junction
                             ! enhancement element.
                             if(Key_Junction_Check==1)then
                                 call Cal_and_Check_Junction_Point(ifra,iter,Tri_New_Seg_B,&
                                 c_Seg_C,c_Seg_D,Tri_Intersec_x,Tri_Intersec_y)
                             endif
                             Yes_HF_met_NF = .True.
                             tem_Intersec = [Tri_Intersec_x,Tri_Intersec_y]
                             Real_Cr_NF = j_C
                             tem_Positive_num = i_C
                             if (i_Tip==1) then
                               ! Correct crack tip coordinates
                              Crack_Coor(i_C,1,1:2)=[Tri_Intersec_x,Tri_Intersec_y]
                               write(*,2013) i_Tip,i_C,j_C
                               write(*,2014) i_C
                             elseif(i_Tip==2) then
                               ! Correct crack tip coordinates
                               num_pt = Each_Cr_Poi_Num(i_C)
                               Crack_Coor(i_C,num_pt,1:2)=[Tri_Intersec_x,Tri_Intersec_y]
                               write(*,2013) i_Tip,i_C,j_C
                               write(*,2014) i_C
                             end if                                
                           end if    
                         end if   
                    end if
                  enddo
                 endif
                enddo
                !''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
                ! If an intersection with a natural fracture occurs, and the natural fracture has not been
                ! previously parasitized:
                ! Then change the length of the natural fractures, change their type to hydraulic fractures, and
                ! allow the fractures
                ! Extension; the length of the new crack is 5 times the element length, with each side being 2.5
                ! times -- 2016-06-09
                !''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
                
                if (Yes_HF_met_NF .eqv. .True.) then
                  if (Cracks_QinS_Stat(Real_Cr_NF)/= 1) then    
                    c_Na_Cr_Num = Cracks_fric_NF_num(Real_Cr_NF)
                    NF_Tip1 = Na_Crack_Coor(c_Na_Cr_Num,1,1:2)
                    NF_Tip2 = Na_Crack_Coor(c_Na_Cr_Num,2,1:2)
                    ! The natural fracture number corresponding to the current parasitic fracture
                    Cracks_NF_JiS_Cr_Num(Real_Cr_NF) = c_Na_Cr_Num
                    ! Tip1 point of the new fracture generated from the HF and NF intersection along the Tip1 direction
                    call Tool_Point_Giv_2P_and_L(        &
                        tem_Intersec,NF_Tip1,HLF*L_New_Crack,&
                        New_Crack_Tip_1(1:2),Statu_Tip1,L_C_Tip1)
                    ! Tip2 point of the new fracture generated from the HF and NF intersection along the Tip2 direction
                    call Tool_Point_Giv_2P_and_L(    &
                        tem_Intersec,NF_Tip2,HLF*L_New_Crack,&
                        New_Crack_Tip_2(1:2),Statu_Tip2,L_C_Tip2)
                    ! Adjust the length of the actual fractures corresponding to the frictional natural fractures
                    ! (participating in the calculation)
                    Crack_Coor(Real_Cr_NF,1,1:2) = New_Crack_Tip_1
                    Crack_Coor(Real_Cr_NF,2,1:2) = New_Crack_Tip_2
                    ! Updated to hydraulic fracturing
                    Cracks_HF_State(Real_Cr_NF) = 1
                    ! Mark as extendable crack
                    Cracks_Allow_Propa(Real_Cr_NF) = 1
                    ! Determine the status of newly generated cracks (each crack tip has): =0, non-parasitic crack; =-1,
                    ! parasitic incomplete; =1, parasitic complete
                    if(Statu_Tip1==1 .or. Statu_Tip1==0)then
                        Cracks_NF_JiS_Stat(Real_Cr_NF,1)  = 1
                    elseif(Statu_Tip1==-1)then
                        Cracks_NF_JiS_Stat(Real_Cr_NF,1)  = -1
                    endif
                    ! Determine the status of the newly generated cracks (for each crack tip): =0, non-parasitic crack;
                    ! =-1, incomplete parasitic; =1, parasitic completed
                    if(Statu_Tip2==1 .or. Statu_Tip2==0)then
                        Cracks_NF_JiS_Stat(Real_Cr_NF,2)  = 1
                    elseif(Statu_Tip2==-1)then
                        Cracks_NF_JiS_Stat(Real_Cr_NF,2)  = -1
                    endif
                    ! Mark the fracture as one that has been eroded by hydraulic fracturing
                    Cracks_QinS_Stat(Real_Cr_NF)  = 1
                  endif
                endif
              !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
              ! case4: General analysis (situations with natural fractures without hydraulic fracturing)
              !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
              else 
                !The new segment.
                New_Seg_A(1:2)=[new_tip(i_Tip,1),new_tip(i_Tip,2)]
                New_Seg_B(1:2)= cr_Tip(i_Tip,1:2) 
              !Loop through each crack except i_C
              do j_C =1,num_Crack    
                 if(j_C .ne. i_C) then
                   ! Number of coordinate points of j_C cracks
                   nPt = Each_Cr_Poi_Num(j_C)  
                     do jPt = 2,nPt
                       !The current segment.
                       c_Seg_C = Crack_Coor(j_C,jPt-1,1:2)
                       c_Seg_D = Crack_Coor(j_C,jPt,1:2)  
                       !Check if the new segment intersects with the current segment.
                       call Tool_Intersection(New_Seg_A,New_Seg_B,c_Seg_C,c_Seg_D,Intersec_x,Intersec_y,Yes_Cross)
                       !If intersect, then:
                       if(Yes_Cross .eqv. .True.) then
                          ! A crucial correction process to ensure that each intersection only requires one Junction
                          ! enhancement element.
                          if(Key_Junction_Check==1)then
                              call Cal_and_Check_Junction_Point(ifra,iter,New_Seg_B,c_Seg_C,c_Seg_D,Intersec_x,Intersec_y)
                          endif
                         if (i_Tip==1) then
                           ! Correct crack tip coordinates
                           Crack_Coor(i_C,1,1:2) = [Intersec_x,Intersec_y]  
                           write(*,1013) i_Tip,i_C,j_C
                           write(*,1014) i_C
                         elseif(i_Tip==2) then
                           ! Correct crack tip coordinates
                           num_pt = Each_Cr_Poi_Num(i_C)
                           Crack_Coor(i_C,num_pt,1:2) = [Intersec_x,Intersec_y]
                           write(*,1013) i_Tip,i_C,j_C
                           write(*,1014) i_C
                         end if
                       !If not intersect, then calculate the distance
                       !from the new tip to the current segment. 
                       elseif (Yes_Cross .eqv..False.) then 
                         c_Seg_CD(1,1:2) = c_Seg_C
                         c_Seg_CD(2,1:2) = c_Seg_D
                         call Cal_Signed_Distance(c_Seg_CD,new_tip(i_Tip,:),Signed_Dis)
                         !If the new tip is too close to the current segment, then try a new tip by 2 times △l.
                         if (abs(Signed_Dis) <= Check_Distance) then
                           ! Multiple attempts were made, each time increasing the crack length until it intersected with the
                           ! crack (then abandoned).
                           !do iii_try =1,20
                            !tried_New_Tip(i_Tip,1) =cr_Tip(i_Tip,1)  + dble(iii_try+1)*delta_L_Growth*cos(theta(i_Tip))
                            !tried_New_Tip(i_Tip,2)= cr_Tip(i_Tip,2)  + dble(iii_try+1)*delta_L_Growth*sin(theta(i_Tip)) 
                             tried_New_Tip(i_Tip,1) =cr_Tip(i_Tip,1)  + TWO*delta_L_Growth*cos(theta(i_Tip))
                             tried_New_Tip(i_Tip,2)= cr_Tip(i_Tip,2)  + TWO*delta_L_Growth*sin(theta(i_Tip)) 
                             !Recheck if the tried new tip segment intersects with the current segment.
                             Tri_New_Seg_A(1:2)=[tried_New_Tip(i_Tip,1),tried_New_Tip(i_Tip,2)]
                             Tri_New_Seg_B(1:2)=cr_Tip(i_Tip,1:2)  
                             call Tool_Intersection(Tri_New_Seg_A,Tri_New_Seg_B, &
                                   c_Seg_C,c_Seg_D,Tri_Intersec_x,Tri_Intersec_y,Tri_Yes_Cross)     
                             !if intersect, then: 
                             if(Tri_Yes_Cross)then
                               ! A crucial correction process to ensure that each intersection only requires one Junction
                               ! enhancement element.
                               if(Key_Junction_Check==1)then
                                 call Cal_and_Check_Junction_Point(ifra,iter,Tri_New_Seg_B, &
                                 c_Seg_C,c_Seg_D,Tri_Intersec_x,Tri_Intersec_y)
                               endif
                               if (i_Tip==1) then
                               ! Correct crack tip coordinates
                                Crack_Coor(i_C,1,1:2)=[Tri_Intersec_x,Tri_Intersec_y]  
                                 write(*,1013) i_Tip,i_C,j_C
                                 write(*,1014) i_C
                               elseif(i_Tip==2) then
                                 ! Correct crack tip coordinates
                                 num_pt = Each_Cr_Poi_Num(i_C)
                                 Crack_Coor(i_C,num_pt,1:2)=[Tri_Intersec_x,Tri_Intersec_y]
                                 write(*,1013) i_Tip,i_C,j_C
                                 write(*,1014) i_C
                               end if
                               exit
                             end if
                           !enddo
                         end if
                       end if
                     end do
                 end if                      
               end do   
              end if
              
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! If there is a circular hole, check whether it is inside the hole or very close to
              ! it.
              !Added on 2017-05-20
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              if(num_Circ_Hole >=1 )then
                  !The new segment.
                  New_Seg_A(1:2)=[new_tip(i_Tip,1),new_tip(i_Tip,2)]
                  New_Seg_B(1:2)= cr_Tip(i_Tip,1:2) 
                  ! Hole Loop
                  loop_i_Hole:do i_Hole =1,num_Circ_Hole
                      x0 =  Hole_Coor(i_Hole,1)
                      y0 =  Hole_Coor(i_Hole,2)
                      R0 =  Hole_Coor(i_Hole,3) 
                      ! Calculate the intersection between the current segment and the current hole
                      call Tool_Intersection_Line_and_Circle(x0,y0,R0,New_Seg_A,New_Seg_B,c_num_Inter,c_State,c_Inter)
                      ! If there is an intersection point, adjust the crack tip to the intersection point (if there are
                      ! two intersection points, adjust it to the first intersection point).
                      if(c_num_Inter >=1) then
                         if (i_Tip==1) then
                           ! Correct crack tip coordinates
                           Crack_Coor(i_C,1,1:2)=c_Inter(1,1:2)
                           Crack_Tip_Type(i_C,i_Tip) = 1 
                           write(*,1023) i_Tip,i_C,i_Hole
                           write(*,1024) i_Tip,i_C
                         elseif(i_Tip==2) then
                           ! Correct crack tip coordinates
                           num_pt = Each_Cr_Poi_Num(i_C)
                           Crack_Coor(i_C,num_pt,1:2)=c_Inter(1,1:2)
                           Crack_Tip_Type(i_C,i_Tip) = 1 
                           write(*,1023) i_Tip,i_C,i_Hole
                           write(*,1024) i_Tip,i_C
                         end if  
                      endif
                      ! If there is no intersection, check if they are very close.
                      if(c_num_Inter ==0) then
                         ! Distance from the new crack tip to the center
                         c_Dis = Tool_Function_2Point_Dis([x0,y0],New_Seg_A)
                         !If the new tip is too close to the hole, then try a new tip by 2 times △l.
                         if (abs(c_Dis-R0) <= Check_Distance) then
                           ! Try multiple times, increasing the crack length each time until it intersects with the hole.
                           do iii_try =1,20
                             tried_New_Tip(i_Tip,1) =cr_Tip(i_Tip,1)+ dble(iii_try+1)*delta_L_Growth*cos(theta(i_Tip))
                             tried_New_Tip(i_Tip,2) =cr_Tip(i_Tip,2)+ dble(iii_try+1)*delta_L_Growth*sin(theta(i_Tip)) 
                             !Recheck if the tried new tip segment intersects with the current hole.
                             Tri_New_Seg_A(1:2)=[tried_New_Tip(i_Tip,1),tried_New_Tip(i_Tip,2)]
                             Tri_New_Seg_B(1:2)=cr_Tip(i_Tip,1:2)  
                            ! Calculate the intersection between the current segment and the current hole
                            call Tool_Intersection_Line_and_Circle(x0,y0,R0,Tri_New_Seg_A,Tri_New_Seg_B,c_num_Inter,c_State,c_Inter)
                            ! If there is an intersection point, adjust the crack tip to the intersection point (if there are
                            ! two intersection points, adjust it to the first intersection point).
                            if(c_num_Inter >=1) then
                            
                            ! Similar to the junction of cracks, this part of the code may also be needed (requiring
                            ! corresponding modifications), otherwise the program may encounter errors in certain situations.
                            ! To be studied in detail
                            
                            ! A crucial correction process to ensure that each intersection only requires one Junction
                            ! enhancement element.
                            !if(Key_Junction_Check==1)then
                            !    call
                            !    Cal_and_Check_Junction_Point(ifra,iter,Tri_New_Seg_B,c_Seg_C,c_Seg_D,Tri_Intersec_x,Tri_Intersec_y)
                            ! Intersec_x and Intersec_y are input and output variables
                            !endif
                               
                             if (i_Tip==1) then
                               ! Correct crack tip coordinates
                               Crack_Coor(i_C,1,1:2)=c_Inter(1,1:2)
                               Crack_Tip_Type(i_C,i_Tip) = 1 
                               write(*,1023) i_Tip,i_C,i_Hole
                               write(*,1024) i_Tip,i_C
                             elseif(i_Tip==2) then
                               ! Correct crack tip coordinates
                               num_pt = Each_Cr_Poi_Num(i_C)
                               Crack_Coor(i_C,num_pt,1:2)=c_Inter(1,1:2)
                               Crack_Tip_Type(i_C,i_Tip) = 1 
                               write(*,1023) i_Tip,i_C,i_Hole
                               write(*,1024) i_Tip,i_C
                             end if
                             exit
                            endif
                          enddo
                         end if
                      endif
                  enddo loop_i_Hole
              endif
              
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              ! If there is an elliptical hole, determine whether it is inside the hole or very
              ! close to the hole.
              !Added on 2020-08-09
              !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              if(num_Ellip_Hole >=1 )then
                  !The new segment.
                  New_Seg_A(1:2)=[new_tip(i_Tip,1),new_tip(i_Tip,2)]
                  New_Seg_B(1:2)= cr_Tip(i_Tip,1:2) 
                  ! Hole Loop
                  loop_i_Hole2:do i_Hole =1,num_Ellip_Hole
                      x0_Hole =  Ellip_Hole_Coor(i_Hole,1)
                      y0_Hole =  Ellip_Hole_Coor(i_Hole,2)
                      a_Hole  =  Ellip_Hole_Coor(i_Hole,3) 
                      b_Hole  =  Ellip_Hole_Coor(i_Hole,4) 
                      theta_Hole =  Ellip_Hole_Coor(i_Hole,5)
                      ! Calculate the intersection between the current segment and the current hole
                      call Tool_Intersection_Line_and_Oblique_Ellipse( &
                          x0_Hole,y0_Hole,a_Hole,b_Hole,theta_Hole, &
                          New_Seg_A,New_Seg_B,c_num_Inter,c_Inter(1,1:2))
                      
                      ! If there is an intersection point, adjust the crack tip to the intersection point (if there are
                      ! two intersection points, adjust it to the first intersection point).
                      if(c_num_Inter >=1) then
                         if (i_Tip==1) then
                           ! Correct crack tip coordinates
                           Crack_Coor(i_C,1,1:2)=c_Inter(1,1:2)
                           Crack_Tip_Type(i_C,i_Tip) = 1 
                           write(*,1023) i_Tip,i_C,i_Hole
                           write(*,1024) i_Tip,i_C
                         elseif(i_Tip==2) then
                           ! Correct crack tip coordinates
                           num_pt = Each_Cr_Poi_Num(i_C)
                           Crack_Coor(i_C,num_pt,1:2)=c_Inter(1,1:2)
                           Crack_Tip_Type(i_C,i_Tip) = 1 
                           write(*,1023) i_Tip,i_C,i_Hole
                           write(*,1024) i_Tip,i_C
                         end if  
                      endif
                      ! If there is no intersection, check if they are very close.
                      if(c_num_Inter ==0) then
                         ! Calculate the distance from the new crack tip to the ellipse
                         call Tool_Intersection_Line_and_Oblique_Ellipse( &
                          x0_Hole,y0_Hole,a_Hole,b_Hole,theta_Hole, &
                          New_Seg_A,[x0_Hole,y0_Hole],c_num_Inter,c_Inter(1,1:2))
                           if(c_num_Inter==1)then
                               c_dis=Tool_Function_2Point_Dis(New_Seg_A,c_Inter(1,1:2))
                           endif
                           ! Calculate
                         !If the new tip is too close to the hole, then try a new tip by 2 times △l.
                         if (c_num_Inter==1  .and. (c_dis<Check_Distance)) then
                           ! Try multiple times, increasing the length of the crack each time until it intersects with the
                           ! hole.
                           do iii_try =1,20
                             tried_New_Tip(i_Tip,1) =cr_Tip(i_Tip,1)  + dble(iii_try+1)*delta_L_Growth*cos(theta(i_Tip))
                             tried_New_Tip(i_Tip,2) =cr_Tip(i_Tip,2)  + dble(iii_try+1)*delta_L_Growth*sin(theta(i_Tip)) 
                             !Recheck if the tried new tip segment intersects with the current hole.
                             Tri_New_Seg_A(1:2)=[tried_New_Tip(i_Tip,1),tried_New_Tip(i_Tip,2)]
                             Tri_New_Seg_B(1:2)=cr_Tip(i_Tip,1:2)  
                            ! Calculate the intersection points between the current segment and the current elliptical hole
                            call Tool_Intersection_Line_and_Oblique_Ellipse(x0_Hole,y0_Hole,a_Hole,b_Hole,theta_Hole, &
                               Tri_New_Seg_A,Tri_New_Seg_B,c_num_Inter,c_Inter(1,1:2))

                            ! If there is an intersection point, adjust the crack tip to the intersection point (if there are
                            ! two intersection points, adjust it to the first intersection point).
                            if(c_num_Inter >=1) then
                            
                            ! Similar to the junction between cracks, this part of the code may also need (corresponding
                            ! modifications), otherwise the program could fail in certain situations.
                            ! To be studied in detail
                            
                            ! !A crucial correction process to ensure that each intersection only requires one Junction
                            ! enhancement element
                            !if(Key_Junction_Check==1)then
                            !    call
                            !    Cal_and_Check_Junction_Point(ifra,iter,Tri_New_Seg_B,c_Seg_C,c_Seg_D,Tri_Intersec_x,Tri_Intersec_y)
                            ! !Intersec_x, Intersec_y are input and output variables
                            !endif
                               
                             if (i_Tip==1) then
                               ! Correct crack tip coordinates
                               Crack_Coor(i_C,1,1:2)=c_Inter(1,1:2)
                               Crack_Tip_Type(i_C,i_Tip) = 1 
                               write(*,1023) i_Tip,i_C,i_Hole
                               write(*,1024) i_Tip,i_C
                             elseif(i_Tip==2) then
                               ! Correct crack tip coordinates
                               num_pt = Each_Cr_Poi_Num(i_C)
                               Crack_Coor(i_C,num_pt,1:2)=c_Inter(1,1:2)
                               Crack_Tip_Type(i_C,i_Tip) = 1 
                               write(*,1023) i_Tip,i_C,i_Hole
                               write(*,1024) i_Tip,i_C
                             end if
                             exit
                            endif
                          enddo
                         end if
                      endif
                  enddo loop_i_Hole2
              endif
          end if
99     continue      
      end do loop_i_Tip
100 continue   
end do loop_i_C

!****************************************************************************************************
! c For the issue of natural fracture intersections during hydraulic fracturing, update the newly
! formed fracture network after the intersections
! c     This situation only occurs when Key_Na_Crack_Type == 1, because:
! c Key_Na_Crack_Type==1 indicates a bilateral extended bonded crack, and only then will new cracks
! form
! c Key_Na_Crack_Type==2 refers to a single-sided extended cemented crack, with the number of cracks
! always being 1
! c Key_Na_Crack_Type==3 represents real friction cracks, which are involved in friction
! calculations from the initial moment, so no new cracks are generated.
!****************************************************************************************************
if((Key_Analysis_Type==3 .or.Key_Analysis_Type==4) .and. num_Na_Crack>=1)then
      if(num_New_Crack>0)then
          print *,'    Updating new created fractures......' 
          do i_New_C = 1,num_New_Crack
              num_Crack =  num_Crack +1
              ! Coordinates of newly generated cracks
              Crack_Coor(num_Crack,1,1:2) = New_Crack(i_New_C,1,1:2)
              Crack_Coor(num_Crack,2,1:2) = New_Crack(i_New_C,2,1:2)
              ! Number of crack coordinate points
              Each_Cr_Poi_Num(num_Crack)  = 2
              ! The natural fracture number corresponding to the newly calculated fracture generated by the
              ! intersection of HF and NF
              Cracks_NF_JiS_Cr_Num(num_Crack)=Na_Cr_num(i_New_C) 
              ! Parasitic state of newly generated cracks (used to determine the direction of crack propagation)
              Cracks_NF_JiS_Stat(num_Crack,1)=New_Crack_Statu(i_New_C,1)
              Cracks_NF_JiS_Stat(num_Crack,2)=New_Crack_Statu(i_New_C,2)
              ! Add new crack marker as water-driven crack
              Cracks_HF_State(num_Crack) = 1
              ! The newly generated crack is marked as an extensible crack
              Cracks_Allow_Propa(num_Crack) = 1
              ! Mark as newly generated T-shaped cracks as passive cracks, and save the corresponding active crack
              ! numbers (used to determine initial water pressure, etc.)
              Cracks_NF_T_Stat(num_Crack,1)= -1
              Cracks_NF_T_Stat(num_Crack,2)=Positive_Cr_num(i_New_C)
          enddo
      endif
endif

RETURN
end subroutine Check_Crack_Grows