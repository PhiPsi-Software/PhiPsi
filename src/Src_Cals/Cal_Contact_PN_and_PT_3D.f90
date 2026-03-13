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
 
subroutine Cal_Contact_PN_and_PT_3D(iter,ifra,Counter_Iter,i_NR_P, &
                      c_Total_Freedom,c_num_freeDOF,&
                       Kn,Kn_Gauss,Kt1_Gauss,Kt2_Gauss,&
                       fric_mu,c_DISP,delta_u_a,&
                       CT_State_Gauss,c_Elem_Conta_Sta,&
                       PC_Gauss_x,PC_Gauss_y,PC_Gauss_z,Num_Indent)
! (1): Transform the contact forces at each Gauss point of the crack surface elements during the
! contact iteration process into the global rectangular coordinate system
! The nodes of the contact elements and the fluid elements completely coincide, which facilitates
! the programming.
! (2): Calculate the contact status of each Gauss point during each contact iteration step, and save
! it in CT_State_Gauss, along with
! The element contact status used for post-processing is saved in c_Elem_Conta_Sta.
! (3): During the iteration process, once a slip occurs, the status of that element is marked as a
! slipping state (c_Elem_Conta_Sta).


!-----------------------------
! Read public variable module
!-----------------------------
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Crack_Common
use Global_Crack_3D
use Global_HF
use Global_Material
use Global_Field_Problem
use Global_Crack
use omp_lib

implicit none
integer, intent(in)::iter,ifra,Counter_Iter,i_NR_P
integer, intent(in)::c_Total_Freedom,c_num_freeDOF
real(kind=FT),intent(in)::c_DISP(c_Total_Freedom)
real(kind=FT),intent(in)::delta_u_a(c_Total_Freedom)
real(kind=FT),intent(in)::Kn,fric_mu
real(kind=FT),intent(inout)::Kn_Gauss(num_Crack,Max_Max_N_FluEl_3D)
real(kind=FT),intent(inout)::Kt1_Gauss(num_Crack,Max_Max_N_FluEl_3D)
real(kind=FT),intent(inout)::Kt2_Gauss(num_Crack,Max_Max_N_FluEl_3D)
integer,intent(inout)::CT_State_Gauss(num_Crack,Max_Max_N_FluEl_3D)             
integer,intent(inout)::c_Elem_Conta_Sta(Num_Elem,num_Crack)
real(kind=FT),intent(inout):: &
            PC_Gauss_x(num_Crack,Max_Max_N_FluEl_3D),&
            PC_Gauss_y(num_Crack,Max_Max_N_FluEl_3D),&
            PC_Gauss_z(num_Crack,Max_Max_N_FluEl_3D)
integer Num_Indent
integer i_C
integer i_CT_Elem
real(kind=FT) c_Aperture
real(kind=FT) all_Apertures(num_Crack,Max_Max_N_FluEl_3D)
real(kind=FT) c_PN,c_PT1,c_PT2
real(kind=FT) c_Kt1,c_Kt2
integer Old_Elem_Conta_Sta(Num_Elem,num_Crack)
real(kind=FT) c_Kn
real(kind=FT) c_Center(3),ori_n(3),T(3,3),tran_T(3,3)
real(kind=FT) PC_Gauss(3)
integer Solid_El
real(kind=FT) Relative_Disp(3)

2001 FORMAT(12X,'Max aperture is:   ',F14.5,' mm')       
2002 FORMAT(12X,'Min aperture is:   ',F14.5,' mm') 
3001 FORMAT(12X,'Max contact stress xx is:   ',F15.4,' MPa')       
3002 FORMAT(12X,'Min contact stress xx is:   ',F15.4,' MPa')
3003 FORMAT(12X,'Max contact stress yy is:   ',F15.4,' MPa')       
3004 FORMAT(12X,'Min contact stress yy is:   ',F15.4,' MPa')
3005 FORMAT(12X,'Max contact stress zz is:   ',F15.4,' MPa')       
3006 FORMAT(12X,'Min contact stress zz is:   ',F15.4,' MPa') 

2101 FORMAT(9X,'Max aperture is:   ',F14.5,' mm')       
2102 FORMAT(9X,'Min aperture is:   ',F14.5,' mm') 
3101 FORMAT(9X,'Max contact stress xx is:   ',F15.4,' MPa')       
3102 FORMAT(9X,'Min contact stress xx is:   ',F15.4,' MPa')
3103 FORMAT(9X,'Max contact stress yy is:   ',F15.4,' MPa')       
3104 FORMAT(9X,'Min contact stress yy is:   ',F15.4,' MPa')
3105 FORMAT(9X,'Max contact stress zz is:   ',F15.4,' MPa')       
3106 FORMAT(9X,'Min contact stress zz is:   ',F15.4,' MPa') 


if(i_NR_P>=2)then
  Old_Elem_Conta_Sta = c_Elem_Conta_Sta
  CT_State_Gauss(1:num_Crack,1:Max_Max_N_FluEl_3D)  =0
  c_Elem_Conta_Sta(1:Num_Elem,1:num_Crack) = 0
endif

PC_Gauss_x(1:num_Crack,1:Max_Max_N_FluEl_3D) = ZR 
PC_Gauss_y(1:num_Crack,1:Max_Max_N_FluEl_3D) = ZR
PC_Gauss_z(1:num_Crack,1:Max_Max_N_FluEl_3D) = ZR
all_Apertures(1:num_Crack,1:Max_Max_N_FluEl_3D) = ZR



!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_C,i_CT_Elem,c_Center, &
!$OMP            ori_n,T,tran_T,Solid_El,c_Aperture,c_Kt1,c_Kt2,c_Kn, &
!$OMP            c_PN,c_PT1,c_PT2,PC_Gauss,Relative_Disp)  &
!$OMP            SCHEDULE(STATIC)   
do i_C=1,num_Crack
  do i_CT_Elem = 1,Cracks_FluidEle_num_3D(i_C)
      c_Center = Cracks_FluidEle_Centroid_3D(i_C)%row(i_CT_Elem,1:3)
      ori_n  = Cracks_FluidEle_Vector_3D(i_C)%row(i_CT_Elem,1:3)
      T=Cracks_FluidEle_LCS_T_3D(i_C)%row(i_CT_Elem,1:3,1:3)
      tran_T = transpose(T)
      Solid_El = Cracks_FluidEle_EleNum_3D(i_C)%row(i_CT_Elem)   
      
      if(Key_Crack_Aperture_Method==1) then
          call Cal_Crack_Point_Aperture_3D(c_DISP,i_C,c_Center,Relative_Disp,i_CT_Elem,0,0,0)
          c_Aperture   = (Relative_Disp(1)*ori_n(1)+ Relative_Disp(2)*ori_n(2)+ Relative_Disp(3)*ori_n(3))/1.0D0
      elseif(Key_Crack_Aperture_Method==2) then
          call Cal_Point_Aperture_3D(i_C,c_Center,c_DISP,ori_n,T,c_Aperture,10)
      endif
      
      all_Apertures(i_C,i_CT_Elem) = c_Aperture
      
      Cracks_FluidEle_Aper_3D(i_C)%row(i_CT_Elem) = c_Aperture
      
      
      
      c_Kt1 = Kt1_Gauss(i_C,i_CT_Elem) 
      c_Kt2 = Kt2_Gauss(i_C,i_CT_Elem)       
      c_Kn  = Kn_Gauss(i_C,i_CT_Elem) 

      
      if(c_Aperture <= ZR)then
          
          
          c_PN = Kn*c_Aperture
          
          
          c_PT1 = ZR
          c_PT2 = ZR
          
          
          if(i_NR_P>=2)then
              CT_State_Gauss(i_C,i_CT_Elem) = 1
              !$OMP CRITICAL
              if(old_Elem_Conta_Sta(Solid_El,i_C)/=2) then
                  c_Elem_Conta_Sta(Solid_El,i_C) = 1
              else
                  c_Elem_Conta_Sta(Solid_El,i_C) = 2
              endif
              !$OMP END CRITICAL
          endif
      
          PC_Gauss(1:3) =MATMUL(tran_T,[c_PT1,c_PT2,c_PN])
          PC_Gauss_x(i_C,i_CT_Elem) =PC_Gauss(1)
          PC_Gauss_y(i_C,i_CT_Elem) =PC_Gauss(2)
          PC_Gauss_z(i_C,i_CT_Elem) =PC_Gauss(3)
          
          
          Kn_Gauss(i_C,i_CT_Elem)  = c_Kn
          Kt1_Gauss(i_C,i_CT_Elem) = c_Kt1
          Kt2_Gauss(i_C,i_CT_Elem) = c_Kt2
      else
          CT_State_Gauss(i_C,i_CT_Elem) = 0
          PC_Gauss_x(i_C,i_CT_Elem) =ZR
          PC_Gauss_y(i_C,i_CT_Elem) =ZR
          PC_Gauss_z(i_C,i_CT_Elem) =ZR                      
          Kn_Gauss(i_C,i_CT_Elem)  = c_Kn
          Kt1_Gauss(i_C,i_CT_Elem) = c_Kt1
          Kt2_Gauss(i_C,i_CT_Elem) = c_Kt2
      endif
  enddo
enddo
!$omp end parallel do        

if(i_NR_P>=2)then
    if (Num_Indent==8) then
        WRITE(*,2001) maxval(all_Apertures)*1000.0D0
        WRITE(*,2002) minval(all_Apertures)*1000.0D0
        if(minval(PC_Gauss_x)<ZR)then
            WRITE(*,3001) maxval(PC_Gauss_x,mask=(PC_Gauss_x<ZR))/1.0D6
            WRITE(*,3002) minval(PC_Gauss_x,mask=(PC_Gauss_x<ZR))/1.0D6
        else
            WRITE(*,3001) maxval(PC_Gauss_x)/1.0D6
            WRITE(*,3002) minval(PC_Gauss_x)/1.0D6        
        endif
        if(minval(PC_Gauss_y)<ZR)then
            WRITE(*,3003) maxval(PC_Gauss_y,mask=(PC_Gauss_y<ZR))/1.0D6
            WRITE(*,3004) minval(PC_Gauss_y,mask=(PC_Gauss_y<ZR))/1.0D6
        else
            WRITE(*,3003) maxval(PC_Gauss_y)/1.0D6
            WRITE(*,3004) minval(PC_Gauss_y)/1.0D6           
        endif
        if(minval(PC_Gauss_z)<ZR)then
            WRITE(*,3005) maxval(PC_Gauss_z,mask=(PC_Gauss_z<ZR))/1.0D6
            WRITE(*,3006) minval(PC_Gauss_z,mask=(PC_Gauss_z<ZR))/1.0D6
        else
            WRITE(*,3005) maxval(PC_Gauss_z)/1.0D6
            WRITE(*,3006) minval(PC_Gauss_z)/1.0D6          
        endif    
    elseif (Num_Indent==5) then
        WRITE(*,2101) maxval(all_Apertures)*1000.0D0
        WRITE(*,2102) minval(all_Apertures)*1000.0D0
        if(minval(PC_Gauss_x)<ZR)then
            WRITE(*,3101) maxval(PC_Gauss_x,mask=(PC_Gauss_x<ZR))/1.0D6
            WRITE(*,3102) minval(PC_Gauss_x,mask=(PC_Gauss_x<ZR))/1.0D6
        else
            WRITE(*,3101) maxval(PC_Gauss_x)/1.0D6
            WRITE(*,3102) minval(PC_Gauss_x)/1.0D6        
        endif
        if(minval(PC_Gauss_y)<ZR)then
            WRITE(*,3103) maxval(PC_Gauss_y,mask=(PC_Gauss_y<ZR))/1.0D6
            WRITE(*,3104) minval(PC_Gauss_y,mask=(PC_Gauss_y<ZR))/1.0D6
        else
            WRITE(*,3103) maxval(PC_Gauss_y)/1.0D6
            WRITE(*,3104) minval(PC_Gauss_y)/1.0D6           
        endif
        if(minval(PC_Gauss_z)<ZR)then
            WRITE(*,3105) maxval(PC_Gauss_z,mask=(PC_Gauss_z<ZR))/1.0D6
            WRITE(*,3106) minval(PC_Gauss_z,mask=(PC_Gauss_z<ZR))/1.0D6
        else
            WRITE(*,3105) maxval(PC_Gauss_z)/1.0D6
            WRITE(*,3106) minval(PC_Gauss_z)/1.0D6          
        endif    
    endif
endif
return 
end SUBROUTINE Cal_Contact_PN_and_PT_3D         
