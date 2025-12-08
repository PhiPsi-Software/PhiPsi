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
 
      subroutine D3_Get_Max_Min_Ave_KIeq_of_HF_Cracks(
     &                    Max_KI_eq_3D,Min_KI_eq_3D,Ave_KI_eq_3D,
     &                    c_Crack_Max_Num,c_Vertex_Max_Num)
      ! Extract the maximum, minimum, and average equivalent stress intensity factors of the HF fracture.
      ! Called by the main program PhiPsi3D_Static_HF_SlipWater.f.
      !2022-07-05.   
      !Updated 2023-05-16.
      
c     ============================
c     Read public variable module
c     ============================
      use Global_Float_Type                                                                
      use Global_Crack_Common
      use Global_Crack_3D
      use Global_Ragged_Array_Real_Classs
      
c     =====================
c     Variable Declaration
c     =====================
      implicit none
      real(kind=FT),intent(out)::Max_KI_eq_3D,Min_KI_eq_3D,Ave_KI_eq_3D
      integer,intent(out)::c_Crack_Max_Num,c_Vertex_Max_Num
      integer i_C
      integer num_Cr_Edges
      real(kind=FT) c_Max,c_Min, Sum_KIeq
      integer Ver_Count
      integer c_Max_location
      
c     ===============
c     Fracture cycle
c     ===============
      Max_KI_eq_3D = -Con_Big_15   
      Min_KI_eq_3D =  Con_Big_15
      Ver_Count = 0
      Sum_KIeq  = ZR
      
      do i_C = 1,num_Crack
          num_Cr_Edges = Crack3D_Meshed_Outline_num(i_C)
          ! If it is an HF crack
          if(Crack_Type_Status_3D(i_C,1)==1) then 
              c_Max = maxval(KI_eq_3D(i_C)%row(1:num_Cr_Edges))
              c_Max_location=MAXLOC(KI_eq_3D(i_C)%row(1:num_Cr_Edges),1)
              c_Min = minval(KI_eq_3D(i_C)%row(1:num_Cr_Edges))
              
              if(c_Max > Max_KI_eq_3D) then 
                  Max_KI_eq_3D = c_Max
                  c_Crack_Max_Num = i_C
                  c_Vertex_Max_Num= c_Max_location
              endif
              if(c_Min < Min_KI_eq_3D) then 
                  Min_KI_eq_3D = c_Min
              endif
              Ver_Count = Ver_Count +num_Cr_Edges
              Sum_KIeq  = Sum_KIeq + 
     &                    sum(KI_eq_3D(i_C)%row(1:num_Cr_Edges))
          endif
      enddo
      
      Ave_KI_eq_3D = Sum_KIeq/Ver_Count
      
      return 
      end SUBROUTINE D3_Get_Max_Min_Ave_KIeq_of_HF_Cracks
