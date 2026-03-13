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
 
      subroutine Tool_Length_CalP_to_InjP(i_Crack,Length_Out,
     &                                  CalP_x,CalP_y,i_CalP)
C     Calculate the distance from the point to the water injection point
c     Calculate the point number based on the water injection point (this method is very convenient)
      use Global_Float_Type      
      use Global_Common
      use Global_Crack
      use Global_Crack_Common
      use Global_HF
      
      implicit none
      integer,intent(in)::i_Crack,i_CalP
      real(kind=FT),intent(in)::CalP_x,CalP_y
      real(kind=FT),intent(out)::Length_Out
      real(kind=FT) delta_L
      integer i
      real(kind=FT) CalP_L_x,CalP_L_y,CalP_R_x,CalP_R_y
      
      if(i_Crack>num_Crack)then
          print *,'     Error:: wrong crack number for '
     &            // 'Tool_Length_of_Crack.f!'
          call Warning_Message('S',Keywords_Blank) 
      endif
      
      Length_Out = ZR
      
      if(Key_Symm_HF==0)then
          if(i_CalP>CalP_num_InjP_Local)then
              do i=CalP_num_InjP_Local,i_CalP-1
                  CalP_L_x=Cracks_CalP_Coors(i_Crack,i,1)
                  CalP_L_y=Cracks_CalP_Coors(i_Crack,i,2)
                  CalP_R_x=Cracks_CalP_Coors(i_Crack,i+1,1)
                  CalP_R_y=Cracks_CalP_Coors(i_Crack,i+1,2)
                  delta_L=sqrt((CalP_L_x-CalP_R_x)**2 +
     &                         (CalP_L_y-CalP_R_y)**2) 
                  Length_Out = Length_Out + delta_L
              end do
          elseif(i_CalP<CalP_num_InjP_Local)then
              do i=i_CalP,CalP_num_InjP_Local-1
                  CalP_L_x=Cracks_CalP_Coors(i_Crack,i,1)
                  CalP_L_y=Cracks_CalP_Coors(i_Crack,i,2)
                  CalP_R_x=Cracks_CalP_Coors(i_Crack,i+1,1)
                  CalP_R_y=Cracks_CalP_Coors(i_Crack,i+1,2)
                  delta_L=sqrt((CalP_L_x-CalP_R_x)**2 +
     &                         (CalP_L_y-CalP_R_y)**2) 
                  Length_Out = Length_Out + delta_L
              end do
          else
              Length_Out=ZR
          end if
      end if
      if(Key_Symm_HF==1)then
          if(i_CalP>0)then
              do i=1,i_CalP-1
                  CalP_L_x=Cracks_CalP_Coors(i_Crack,i,1)
                  CalP_L_y=Cracks_CalP_Coors(i_Crack,i,2)
                  CalP_R_x=Cracks_CalP_Coors(i_Crack,i+1,1)
                  CalP_R_y=Cracks_CalP_Coors(i_Crack,i+1,2)
                  delta_L=sqrt((CalP_L_x-CalP_R_x)**2 +
     &                         (CalP_L_y-CalP_R_y)**2) 
                  Length_Out = Length_Out + delta_L
              end do
          else
              Length_Out=ZR
          end if
      end if
      
      return 
      end SUBROUTINE Tool_Length_CalP_to_InjP                         
