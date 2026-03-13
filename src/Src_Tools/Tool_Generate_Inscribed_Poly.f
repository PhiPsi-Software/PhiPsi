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
 
      subroutine Tool_Generate_Inscribed_Poly(
     &                 M_x_min,M_x_max,M_y_min,M_y_max,
     &                 c_Ave_R,c_Delta_R,
     &                 num_Rand_Poly_Incl,num_Vert_Poly_Incl,
     %                 check_R,
     &                 Max_Dis_to_Edge,
     &                 Poly_Incl_Coor_x,
     &                 Poly_Incl_Coor_y)
c     Randomly generate a circle, then generate an inscribed polygon within the circle
c     Algorithm: Generate one by one, and the midpoint of the circle generated later should not fall within the radius of the midpoints of other cracks.
c     inside the circle of radius r

      !***************************
      ! Variable Type Declaration
      !***************************
      use Global_Float_Type
      use Global_Common
      use Global_Model
      
      implicit none
      integer,intent(in)::num_Rand_Poly_Incl
      integer,intent(in)::num_Vert_Poly_Incl
      real(kind=FT),intent(in)::M_x_min,M_x_max,M_y_min,M_y_max
      real(kind=FT),intent(in)::c_Ave_R
      real(kind=FT),intent(in)::c_Delta_R
      real(kind=FT),intent(in)::check_R
      real(kind=FT),intent(out)::  Poly_Incl_Coor_x(num_Rand_Poly_Incl,
     &                                              num_Vert_Poly_Incl)   
      real(kind=FT),intent(out)::  Poly_Incl_Coor_y(num_Rand_Poly_Incl,
     &                                              num_Vert_Poly_Incl)   
      real(kind=FT),intent(in)::Max_Dis_to_Edge
      
      integer num_Cir,i_Try,i_Check
      real(kind=FT) c_x,c_y,Ran_Num,check_x,check_y,tem_Dis
      logical Yes_In_Circle
      real(kind=FT) U_1,U_2,UUU1
      real(kind=FT)  try_R,try_x,try_y
      logical Yes_In_Model_x,Yes_In_Model_y
      real(kind=FT) c_Center_Cir(num_Rand_Poly_Incl,2),
     &              c_R_Cir(num_Rand_Poly_Incl)
      integer i_Incl,i_point
      real(kind=FT) delta_Theta,Start_Theta
      print *,'    Generating random polygon......'
      
      Poly_Incl_Coor_x(1:num_Rand_Poly_Incl,1:num_Vert_Poly_Incl)  = ZR
      Poly_Incl_Coor_y(1:num_Rand_Poly_Incl,1:num_Vert_Poly_Incl)  = ZR
      c_Center_Cir(1:num_Rand_Poly_Incl,1:2) = ZR
      c_R_Cir(1:num_Rand_Poly_Incl)           = ZR
      num_Cir = 0
      
      if(key_Random==1)then
          call Init_Random_Seed()
      endif
      
      do i_Try = 1,50000000
          if(key_Random==1 .or. key_Random==0)then
              call random_number (Ran_Num)
          elseif(key_Random==2 .or. key_Random==3)then
              call Tool_Generate_Random_Number(Ran_Num)
          endif
          
          
          c_x =  M_x_min +(M_x_max-M_x_min)*Ran_Num
          
          if(key_Random==1 .or. key_Random==0)then
              call random_number (Ran_Num)
          elseif(key_Random==2 .or. key_Random==3)then
              call Tool_Generate_Random_Number(Ran_Num)
          endif
          
          c_y =  M_y_min +(M_y_max-M_y_min)*Ran_Num
          Yes_In_Circle = .False.
          do i_Check=1,num_Cir
              check_x = c_Center_Cir(i_Check,1)
              check_y = c_Center_Cir(i_Check,2)
              tem_Dis = sqrt((c_x-check_x)**2 + (c_y-check_y)**2)
              if(tem_Dis<=check_R)then
                  Yes_In_Circle = .True.
                  exit
              endif
          enddo
          if(Yes_In_Circle .eqv. .False.)then
              if(key_Random==1 .or. key_Random==0)then
                  call random_number (U_1) 
                  call random_number (U_2)
              elseif(key_Random==2 .or. key_Random==3)then
                  call Tool_Generate_Random_Number(U_1)
                  call Tool_Generate_Random_Number(U_2)
              endif
              
              UUU1 = sqrt(-TWO*log(U_1)) * cos(TWO*PI*U_2)
              if(UUU1 >THR)UUU1=THR
              if(UUU1 <-THR)UUU1=-THR
              try_R = c_Ave_R + c_Delta_R*UUU1/THR
              Try_x = c_x ; Try_y = c_y 
              Yes_In_Model_x = .False.
              Yes_In_Model_y = .False.
              if (Try_x+try_R<M_x_max-Max_Dis_to_Edge .and.
     &            Try_x-try_R > M_x_min+Max_Dis_to_Edge) then
                  Yes_In_Model_x = .True.
              endif
              if (Try_y+try_R<M_y_max -Max_Dis_to_Edge .and.
     &            Try_y-try_R > M_y_min +Max_Dis_to_Edge) then
                  Yes_In_Model_y = .True.
              endif
              if(Yes_In_Model_x .and. Yes_In_Model_y)then
                  num_Cir = num_Cir +1
                  c_Center_Cir(num_Cir,1) = c_x 
                  c_Center_Cir(num_Cir,2) = c_y
                  c_R_Cir(num_Cir)        = try_R
              endif
          endif
          if(num_Cir>=num_Rand_Poly_Incl)then
              exit
          endif
      enddo
      if(num_Cir < num_Rand_Poly_Incl)then
          print *,'    Error :: Inscribed polygons generation failed!'
          call Warning_Message('S',Keywords_Blank)
      endif
      
      
      delta_Theta = TWO*pi/num_Vert_Poly_Incl
      do i_Incl =1,num_Rand_Poly_Incl
          if(key_Random==1 .or. key_Random==0)then
              call random_number (Ran_Num)
          elseif(key_Random==2 .or. key_Random==3)then
              call Tool_Generate_Random_Number(Ran_Num)
          endif
          
          Start_Theta =  TWO*pi*Ran_Num
          do i_Point =1,num_Vert_Poly_Incl
              Poly_Incl_Coor_x(i_Incl,i_Point) = 
     &          c_Center_Cir(i_Incl,1) + 
     &          c_R_Cir(i_Incl)*cos(Start_Theta+delta_Theta*(i_Point-1))
              Poly_Incl_Coor_y(i_Incl,i_Point) = 
     &          c_Center_Cir(i_Incl,2) + 
     &          c_R_Cir(i_Incl)*sin(Start_Theta+delta_Theta*(i_Point-1))
          enddo
      enddo
      
      return 
      end SUBROUTINE Tool_Generate_Inscribed_Poly 
      

