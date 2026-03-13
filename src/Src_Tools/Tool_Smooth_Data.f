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
 
      subroutine Tool_Smooth_Data(Data_Points,num_points,Smooth_method,
     &                            Moving_Average_Method_n,logical_Close)
c     Used for smoothing data. One-dimensional curve data.
c     Added on 2022-04-25.
c     Smooth_method = 1, moving average method, using the values from previous and subsequent moments
c     A total of 2n + 1 values are averaged (\theory_documents\029 Moving Average Smoothing Method-2022-04-25.docx).
c     

      !......................
      ! Variable Declaration
      !......................
      use Global_Float_Type
      use Global_Common
      implicit none
      integer,intent(in)::num_points,Smooth_method
      integer,intent(in)::Moving_Average_Method_n
      logical,intent(in)::logical_Close
      real(kind=FT),intent(inout)::Data_Points(num_points)
      real(kind=FT) in_Data(num_points),out_Data(num_points)
      integer n
      integer i,j,k
      integer twice_n1
      real(kind=FT),ALLOCATABLE::tem_Data(:)
      real(kind=FT) c_average
      
      in_Data = Data_Points
      n       = Moving_Average_Method_n
      
      
      twice_n1 = 2*n+1
      
      if(twice_n1 >= num_points) then
        print *,'    Error :: n is too big or num_points is too small!'
        print *,'             in Tool_Smooth_Data.f!'
        print *,'             n:',n
        print *,'             2*n+1:',2*n+1
        print *,'             num_points:',num_points
        call Warning_Message('S',Keywords_Blank)   
      endif
      
      if (Smooth_method==1) then
          ALLOCATE(tem_Data(-num_points:2*num_points))
          if (logical_Close .eqv. .False.) then
              tem_Data(-num_points:0) = ZR
              tem_Data(1:num_points)  = in_Data(1:num_points)
              tem_Data((num_points+1):2*num_points) = ZR
          elseif (logical_Close .eqv. .True.) then
              tem_Data(-num_points+1:0)  = in_Data(1:num_points)
              tem_Data(1:num_points)     = in_Data(1:num_points)
              tem_Data((num_points+1):2*num_points) =  
     &                                   in_Data(1:num_points)
          endif
          do i=1,num_points
             c_average   = sum(tem_Data((i-n):(i+n)))/dble(twice_n1)
             out_Data(i) = c_average 
          enddo
          DEALLOCATE(tem_Data)
      endif
      
      Data_Points = out_Data
      
      return 
      end SUBROUTINE Tool_Smooth_Data             
