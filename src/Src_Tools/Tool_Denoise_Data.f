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
 
      subroutine Tool_Denoise_Data(Data_Points,num_points,
     &                            Denoise_method, n_Sigma,
     &                            logical_Close)   
     
c     Used for denoising data. One-dimensional curve data.
c     Added on 2022-04-26.
c     Denoise = 1, standard deviation method.
c                  Ref:https://medium.com/geekculture/curve-smoothing-and-outliers-removal-using-c-d3bf6e2fbc78 
c                  or
c     \theory_documents\030 Curve Denoising (Deburring) and Smoothing Treatment-2022-04-26.pdf

      !......................
      ! Variable Declaration
      !......................
      use Global_Float_Type
      use Global_Common
      implicit none
      integer,intent(in)::num_points,Denoise_method
      integer,intent(in)::n_Sigma
      logical,intent(in)::logical_Close
      real(kind=FT),intent(inout)::Data_Points(num_points)
      real(kind=FT) in_Data(num_points)
      integer i
      real(kind=FT) c_average,SD,c_x,low_limit,up_limit
      integer c_count
          
      in_Data = Data_Points

      
      if (Denoise_method==1) then
          c_average = sum(in_Data)/dble(num_points)
          call Tool_Standard_Deviation(in_Data,num_points,SD)
          
          
          c_count = 0 
          do i=1,num_points
             c_x = in_Data(i)
             low_limit = c_average - n_Sigma*SD
             up_limit  = c_average + n_Sigma*SD 
             if((c_x < low_limit) .or. (c_x > up_limit)) then
                 if (c_x > c_average ) then 
                     Data_Points(i) = c_average + SD
                 else
                     Data_Points(i) = c_average - SD
                 endif
                 c_count =  c_count +1
             endif
          enddo  
          
      endif
      
      return 
      end SUBROUTINE Tool_Denoise_Data             
