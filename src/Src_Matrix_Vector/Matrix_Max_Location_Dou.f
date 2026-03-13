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
 
      SUBROUTINE Matrix_Max_Location_Dou(m,n,Matrix,pos_i,pos_j)
C     Get the position of the maximum element in the matrix. Double precision.
c     2022-04-30.
c     NEWFTU2022043005.
c     Usage of the MAXLOC function: https://www.cita.utoronto.ca/~merz/intel_f10b/main_for/mergedProjects/lref_for/source_files/rfmaxloc.htm

      use Global_Float_Type
      implicit none
      integer,intent(in)::m,n
      real(kind=FT),intent(in)::Matrix(m,n)       
      integer,intent(out)::pos_i,pos_j
      integer Location_Cols(m)
      real(kind=FT) Max_Value_Cols(m)
      
      Max_Value_Cols = MAXVAL(Matrix, DIM=2)
      Location_Cols  = MAXLOC(Matrix, DIM=2)
      
      pos_i = MAXLOC(Max_Value_Cols, DIM=1)
      pos_j = Location_Cols(pos_i)

      return
      END SUBROUTINE Matrix_Max_Location_Dou



