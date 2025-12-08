 
      subroutine Tool_Fit_3D_Points_to_Straight_Line(num_Point,
     &                                    In_Points,Out_Points)
      use Global_Float_Type     
      use Global_Model
      
      implicit none
      integer,intent(in)::num_Point
      real(kind=FT),intent(in)::In_Points(num_Point,3)
      real(kind=FT),intent(out)::Out_Points(num_Point,3)
      real(kind=FT) Ave_Location(3)
      real(kind=FT) var_X,var_Y,var_Z
      real(kind=FT) cov_XY,cov_YZ,cov_XZ,cov_Matrix(3,3)
      real(kind=FT) abserr,Eigen_Vector(3,3),Eigen_Value(3)
      integer max_eigen_value_number
      real(kind=FT) Line_n_Vector(1:3),Line_P1(3),Line_P2(3)
      real(kind=FT) Nearest_Point(3),c_dist,c_t
      
      
      integer i_P
      
      
      Ave_Location(1) = sum(In_Points(1:num_Point,1))/dble(num_Point)
      Ave_Location(2) = sum(In_Points(1:num_Point,2))/dble(num_Point)
      Ave_Location(3) = sum(In_Points(1:num_Point,3))/dble(num_Point)
      
      var_X=ZR
      var_Y=ZR
      var_Z=ZR
      cov_XY=ZR
      cov_YZ=ZR
      cov_XZ=ZR
      do i_P = 1,num_Point
          var_X = var_X + (In_Points(i_P,1) -Ave_Location(1))**2
          var_Y = var_Y + (In_Points(i_P,2) -Ave_Location(2))**2
          var_Z = var_Z + (In_Points(i_P,3) -Ave_Location(3))**2
          cov_XY= cov_XY + (In_Points(i_P,1) -Ave_Location(1))*
     &                     (In_Points(i_P,2) -Ave_Location(2))
          cov_YZ= cov_YZ + (In_Points(i_P,2) -Ave_Location(2))*
     &                     (In_Points(i_P,3) -Ave_Location(3))
          cov_XZ= cov_XZ + (In_Points(i_P,1) -Ave_Location(1))*
     &                     (In_Points(i_P,3) -Ave_Location(3))     
      enddo
      
      var_X = var_X/dble(num_Point-1)
      var_Y = var_Y/dble(num_Point-1)
      var_Z = var_Z/dble(num_Point-1)
      cov_XY = cov_XY/dble(num_Point-1)
      cov_YZ = cov_YZ/dble(num_Point-1)
      cov_XZ = cov_XZ/dble(num_Point-1)
      
      cov_Matrix(1,1:3) = [var_X,  cov_XY, cov_XZ]
      cov_Matrix(2,1:3) = [cov_XY, var_Y,  cov_YZ]
      cov_Matrix(3,1:3) = [cov_XZ, cov_YZ,  var_Z]
      
      abserr = Tol_9
      call Matrix_Eigenvalues_and_Eigenvectors(cov_Matrix(1:3,1:3),
     &                                         Eigen_Vector,abserr,3)
      Eigen_Value(1) = cov_Matrix(1,1)
      Eigen_Value(2) = cov_Matrix(2,2)
      Eigen_Value(3) = cov_Matrix(3,3)
      
      
      max_eigen_value_number = MAXLOC(Eigen_Value,1)
      
      Line_n_Vector(1:3) = Eigen_Vector(1:3,max_eigen_value_number)
      Line_P1  = Ave_Location + Max_Model_Range*Line_n_Vector
      Line_P2  = Ave_Location - Max_Model_Range*Line_n_Vector
      
      do i_P = 1,num_Point
          call line_exp_point_near_3d(Line_P1,Line_P2,
     &                                In_Points(i_P,1:3),
     &                                Nearest_Point(1:3), c_dist, c_t)
          Out_Points(i_P,1:3) = Nearest_Point(1:3)
      enddo
      
      return 
      end SUBROUTINE Tool_Fit_3D_Points_to_Straight_Line                  
