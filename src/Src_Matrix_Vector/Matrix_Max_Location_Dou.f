 
      SUBROUTINE Matrix_Max_Location_Dou(m,n,Matrix,pos_i,pos_j)

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



