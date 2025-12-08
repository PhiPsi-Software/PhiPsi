 
      SUBROUTINE Matrix_Count_Row_Int_m_x_2(m,Matrix,
     &                                      Count_Row,Num_Once,
     &                                      key_process)   

      use Global_Float_Type    
      use CMD_Progress
      use Global_Common
#ifndef Silverfrost
      use omp_lib
#endif 
      implicit none
      
      integer,intent(in)::m
      integer,intent(in)::key_process
      integer,intent(in)::Matrix(m,2)
      integer,intent(out)::Count_Row(m)
      integer,intent(out)::Num_Once
      
      integer i
      
      type(CLS_CMD_Progress)::Progress

      Count_Row(1:m) = 0
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(static)  
      do i=1,m    
          Count_Row(i) =count((Matrix(1:m,1) .EQ. Matrix(i,1)) .and.
     &                        (Matrix(1:m,2) .EQ. Matrix(i,2)))
      end do      
      !$omp end parallel do 
      
      Num_Once= 0
      Num_Once=count(Count_Row .eq. 1)
      
      return
      END SUBROUTINE Matrix_Count_Row_Int_m_x_2
    


