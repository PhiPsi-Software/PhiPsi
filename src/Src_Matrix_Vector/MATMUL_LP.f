 
      function MATMUL_LP(A,B) result (Output)  
      
      use Global_Float_Type
      implicit none 
      real(kind=FT),dimension(:,:),intent(in):: A,B
      real(kind=FT),dimension(size(A,1),size(B,2))::Output
      integer::m,k,n
      
      m = size(A,1)
      n = size(B,2)
      k = size(A,2)
      
      if(FT==8)then
          CALL DGEMM('N','N',m,n,k,ONE,A,m,B,k,ZR,Output,m)
      elseif(FT==4)then
          CALL SGEMM('N','N',m,n,k,ONE,A,m,B,k,ZR,Output,m)
      endif

      END function MATMUL_LP
    


