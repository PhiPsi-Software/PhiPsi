 
      function MVMUL_LP(A,V) result (Output)  
      
      use Global_Float_Type
      implicit none 
      real(kind=FT),dimension(:,:),intent(in):: A
      real(kind=FT),dimension(:),  intent(in):: V
      real(kind=FT),dimension(size(A,1))::Output
      integer::m,n
      
      m = size(A,1)
      n = size(V)
      
      if(FT==8)then
          call dgemv('N', m, n, ONE, A, m, V, 1, ZR, Output, 1)
      elseif(FT==4)then
          call sgemv('N', m, n, ONE, A, m, V, 1, ZR, Output, 1)
      endif
      END function MVMUL_LP
    


