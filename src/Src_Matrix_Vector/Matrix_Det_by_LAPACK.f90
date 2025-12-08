 
SUBROUTINE Matrix_Det_by_LAPACK(n,Matrix,det,INFO)   

use Global_Float_Type      
use Global_Common
use Global_Model
use Global_Crack
use Global_Crack_3D
use Global_Crack_Common 

implicit none
integer n
real(kind=FT),intent(in)::Matrix(n,n)
real(kind=FT),intent(out)::det
integer,intent(out)::INFO

integer ipiv(n),lda
real(kind=FT) A(n,n)
real(kind=FT) det_P,det_L,det_U
integer i

integer Corresp_Dof,Corresp_Crack,Corresp_Node,Corresp_Direction
integer i_N,i_C



A = Matrix

INFO = 0

lda = n
call dgetrf( n, n, A,lda, ipiv, info)
if (info>0) then
       print *,'                 ZERO det found!'
       print *,'                 Return info:',INFO
       print *,'                 Look into matrix lines around',INFO
      det = ZR
endif
if (info<0) then
      print *,'    ERROR-2022081801 :: in Matrix_Det_by_LAPACK.f'
      print *,'                        Failed to get det!'          
      print *,'                        return info:',INFO
      call Warning_Message('S',Keywords_Blank)
endif    
  

det_U = ONE
do i=1,n
      det_U =  det_U*A(i,i)
enddo

det_L = ONE


det_P = ONE
do i=1,n
      if (ipiv(i)/=i)then
          det_P = -det_P
      endif
enddo

det = det_P*det_L*det_U

print *,'                 Determinant of K:', real(det)


101 FORMAT(18X,'Corresponding node number:  ',I7)     
102 FORMAT(18X,'Corresponding crack number: ',I7)     
103 FORMAT(18X,'Corresponding direction:    ',I7)     
104 FORMAT(18X,'Node enrichment type:       ',I7)    
105 FORMAT(18X,'Dis_Node_to_FS:             ',F12.5,' m')    
106 FORMAT(18X,'Surounding elements of node:',8I7)  
if(Key_Dimension==3) then
      if(info>=1)then
          print *,'                 sum(abs(K(det_INFO,  :))): ',real(sum(abs(Matrix(info,  1:n))))
          print *,'                 sum(abs(K(det_INFO-2,:))): ',real(sum(abs(Matrix(info-2,1:n))))                           
          print *,'                 sum(abs(K(det_INFO-1,:))): ',real(sum(abs(Matrix(info-1,1:n))))
          print *,'                 sum(abs(K(det_INFO+1,:))): ',real(sum(abs(Matrix(info+1,1:n))))
          print *,'                 sum(abs(K(det_INFO+2,:))): ',real(sum(abs(Matrix(info+2,1:n))))          
          Corresp_Dof = freeDOF_for_Det(info)
          print *,'                 Corresponding DOF number:  ',Corresp_Dof
          
          Corresp_Node  = -999
          Corresp_Crack = -999
          do i_N = 1,num_Node
              do i_C = 1,num_Crack
                  if(c_POS_3D(i_N,i_C)*3 == Corresp_Dof)then
                      Corresp_Node  = i_N 
                      Corresp_Crack = i_C
                      Corresp_Direction = 3
                      goto 665
                  endif
                  if((c_POS_3D(i_N,i_C)*3-1) == Corresp_Dof)then
                      Corresp_Node  = i_N 
                      Corresp_Crack = i_C
                      Corresp_Direction = 2
                      goto 665
                  endif    
                  if((c_POS_3D(i_N,i_C)*3-2) == Corresp_Dof)then
                      Corresp_Node  = i_N 
                      Corresp_Crack = i_C
                      Corresp_Direction = 1
                      goto 665
                  endif                     
              enddo
          enddo
          
          665 continue   
          write(*,101) Corresp_Node
          write(*,102) Corresp_Crack
          write(*,103) Corresp_Direction  
          write(*,104) Enriched_Node_Type_3D(i_N,i_C)  
          write(*,105) Dis_Node_to_FS(i_N)%row(i_C)
          write(*,106) Node_Elements(i_N,1:num_Node_Elements(i_N))
          

      endif
endif    
      
return
END SUBROUTINE Matrix_Det_by_LAPACK



