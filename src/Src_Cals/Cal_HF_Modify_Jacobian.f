 
      subroutine Cal_HF_Modify_Jacobian(Counter,NR_Deri,Total_FD,
     &                              num_free_CalP,
     &                              freeDOF_HF,
     &                              Local_freeDOF_HF,delta_Time,H)   
      
      use Global_Float_Type
      use Global_Crack
      use Global_Crack_Common

      implicit none
      integer, intent(in)::Counter,Total_FD,
     &                     num_free_CalP
      real(kind=FT),intent(in)::H(num_Tol_CalP_Water,
     &                             num_Tol_CalP_Water),
     &                             delta_Time
      integer, intent(in)::freeDOF_HF(num_Tol_CalP_Water),
     &                     Local_freeDOF_HF(num_Tol_CalP_Water)
      real(kind=FT),intent(inout)::
     &  NR_Deri(Total_FD+num_Tol_CalP_Water,Total_FD+num_Tol_CalP_Water)
      real(kind=FT) k,max_K(3)
      
      integer i
      
      print *,'    Modifying the Jacobian matrix......'
     
      NR_Deri(freeDOF_HF(1:num_free_CalP),freeDOF_HF(1:num_free_CalP))=
     &          -H(Local_freeDOF_HF(1:num_free_CalP),
     &             Local_freeDOF_HF(1:num_free_CalP))*delta_Time
     
      return 
      end SUBROUTINE Cal_HF_Modify_Jacobian      
