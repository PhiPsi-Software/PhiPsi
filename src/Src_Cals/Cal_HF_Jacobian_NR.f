 
      subroutine Cal_HF_Jacobian_NR(Counter,NR_Deri,Total_FD,
     &                              num_FreeD,
     &                              num_free_CalP,
     &                              globalK,Coupled_Q,
     &                              freeDOF,freeDOF_HF,
     &                              Local_freeDOF_HF,delta_Time,H)   
      
      use Global_Float_Type
      use Global_Crack
      use Global_Crack_Common

      implicit none
      integer, intent(in)::Counter,Total_FD,num_FreeD,
     &                     num_free_CalP
      real(kind=FT),intent(in)::
     &              Coupled_Q(Total_FD,num_Tol_CalP_Water),
     &              globalK(Total_FD,Total_FD),delta_Time,
     &              H(num_Tol_CalP_Water,num_Tol_CalP_Water)
      integer, intent(in)::freeDOF(Total_FD),
     &                     freeDOF_HF(num_Tol_CalP_Water),
     &                     Local_freeDOF_HF(num_Tol_CalP_Water)
      real(kind=FT),intent(out)::
     &  NR_Deri(Total_FD+num_Tol_CalP_Water,Total_FD+num_Tol_CalP_Water)
      print *,'    Calculating the Jacobian matrix of ' 
     &             //'NR interation......'
      
      NR_Deri(1:Total_FD+num_Tol_CalP_Water,
     &        1:Total_FD+num_Tol_CalP_Water)=ZR
      
      NR_Deri(freeDOF(1:num_FreeD),freeDOF(1:num_FreeD)) = 
     &          globalK(freeDOF(1:num_FreeD),freeDOF(1:num_FreeD))
      NR_Deri(freeDOF(1:num_FreeD),freeDOF_HF(1:num_free_CalP)) = 
     &           -Coupled_Q(freeDOF(1:num_FreeD),
     &                      Local_freeDOF_HF(1:num_free_CalP))
      NR_Deri(freeDOF_HF(1:num_free_CalP),freeDOF(1:num_FreeD)) = 
     &         -transpose(Coupled_Q(freeDOF(1:num_FreeD),
     &                             Local_freeDOF_HF(1:num_free_CalP)))
      NR_Deri(freeDOF_HF(1:num_free_CalP),freeDOF_HF(1:num_free_CalP))=
     &          -H(Local_freeDOF_HF(1:num_free_CalP),
     &             Local_freeDOF_HF(1:num_free_CalP))*delta_Time
      
      return 
      end SUBROUTINE Cal_HF_Jacobian_NR          
