 
      subroutine Cal_HF_Resid(ifra,iter,Counter_Iter,R,Total_FD,
     &                              num_FreeD,
     &                              num_free_CalP,
     &                              globalK,Coupled_Q,F_U,
     &                              Last_DISP,Last_Last_DISP,
     &                              freeDOF,freeDOF_HF,Local_freeDOF_HF,
     &                              Last_CalP_Pres,delta_Time,
     &                              H,c_S)
      
      use Global_Float_Type
      use Global_Crack
      use Global_Crack_Common
      use Global_HF
      use Global_Material
      
      implicit none
      integer, intent(in)::ifra,iter,Counter_Iter,Total_FD,num_FreeD,
     &                     num_free_CalP
      real(kind=FT),intent(in)::
     &              Coupled_Q(Total_FD,num_Tol_CalP_Water),
     &              globalK(Total_FD,Total_FD),
     &              F_U(Total_FD),Last_DISP(Total_FD),
     &              Last_Last_DISP(Total_FD),
     &              Last_CalP_Pres(num_Tol_CalP_Water),delta_Time,
     &              H(num_Tol_CalP_Water,num_Tol_CalP_Water),
     &              c_S(num_Tol_CalP_Water)
      integer, intent(in)::freeDOF(Total_FD),
     &                     freeDOF_HF(num_Tol_CalP_Water),
     &                     Local_freeDOF_HF(num_Tol_CalP_Water)
      real(kind=FT),intent(out)::R(Total_FD+num_Tol_CalP_Water)
      real(kind=FT) delta_DISP(Total_FD)
      
      R(1:Total_FD+num_Tol_CalP_Water) = ZR

      R(freeDOF(1:num_FreeD))  =
     &    MATMUL(globalK(freeDOF(1:num_FreeD),freeDOF(1:num_FreeD)),
     &            Last_DISP(freeDOF(1:num_FreeD)))
     &   -MATMUL(Coupled_Q(freeDOF(1:num_FreeD),
     &                     Local_freeDOF_HF(1:num_free_CalP)),
     &  thick(1)*Last_CalP_Pres(Local_freeDOF_HF(1:num_free_CalP)))     
     &   -F_U(freeDOF(1:num_FreeD))

      if(Key_IniPre_PassOn==0) then
          if (iter == 1) then  
              delta_DISP = Last_DISP
          elseif (iter > 1) then
              delta_DISP(freeDOF(1:num_FreeD)) = 
     &                       Last_DISP(freeDOF(1:num_FreeD))- 
     &                       Last_Last_DISP(freeDOF(1:num_FreeD))
          end if
      elseif(Key_IniPre_PassOn==1) then
          if (iter == 1) then  
              if (ifra == 1) then
                  delta_DISP = Last_DISP
              elseif(ifra > 1) then
                  delta_DISP = Last_DISP
              endif          
          elseif (iter > 1) then
              delta_DISP(freeDOF(1:num_FreeD)) = 
     &                       Last_DISP(freeDOF(1:num_FreeD))- 
     &                       Last_Last_DISP(freeDOF(1:num_FreeD))
          end if              
      endif      
      R(freeDOF_HF(1:num_free_CalP)) = -MATMUL(transpose(
     &    Coupled_Q(freeDOF(1:num_FreeD),
     &              Local_freeDOF_HF(1:num_free_CalP))),
     &    delta_DISP(freeDOF(1:num_FreeD)))
     &    - delta_Time*MATMUL(H(Local_freeDOF_HF(1:num_free_CalP),
     &               Local_freeDOF_HF(1:num_free_CalP)),
     &      Last_CalP_Pres(Local_freeDOF_HF(1:num_free_CalP)))
     &    + c_S(Local_freeDOF_HF(1:num_free_CalP))*delta_Time
      return 
      end SUBROUTINE Cal_HF_Resid          
