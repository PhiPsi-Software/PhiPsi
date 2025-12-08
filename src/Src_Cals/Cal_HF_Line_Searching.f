 
      subroutine Cal_HF_Line_Searching(ifra,iter,Counter_Iter,
     &     num_ALlDOF,c_Total_FD,num_FreeD,num_free_CalP,Total_freeDOF,
     &     c_num_Tol_CalP_Water,num_Total_FD,freeDOF_HF,
     &       Local_freeDOF_HF,NR_Deri,
     &       delta_x,
     &       c_R,Initial_DISP,
     &       c_DISP,c_freeDOF,c_CalP_Pres,c_globalK,Coupled_Q,F_U,
     &       c_Temp_Cr_CalP_Aper,c_Temp_total_Time,
     &       Num_Line_Search,total_Time,x)
      use Global_Float_Type
      use Global_Common
      use Global_Crack
      use Global_HF
      use Global_Material
      
      
      implicit none
      integer,intent(in)::ifra,iter,Counter_Iter,c_Total_FD,num_FreeD,
     &               num_free_CalP,num_ALlDOF,c_num_Tol_CalP_Water,
     &               num_Total_FD
      integer,intent(in)::Local_freeDOF_HF(c_num_Tol_CalP_Water)
      integer,intent(in)::Total_freeDOF(num_ALlDOF)
      real(kind=FT),intent(in)::c_R(num_ALlDOF)
      real(kind=FT),intent(in)::
     &              Coupled_Q(c_Total_FD,c_num_Tol_CalP_Water),
     &              c_globalK(c_Total_FD,c_Total_FD),F_U(c_Total_FD)
      real(kind=FT),intent(in)::delta_x(num_FreeD+num_free_CalP)
      real(kind=FT),intent(in)::c_DISP(c_Total_FD)
      real(kind=FT),intent(in)::Initial_DISP(c_Total_FD)
      integer,intent(in)::c_freeDOF(c_Total_FD)
      integer,intent(in)::freeDOF_HF(c_num_Tol_CalP_Water)
      real(kind=FT),intent(in)::c_CalP_Pres(c_num_Tol_CalP_Water)
      real(kind=FT),intent(in)::
     &               c_Temp_Cr_CalP_Aper(Max_Num_Cr,Max_Num_Cr_CalP)
      real(kind=FT),intent(in)::c_Temp_total_Time
      real(kind=FT),intent(in)::
     &  NR_Deri(c_Total_FD+c_num_Tol_CalP_Water,
     &          c_Total_FD+c_num_Tol_CalP_Water)
      integer,intent(inout)::Num_Line_Search
      real(kind=FT),intent(inout)::total_Time
      real(kind=FT),intent(out)::x(num_FreeD+num_free_CalP)
      real(kind=FT) Try_DISP(c_Total_FD)    
      real(kind=FT) Try_CalP_Pres(c_num_Tol_CalP_Water)    
      real(kind=FT) g(num_FreeD+num_free_CalP)
      real(kind=FT) p(num_FreeD+num_free_CalP)
      real(kind=FT) x_old(num_FreeD+num_free_CalP)
      real(kind=FT) f_vector(num_FreeD+num_free_CalP)
      real(kind=FT) H(c_num_Tol_CalP_Water,c_num_Tol_CalP_Water)
      real(kind=FT) c_S(c_num_Tol_CalP_Water)
      real(kind=FT) Last_DISP(c_Total_FD)
      real(kind=FT) F_c,F_old,F_2,c_Sum 
      integer i,j,c_i,c_j,Lnsrch_Count
      real(kind=FT) STPMX,stpmax,slope,alam,temp,test,alamin
      real(kind=FT) ALF,TOLX,a,b,disc,rhs1,rhs2,delta_Time,alam2
      real(kind=FT) tmplam
      real(kind=FT) c_R_new(num_ALlDOF)
      
 2022 FORMAT(11X,'>>> Line search and backtracking ',I2) 
 3001 FORMAT(11X,'value of F_old:    ', E13.6)
 3002 FORMAT(11X,'    value of F:    ', E13.6)
 
      ALF = 1.0D-4
      TOLX= 1.0D-7
      
      print *,'    Starting Line searches '
     &                               // 'and backtracking...'
      print *,'      ~~~~Prepare data for Line searches '
     &                               // 'and backtracking...'
      f_vector(1:num_FreeD+num_free_CalP) =
     &            c_R(Total_freeDOF(1:num_Total_FD))
      F_old = HLF*dot_product(f_vector,f_vector)
      write(*,3001) F_old

      x_old(1:num_FreeD) = c_DISP(c_freeDOF(1:num_FreeD)) 
      x_old(num_FreeD+1:num_FreeD+num_free_CalP) = 
     &        c_CalP_Pres(Local_freeDOF_HF(1:num_free_CalP))
      do i=1,num_Total_FD
          c_Sum = ZR
          do j=1,num_Total_FD
              c_i = Total_freeDOF(i)
              c_j = Total_freeDOF(j)
              c_Sum = c_Sum + NR_Deri(c_j,c_i)*f_vector(j)
          end do
          g(i) = c_Sum   
      end do
      p = -delta_x
      STPMX = 100.0D0
      c_Sum = DOT_PRODUCT(x_old,x_old)
      stpmax = STPMX*max(sqrt(c_Sum),real(num_Total_FD))
      
      print *,'      ~~~~Line searching and backtracking...'
      c_Sum = sqrt(DOT_PRODUCT(p,p))
      if(c_Sum > stpmax) then
          print *,"          !!!!!!!!!!!!!!!!!!!!!!!!!"
          print *,"          WARNING :: Sum > stpmax     "
          print *,"          !!!!!!!!!!!!!!!!!!!!!!!!!"
          p = p*stpmax/c_Sum    
      end if
      slope = DOT_PRODUCT(g,p)
      print *,'          Value of slope:',slope
      if (slope > ZR) then
          print *,'          Wrong slope value in line '
     &                                        //'searching process!'
          call Warning_Message('c_S',Keywords_Blank)     
      end if
      test = ZR
      do i=1,num_Total_FD
          temp = abs(p(i))/max(abs(x_old(i)),ONE)
          if (temp > test) test = temp
      end do
      alamin = TOLX/test
      print *,'          λ_min:',alamin
      Lnsrch_Count = 0
      alam = ONE
      
8888  continue
          Lnsrch_Count = Lnsrch_Count +1
          Num_Line_Search = Num_Line_Search +1
          write(*,2022) Lnsrch_Count
          if (Lnsrch_Count > Max_Num_Lnsrch)then
              print *,'          +++++++++++++++++++++++++'
              print *,'            Line searching failed!'
              print *,'          +++++++++++++++++++++++++'
              call Warning_Message('S',Keywords_Blank)  
              goto 9999
          end if
          x = x_old + alam*p
          Try_DISP(1:c_Total_FD)  = ZR
          Try_DISP(c_freeDOF(1:num_FreeD)) = x(1:num_FreeD)
          Try_CalP_Pres(Local_freeDOF_HF(1:num_free_CalP))= 
     &                            x(num_FreeD+1:num_Total_FD)
          call Cal_Crack_Aperture(Counter_Iter,Try_DISP)
          call Cal_HF_delta_Time_Linear(Counter_Iter,
     &                    c_Temp_Cr_CalP_Aper,total_Time,delta_Time) 
          total_Time = c_Temp_total_Time + delta_Time
          call Cal_HF_Matrix_H_Linear(ifra,Counter_Iter,H,
     &                       Cracks_CalP_Aper,total_Time) 
          call Cal_HF_S(Counter_Iter,c_S,total_Time) 
          if(iter==1)then
              Last_DISP(1:c_Total_FD) = Initial_DISP
          else
              call Read_Disp_to(Counter_Iter-1,Last_DISP,
     &                          c_Total_FD)
          endif
          call Cal_HF_Resid(ifra,2,Counter_Iter,c_R_new,
     &                    c_Total_FD,num_FreeD,num_free_CalP,c_globalK,
     &                    Coupled_Q,F_U,Try_DISP,Last_DISP,
     &                    c_freeDOF,freeDOF_HF,Local_freeDOF_HF,
     &                    Try_CalP_Pres,delta_Time,H,c_S)  
          f_vector(1:num_Total_FD) =
     &                    c_R_new(Total_freeDOF(1:num_Total_FD))
          F_c = HLF*dot_product(f_vector,f_vector)          
          write(*,3002) F_c
          if(F_c < F_old + ALF*alam*slope) then
              print *,'          ++++++++++++++++++++++'
              print *,'           Line searching done! '
              print *,'          ++++++++++++++++++++++'
              goto 9999
          else
              print *,'              value of F is too big,'
     &                       //' goto the next backtracking'
          if(Lnsrch_Count==1) then
              tmplam = -slope/(TWO*(F_c-F_old-slope))
          else
              rhs1=F_c-F_old-alam*slope
              rhs2=F_2-F_old-alam2*slope
              a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
              b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
              if(a==ZR)then
                  tmplam = -slope/(TWO*b)
              else
                  disc=b*b-THR*a*slope
                  if(disc < ZR)then
                      tmplam =HLF*alam
                  elseif(b <= ZR)then
                      tmplam = (-b+sqrt(disc))/(THR*a)
                  else
                      tmplam = -slope/(b+sqrt(disc))
                  endif
              end if
              if(tmplam>HLF*alam)tmplam=HLF*alam
          end if
          end if
          alam2 = alam
          F_2 = F_c
          alam = max(tmplam,ZP1*alam)  
          print *,'              Value of λ:',alam
          goto 8888    
                      
 9999 continue 
      return 
      end SUBROUTINE Cal_HF_Line_Searching          
