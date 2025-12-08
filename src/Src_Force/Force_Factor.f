 
      SUBROUTINE Force_Factor(Lambda,isub,Yes_Last_Growth)

      use Global_Float_Type
      use Global_Common
      use Global_Dynamic
      
      implicit none
      logical Yes_Last_Growth
      integer isub
      real(kind=FT) Lambda
      
      select case(Key_Force_Control)
      case(1)
            if (Key_Analysis_Type==1 .or.
     &        Key_Analysis_Type==3 .or. Key_Analysis_Type==4)then
              Lambda = ONE
          elseif (Key_Analysis_Type.eq.2) then
              if (isub .le. IDy_Num_force_Itr) then
                  Lambda = ONE
            else 
            Lambda = ZR
            end if
          elseif (Key_Analysis_Type.eq.6) then
              if (isub .le. EDy_Num_force_Itr) then
                  Lambda = ONE
            else 
            Lambda = ZR
            end if
          end if
      case(2)
            if (Key_Analysis_Type==1 .or.
     &        Key_Analysis_Type==3 .or. Key_Analysis_Type==4)then
              Lambda = Lambda + ONE/Num_Substeps
          elseif (Key_Analysis_Type.eq.2) then
              if (isub .le. IDy_Num_force_Itr) then
                  Lambda = ONE
            else 
            Lambda = ZR
            end if
          elseif (Key_Analysis_Type.eq.6) then
              if (isub .le. EDy_Num_force_Itr) then
                  Lambda = ONE
            else 
            Lambda = ZR
            end if
          end if
      case(3)    
            if (Key_Analysis_Type==1 .or.
     &        Key_Analysis_Type==3 .or. Key_Analysis_Type==4)then
              if (Yes_Last_Growth .eqv. .False.) then
                  Lambda = Lambda + ONE/Num_Force_Divs
              else
            Lambda = 1/Num_Force_Divs
               end if
          elseif (Key_Analysis_Type.eq.2) then
              if (isub .le. IDy_Num_force_Itr) then
                  Lambda = ONE
            else 
            Lambda = ZR
            end if
          elseif (Key_Analysis_Type.eq.6) then
              if (isub .le. EDy_Num_force_Itr) then
                  Lambda = ONE
            else 
            Lambda = ZR
            end if
          end if
      case(4)
          Lambda = ONE
      end select
      
      RETURN
      END SUBROUTINE Force_Factor
