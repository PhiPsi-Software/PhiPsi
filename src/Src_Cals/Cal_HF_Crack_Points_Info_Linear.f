 
      subroutine Cal_HF_Crack_Points_Info_Linear(isub)
          
      
      
      use Global_Float_Type                                                                
      use Global_Crack
      use Global_Crack_Common
      use Global_Model
      use Global_Elem_Area_Vol
      use Global_Common
      use Global_HF
      use Global_Stress
      use module_INTERFACE_Vector_belongs_Matrix_Is_Dou

      implicit none
      integer,intent(in)::isub
      integer i_C,i_S,i_P
      real(kind=FT) Cr_p1(2),Cr_p2(2),L_Seg
      integer Num_Div_Points
      real(kind=FT) Div_Points(Max_Num_Cr_CalP,2)   
      integer Div_Points_Seg_Num(Max_Num_Cr_CalP)       
      integer Num_S_Div_Points
      real(kind=FT) S_Div_Points(Max_Num_Seg_CalP,2) 
      real(kind=FT)  x1,y1,x2,y2,HF_Length
      logical Yes_Tip1_Edge_Crack,Yes_Tip2_Edge_Crack
      real(kind=FT) Line_AB(2,2),delta_L,
     &                 new_Point(2),new_Line_AB(2,2)
      integer i_HF_Elem
      integer i_Jun,n_Crs,num_Jun,i_new_Seg
      real(kind=FT) x_Jun,y_Jun
      logical Yes_SegHasJun,Yes_Jun_OnSeg
      real(kind=FT) New_Seg_withJun(Max_Num_Cone_Cr+1,2)
      real(kind=FT) New_Seg_P1(2),New_Seg_P2(2) 
      real(kind=FT) New_S_Div_Points(Max_Num_Seg_CalP,2) 
      integer New_Num_S_Div_Points
      integer Total_num_Jun,c_Jun_Num,counter_CalP,i_CalP,Other_Crack
      integer j_CalP,i_Pair
      real(kind=FT) Tiny_L
      integer j_Jun,Proba_Ele,num_Jun_Crack
      integer c_Crack,Ele,c_CalP
      real(kind=FT) Inj_Poin_x,Inj_Poin_y,Min_Distance,c_Distance,
     &                 c_CalP_x,c_CalP_y
      integer i_tt
      logical Yes_Has
      integer Location
      real(kind=FT) c_D_P_x,c_D_P_y
      real(kind=FT) c_theta,Orient_L,Orient_R
      integer i_Cr_P,i_Coor
      real(kind=FT) c_Cr_Point(2),c_Div_P(2),c_Tiny
      real(kind=FT) Old_Div_Points(Max_Num_Cr_CalP,2)   
      integer Old_Num_Div_Points
      logical Yes_Overlap
      real(kind=FT) c_InSitu_x,c_InSitu_y
      integer c_Ele
      real(kind=FT) c_angle,Contact_stress_Normal,Contact_stress_Tan
      
      print *,'    Calculating points of each crack...'
      
      Cracks_CalP_Coors(1:Max_Num_Cr,1:Max_Num_Cr_CalP,1:2) = ZR
      Cracks_CalP_Orient(1:Max_Num_Cr,1:Max_Num_Cr_CalP)    = ZR
      Cracks_CalP_Seg(1:Max_Num_Cr,1:Max_Num_Cr_CalP)       = 0
      Cracks_CalP_Elem(1:Max_Num_Cr,1:Max_Num_Cr_CalP)      = 0
      Cracks_HF_Ele_L(1:Max_Num_Cr,1:Max_Num_Cr_CalP-1)     = 0
      Cracks_CalP_Num(1:Max_Num_Cr)                         = 0
      Cracks_CalP_Type(1:Max_Num_Cr,1:Max_Num_Cr_CalP,1:2)  = 0
      Cracks_TipJunCalpNum(1:Max_Num_Cr,1:2)                = 0 
      Cracks_MidJunCalpNum(1:Max_Num_Cr,1:Max_Num_Cr_CalP)  = 0
      Cracks_GloNumCalP(1:Max_Num_Cr,1:Max_Num_Cr_CalP)     = 0
      Cracks_LocalNumCalP(1:Max_Num_Cr*Max_Num_Cr_CalP,1:2) = 0
      Div_Points(1:Max_Num_Cr_CalP,1:2)                     = ZR
      Cracks_CalP_Contact_Strs(1:Max_Num_Cr,1:Max_Num_Cr_CalP,1:2) = ZR
        Cracks_CalP_Contact_Force_x(1:Max_Num_Cr,1:Max_Num_Cr_CalP)=ZR
      Cracks_CalP_Contact_Force_y(1:Max_Num_Cr,1:Max_Num_Cr_CalP)=ZR
      Num_JunPair = 0
      Cracks_JunPair(1:Max_Num_Cr,1:2) = 0
      
      delta_L = Delta_Factor_Junc*Ave_Elem_L_Enrich
      
      do i_C = 1,num_Crack
          Yes_Tip1_Edge_Crack = .False.
          Yes_Tip2_Edge_Crack = .False.
          Num_Div_Points = 0
          Total_num_Jun  = 0
          
          do i_S = 1,Each_Cr_Poi_Num(i_C)-1
              New_Seg_withJun(1:Max_Num_Cone_Cr-2+1,1:2)= ZR
              Cr_p1 = [Edge_Disposed_Crack(i_C,i_S,1),
     &                 Edge_Disposed_Crack(i_C,i_S,2)]
              Cr_p2 = [Edge_Disposed_Crack(i_C,i_S+1,1),
     &                 Edge_Disposed_Crack(i_C,i_S+1,2)]

              L_Seg=sqrt((Cr_p2(1)-Cr_p1(1))**2+(Cr_p2(2)-Cr_p1(2))**2) 
              Yes_SegHasJun = .False.
              n_Crs = Cracks_Cone_NumMidCr(i_C)
              New_Seg_withJun(1,1:2) = Cr_p1
              num_Jun = 0
              do i_Jun=1,n_Crs
                  x_Jun = Cracks_Cone_MidJuCor(i_C,i_Jun,1)
                  y_Jun = Cracks_Cone_MidJuCor(i_C,i_Jun,2)
                  call Tool_Yes_On_Line(x_Jun,y_Jun,Cr_p1,Cr_p2,
     &                                  Yes_Jun_OnSeg)
                  if (Yes_Jun_OnSeg.eqv. .True.) then
                      Yes_SegHasJun = .True.
                      num_Jun =  num_Jun +1
                      New_Seg_withJun(num_Jun+1,1:2) = [x_Jun,y_Jun] 
                  endif
              end do
              New_Seg_withJun(num_Jun+2,1:2) = Cr_p2
              
              if (num_Jun>=2) then
                  call Tool_Rank_LineSeg_Points(
     &                      New_Seg_withJun(1:num_Jun+2,1:2),num_Jun+2) 
              end if
              if (Yes_SegHasJun .eqv. .False.) then          
                  call Cal_HF_HFNode_Location(Cr_p1,Cr_p2,
     &                                    S_Div_Points,Num_S_Div_Points)
                  if (i_S  < Each_Cr_Poi_Num(i_C)-1) then
                      Num_S_Div_Points = Num_S_Div_Points-1
                  end if
                  do i_tt = 1,Num_S_Div_Points
                      call Vector_belongs_Matrix_Is_Dou(
     &                               Num_Div_Points,2,
     &                               Div_Points(1:Num_Div_Points,1:2),
     &                               S_Div_Points(i_tt,1:2),
     &                               Location,Yes_Has)
                      if(Yes_Has.eqv..False.)then
                          Num_Div_Points = Num_Div_Points+1
                          Div_Points(Num_Div_Points,1:2)
     &                                   = S_Div_Points(i_tt,1:2)
                          Div_Points_Seg_Num(Num_Div_Points)=i_S
                      endif
                  enddo
              elseif(Yes_SegHasJun .eqv. .True.) then
                  do i_new_Seg =1,num_Jun+1
                      New_Seg_P1(1:2)=New_Seg_withJun(i_new_Seg,1:2)
                      New_Seg_P2(1:2)=New_Seg_withJun(i_new_Seg+1,1:2)
                      if (i_new_Seg==1)then     
                          Num_Div_Points  = Num_Div_Points +1
                          Div_Points(Num_Div_Points,1:2) = 
     &                                         New_Seg_withJun(1,1:2)
                          Div_Points_Seg_Num(Num_Div_Points) = i_S
                      end if
                      call Cal_HF_HFNode_Location(New_Seg_P1,New_Seg_P2,
     &                           New_S_Div_Points,New_Num_S_Div_Points)
                      
                      do i_tt = 1,New_Num_S_Div_Points
                          call Vector_belongs_Matrix_Is_Dou(
     &                               Num_Div_Points,2,
     &                               Div_Points(1:Num_Div_Points,1:2),
     &                               New_S_Div_Points(i_tt,1:2),
     &                               Location,Yes_Has)
                          if(Yes_Has.eqv..False.)then
                              Num_Div_Points = Num_Div_Points+1
                              Div_Points(Num_Div_Points,1:2)
     &                                   = New_S_Div_Points(i_tt,1:2)
                              Div_Points_Seg_Num(Num_Div_Points)=i_S
                              x_Jun = New_Seg_withJun(i_new_Seg+1,1)
                              y_Jun = New_Seg_withJun(i_new_Seg+1,2)
                              c_D_P_x = Div_Points(Num_Div_Points,1)
                              c_D_P_y = Div_Points(Num_Div_Points,2)
                              if(x_Jun==c_D_P_x.and.y_Jun==c_D_P_y)then
                                c_Jun_Num = Total_num_Jun+i_new_Seg
                                Cracks_MidJunCalpNum(i_C,c_Jun_Num)
     &                                         = Num_Div_Points
                                Cracks_CalP_Type(i_C,Num_Div_Points,1)=3   
                                Tiny_L = 1.0D-9*sqrt(Ave_Elem_L_Enrich)
                                num_Jun_Crack=Cracks_Cone_NumMidCr(i_C)
                                do j_Jun = 1,num_Jun_Crack
                                  Proba_Ele=Cracks_Cone_MidJuEle(i_C,
     &                                                           j_Jun) 
                                  if((abs(Cracks_Cone_MidJuCor(i_C,
     &                                 j_Jun,1)- x_Jun)<=Tiny_L) .and.
     &                              (abs(Cracks_Cone_MidJuCor(i_C,
     &                               j_Jun,2)- y_Jun)<=Tiny_L)) then
                                  Cracks_CalP_Type(i_C,Num_Div_Points,2)
     &                              = Cracks_Cone_MidCrNum(i_C,j_Jun)
                                   exit
                                 end if
                                end do         
                              end if
                          end if
                      end do 
                  end do
              endif
              Total_num_Jun = Total_num_Jun + num_Jun
          end do
          
          
          if (Crack_Tip_Type(i_C,1).eq.1) then 
              Line_AB(1,:) =[Edge_Disposed_Crack(i_C,1,1),
     &                       Edge_Disposed_Crack(i_C,1,2)]
              Line_AB(2,:) =[Edge_Disposed_Crack(i_C,2,1),
     &                       Edge_Disposed_Crack(i_C,2,2)]
              call Tool_Shorten_or_Extend_Line(Line_AB,-delta_L,'A',
     &                                      new_Line_AB,new_Point)
              Div_Points(1,:) = new_Point
          end if
          if (Crack_Tip_Type(i_C,2).eq.1) then   
              Line_AB(1,:) =[Edge_Disposed_Crack
     &                          (i_C,Each_Cr_Poi_Num(i_C)-1,1),
     &                       Edge_Disposed_Crack
     &                          (i_C,Each_Cr_Poi_Num(i_C)-1,2)]
              Line_AB(2,:) =[Edge_Disposed_Crack
     &                          (i_C,Each_Cr_Poi_Num(i_C),1),
     &                       Edge_Disposed_Crack
     &                          (i_C,Each_Cr_Poi_Num(i_C),2)]
              call Tool_Shorten_or_Extend_Line(Line_AB,-delta_L,'B',
     &                                      new_Line_AB,new_Point)
              Div_Points(Num_Div_Points,:) = new_Point 
          end if
          
          if (Key_Kink_Point==0)then
            Old_Div_Points(1:Max_Num_Cr_CalP,1:2)=
     &                        Div_Points(1:Max_Num_Cr_CalP,1:2)
            Old_Num_Div_Points = Num_Div_Points
            c_Tiny = Tol_11*Ave_Elem_L_Enrich   
            Num_Div_Points = 0
            do i_Cr_P = 1,Old_Num_Div_Points
              Yes_Overlap = .False.
              c_Div_P = Old_Div_Points(i_Cr_P,1:2) 
              if(i_Cr_P==1 .or. i_Cr_P == Old_Num_Div_Points)then
                  Num_Div_Points = Num_Div_Points + 1
                  Div_Points(Num_Div_Points,1:2) = c_Div_P
              else
                  if(Each_Cr_Poi_Num(i_C)>=3)then
                      do i_Coor = 2,Each_Cr_Poi_Num(i_C)-1
                          c_Cr_Point = Crack_Coor(i_C,i_Coor,1:2)
                          if((abs(c_Cr_Point(1)-c_Div_P(1))<c_Tiny).and.
     &                       (abs(c_Cr_Point(2)-c_Div_P(2))<c_Tiny))
     &                                                       then
                              Yes_Overlap = .True.
                              exit
                          endif
                      enddo
                      if(Yes_Overlap.eqv..False.)then
                          Num_Div_Points = Num_Div_Points + 1
                              Div_Points(Num_Div_Points,1:2) = c_Div_P 
                      endif
                  else
                      Num_Div_Points = Num_Div_Points + 1
                      Div_Points(Num_Div_Points,1:2) = c_Div_P
                  endif
              endif
            enddo
          endif
          if (Crack_Tip_Type(i_C,1).ne.-2) then
              Cracks_CalP_Type(i_C,1,1) = 5
          endif
          if (Crack_Tip_Type(i_C,2).ne.-2) then
              Cracks_CalP_Type(i_C,Num_Div_Points,1) = 5
          endif
          if (Crack_Tip_Type(i_C,1).eq.1) then
              if (Cracks_Cone_TipCrNum(i_C,1).ne.0) then
                  Cracks_CalP_Type(i_C,1,1) = 1
                  Cracks_CalP_Type(i_C,1,2)=Cracks_Cone_TipCrNum(i_C,1)
                  Cracks_TipJunCalpNum(i_C,1) = 1
              end if
          end if
          if (Crack_Tip_Type(i_C,2).eq.1) then
              if (Cracks_Cone_TipCrNum(i_C,2).ne.0) then
                  Cracks_CalP_Type(i_C,Num_Div_Points,1) = 2
                  Cracks_CalP_Type(i_C,Num_Div_Points,2) = 
     &                            Cracks_Cone_TipCrNum(i_C,2)
                  Cracks_TipJunCalpNum(i_C,2) = Num_Div_Points
              end if
          end if
          
          Cracks_CalP_Coors(i_C,1:Num_Div_Points,1:2) = 
     &                     Div_Points(1:Num_Div_Points,1:2)
          Cracks_CalP_Seg(i_C,1:Num_Div_Points)       = 
     &                           Div_Points_Seg_Num(1:Num_Div_Points)
          Cracks_CalP_Num(i_C)  = Num_Div_Points
          
          do i_P = 1,Num_Div_Points
              call Cal_Ele_Num_by_Coors(Div_Points(i_P,1),
     &                                  Div_Points(i_P,2),
     &                                  Cracks_CalP_Elem(i_C,i_P))
              if(i_P==1)then
                  Cr_p1 = [Cracks_CalP_Coors(i_C,1,1),
     &                     Cracks_CalP_Coors(i_C,1,2)]
                  Cr_p2 = [Cracks_CalP_Coors(i_C,2,1),
     &                     Cracks_CalP_Coors(i_C,2,2)]
                  Cracks_CalP_Orient(i_C,i_P) = 
     &                 atan2(Cr_p2(2)-Cr_p1(2),Cr_p2(1)-Cr_p1(1))
              elseif(i_P==Num_Div_Points)then
                  Cr_p1 = [Cracks_CalP_Coors(i_C,Num_Div_Points-1,1),
     &                     Cracks_CalP_Coors(i_C,Num_Div_Points-1,2)]
                  Cr_p2 = [Cracks_CalP_Coors(i_C,Num_Div_Points,1),
     &                     Cracks_CalP_Coors(i_C,Num_Div_Points,2)]
                  Cracks_CalP_Orient(i_C,i_P) = 
     &                 atan2(Cr_p2(2)-Cr_p1(2),Cr_p2(1)-Cr_p1(1))
              else
                  Cr_p1 = [Cracks_CalP_Coors(i_C,i_P-1,1),
     &                     Cracks_CalP_Coors(i_C,i_P-1,2)]
                  Cr_p2 = [Cracks_CalP_Coors(i_C,i_P,1),
     &                     Cracks_CalP_Coors(i_C,i_P,2)]
                  Orient_L= atan2(Cr_p2(2)-Cr_p1(2),Cr_p2(1)-Cr_p1(1))
                  Cr_p1 = [Cracks_CalP_Coors(i_C,i_P,1),
     &                     Cracks_CalP_Coors(i_C,i_P,2)]
                  Cr_p2 = [Cracks_CalP_Coors(i_C,i_P+1,1),
     &                     Cracks_CalP_Coors(i_C,i_P+1,2)]
                  Orient_R= atan2(Cr_p2(2)-Cr_p1(2),Cr_p2(1)-Cr_p1(1))
                  Cracks_CalP_Orient(i_C,i_P)=(Orient_L+Orient_R)*HLF
              endif
          end do
          do i_HF_Elem = 1,Num_Div_Points-1
            x1= Cracks_CalP_Coors(i_C,i_HF_Elem,1)
            y1= Cracks_CalP_Coors(i_C,i_HF_Elem,2)
            x2= Cracks_CalP_Coors(i_C,i_HF_Elem+1,1)
            y2= Cracks_CalP_Coors(i_C,i_HF_Elem+1,2)
              HF_Length = sqrt((x2-x1)**2 + (y2-y1)**2)
              Cracks_HF_Ele_L(i_C,i_HF_Elem) = HF_Length
          end do
          
      end do

      counter_CalP = 0
      do i_C = 1,num_Crack
          do i_CalP = 1,Cracks_CalP_Num(i_C) 
              counter_CalP = counter_CalP + 1
              Cracks_GloNumCalP(i_C,i_CalP) = counter_CalP
              Cracks_LocalNumCalP(counter_CalP,1) = i_C
              Cracks_LocalNumCalP(counter_CalP,2) = i_CalP
          end do
      end do
      num_Tol_CalP_All = counter_CalP
      
      counter_CalP = 0
      do i_C = 1,num_Crack
          if (Cracks_HF_State(i_C) == 1) then  
              do i_CalP = 1,Cracks_CalP_Num(i_C) 
                  counter_CalP = counter_CalP + 1
                  Cracks_GloNumCalP_W(i_C,i_CalP) = counter_CalP
                  Cracks_LocalNumCalP_W(counter_CalP,1) = i_C
                  Cracks_LocalNumCalP_W(counter_CalP,2) = i_CalP
              end do
          end if
      end do
      num_Tol_CalP_Water = counter_CalP
      
      if (Key_Contact /= 0) then
          Ele_NumCalP(1:Num_Elem) = 0
          Ele_CalPNum(1:Num_Elem,1:Max_Num_Ele_CalP) = 0
          do i_CalP = 1,num_Tol_CalP_All
              c_Crack = Cracks_LocalNumCalP(i_CalP,1)
              c_CalP = Cracks_LocalNumCalP(i_CalP,2)
              Ele = Cracks_CalP_Elem(c_Crack,c_CalP)
              Ele_NumCalP(Ele) = Ele_NumCalP(Ele) +1
              Ele_CalPNum(Ele,Ele_NumCalP(Ele)) = i_CalP
          enddo 
      end if
      
      do i_C = 1,num_Crack
          do i_CalP = 1,Cracks_CalP_Num(i_C) 
              if ((Cracks_CalP_Type(i_C,i_CalP,1) == 1) .or.
     &            (Cracks_CalP_Type(i_C,i_CalP,1) == 2)) then
                  Other_Crack = Cracks_CalP_Type(i_C,i_CalP,2)
                  if((Cracks_HF_State(i_C) == 1) .and.
     &               (Cracks_HF_State(Other_Crack) == 1) )then
                   do j_CalP=1,Cracks_CalP_Num(Other_Crack)
                     if((Cracks_CalP_Type(Other_Crack,j_CalP,1)==3).and.
     &                  (Cracks_CalP_Type(Other_Crack,j_CalP,2)==i_C))
     &                                                             then

                      Num_JunPair = Num_JunPair +1
                      Cracks_JunPair(Num_JunPair,1) =
     &                         Cracks_GloNumCalP(i_C,i_CalP)
                      Cracks_JunPair(Num_JunPair,2) =
     &                         Cracks_GloNumCalP(Other_Crack,j_CalP)
                    endif
                   end do
                  endif
              end if
          end do
      end do
      
      if (Num_JunPair>0) then
          print *,'         >> Junction pairs of CalP were found:'
      endif
 5001 FORMAT(15X,'>> Pair ',I3,' ---',I5,' and',I5)      
      
      do i_Pair =1,Num_JunPair 
          write(*,5001)i_Pair,Cracks_JunPair(i_Pair,1),
     &                        Cracks_JunPair(i_Pair,2)
      end do   
      if(Key_Analysis_Type == 3 .or. Key_Analysis_Type == 4)then
          Inj_Poin_x = Inj_Point_Loc(1)
          Inj_Poin_y = Inj_Point_Loc(2)
          do i_C=1,num_Crack
              if (i_C == Inject_Crack_Num) then
                   Min_Distance = Max_X_Coor-Min_X_Coor
                   do i_CalP = 1,Cracks_CalP_Num(i_C)
                       c_CalP_x = Cracks_CalP_Coors(i_C,i_CalP,1) 
                       c_CalP_y = Cracks_CalP_Coors(i_C,i_CalP,2) 
                       c_Distance = sqrt((Inj_Poin_x-c_CalP_x)**2 +
     &                                   (Inj_Poin_y-c_CalP_y)**2)
                       if(c_Distance < Min_Distance) then
                           Min_Distance = c_Distance
                           CalP_num_InjP_Local = i_CalP
                       endif
                   end do
              end if
          end do
      end if
      
      do i_C=1,num_Crack
          do i_CalP = 1,Cracks_CalP_Num(i_C)
              c_theta = Cracks_CalP_Orient(i_C,i_CalP) 
              Cracks_CalP_Remo_Strs(i_C,i_CalP) = 
     &              InSitu_y*cos(c_theta)**2 + InSitu_x*sin(c_theta)**2
          end do
      end do 
      
      if (Key_InSitu_Strategy == 2) then
          do i_C=1,num_Crack
              do i_CalP = 1,Cracks_CalP_Num(i_C)
                  c_Ele = Cracks_CalP_Elem(i_C,i_CalP) 
                  c_InSitu_x = 
     &               sum(InSitu_Strs_Gaus_xx(c_Ele,1:Num_Gauss_P_FEM))/
     &                                    Num_Gauss_P_FEM
                  c_InSitu_y = 
     &               sum(InSitu_Strs_Gaus_yy(c_Ele,1:Num_Gauss_P_FEM))/
     &                                    Num_Gauss_P_FEM    
                  c_theta = Cracks_CalP_Orient(i_C,i_CalP) 
     
                  Contact_stress_Normal = c_InSitu_x*sin(c_theta)**2 +
     &                                    c_InSitu_y*cos(c_theta)**2 
                  Contact_stress_Tan = 
     &                (c_InSitu_x-c_InSitu_y)*sin(c_theta)*cos(c_theta)
     
                  Cracks_CalP_Contact_Strs(i_C,i_CalP,1) = 
     &                                            Contact_stress_Normal
                  Cracks_CalP_Contact_Strs(i_C,i_CalP,2) = 
     &                                            Contact_stress_Tan
              end do
          end do       
      endif
      return 
      end SUBROUTINE Cal_HF_Crack_Points_Info_Linear
