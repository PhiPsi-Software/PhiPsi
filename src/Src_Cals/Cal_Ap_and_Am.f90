!-----------------------------------------------------------
! Brief: Determine the A+ and A- sub-domains of a Heaviside node.
!
! Parameters:
!   Input:  c_Crack, c_Cr_Point      - crack and crack-point count
!           DOMAIN_Outline, m_DOMAIN_Outline - support polygon
!           Domain_El, n_Domain_El   - elements in the support
!           c_E, c_Node              - element and node indices
!   Output: Yes_Cutthrough           - crack cuts the support
!           Yes_Ap_has_GP, Yes_Am_has_GP - sub-domains contain GPs
!
! Notes:   Used to decide whether a Heaviside-enriched node should
!          be removed: if one side has no Gauss points it is dropped.
!-----------------------------------------------------------

subroutine Cal_Ap_and_Am(c_Crack,c_Cr_Point, DOMAIN_Outline,m_DOMAIN_Outline, Domain_El,n_Domain_El,c_E,c_Node, &
Yes_Cutthrough, Yes_Ap_has_GP,Yes_Am_has_GP)
!     This function calculates the A+ and A- domain of the support domain of node, used for
!     the selection of the Heaviside node. If either A+ or A- does not contain any gauss point,
!     then the Heaviside node should be removed. Yes_Cutthrough equals to 1 if the support domain
!     is cut through by the crack. 

!     c_Cr_Point: Current number of crack coordinate points
use Global_Float_Type
use Global_Model
use Global_Elem_Area_Vol 
use Global_Common
use Global_Crack

implicit none     

integer,intent(in)::c_E,c_Node,c_Cr_Point
integer,intent(in)::m_DOMAIN_Outline,n_Domain_El
integer,intent(in)::DOMAIN_Outline(m_DOMAIN_Outline,2)
integer,intent(in)::Domain_El(n_Domain_El)
real(kind=FT),intent(in)::c_Crack(Max_Num_Cr_P,2)

logical,intent(out)::Yes_Cutthrough,Yes_Ap_has_GP,Yes_Am_has_GP

integer i_Seg,i_Side,i,j,k,i_E,i_G
integer tem_Num_Node,tem_count
integer :: Max_IVP
real(kind=FT),allocatable:: Inters_and_Vex_point(:,:)
real(kind=FT),allocatable:: Inters_point(:,:)
integer,allocatable:: new_DOMAIN_Outline(:,:)
integer :: m_new_DOMAIN_Outline
integer,allocatable:: Uniq_new_DOMAIN_Outline(:,:)
integer Num_In_and_Vex_Point,Num_Inters_point
integer,allocatable:: Cut_Segments(:,:)
integer :: m_Cut_Segments
integer,allocatable:: EndByEnd_Cut_Segments(:,:)
integer :: cou
integer,allocatable:: changed_DOMAIN_Side_Flag(:)
integer :: num_Cha_DOM_S_Flag
integer,allocatable:: Uniq_changed_DOMAIN_S_F(:)
real(kind=FT)  point_seg_1(2),point_seg_2(2)
real(kind=FT)  point_side_1(2),point_side_2(2),Ver_point(2)
real(kind=FT)  X,Y,Ap_or_Am
integer :: Uniq_num_Cha_DOM_S_Flag
logical Yes_Cross,Yes_point_seg_1_in,Yes_point_seg_2_in
integer :: EndByEnd_DOMAIN_Outline(m_DOMAIN_Outline+1)
real(kind=FT),allocatable:: tem_ALL_Coor(:,:)
logical :: Yes
integer,allocatable:: Uni_Mat_Count(:)
integer :: Uniq_m_new_DOMAIN_Outline
integer num_node_1,num_node_2,num_node_3,Location
integer PP(3),Uniqued_PP(3),n_Uniqued_PP
integer,allocatable:: Ap_Outline(:,:)
integer,allocatable:: Am_Outline(:,:)
integer,allocatable:: Sorted_Ap_Outline(:,:)
integer,allocatable:: Sorted_Am_Outline(:,:)
integer m_Ap_Outline,m_Am_Outline,c_point
real(kind=FT),allocatable:: Signed_Distance(:,:)
real(kind=FT) tol,c_point_Coor(2)
real(kind=FT) Line_AB(2,2)
integer :: cou_1
real(kind=FT),allocatable:: Coor_Ap_Closed(:,:), Coor_Am_Closed(:,:)
real(kind=FT) kesi(Num_Gauss_Points), yita(Num_Gauss_Points), weight(Num_Gauss_Points)
logical In_Ap,In_Am,c_Yes_ON,Yes_Special
integer :: NODES_iE(4)
real(kind=FT) X_NODES(4),Y_NODES(4)
integer c_point_1,c_point_2
real(kind=FT) c_point_Coor_1(2),c_point_Coor_2(2), c_Mid_Point(2),xpol(3),ypol(3)
real(kind=FT) :: c_Signed_Distance
real(kind=FT) :: VBA_area
real(kind=FT) :: temp_vector_1(2),temp_vector_2(2),temp_vector_3(2)

Max_IVP = (c_Cr_Point-1)*m_DOMAIN_Outline + 4
if (Max_IVP < 10) Max_IVP = 10

allocate(Inters_and_Vex_point(Max_IVP,2))
allocate(Inters_point(Max_IVP,2))
allocate(new_DOMAIN_Outline(2*Max_IVP + m_DOMAIN_Outline, 2))
allocate(Uniq_new_DOMAIN_Outline(2*Max_IVP + m_DOMAIN_Outline, 2))
allocate(Cut_Segments(Max_IVP,2))
allocate(EndByEnd_Cut_Segments(Max_IVP,2))
allocate(changed_DOMAIN_Side_Flag(Max_IVP))
allocate(Uniq_changed_DOMAIN_S_F(Max_IVP))
allocate(Uni_Mat_Count(2*Max_IVP + m_DOMAIN_Outline))
allocate(Ap_Outline(2*Max_IVP + m_DOMAIN_Outline, 2))
allocate(Am_Outline(2*Max_IVP + m_DOMAIN_Outline, 2))
allocate(Sorted_Ap_Outline(2*Max_IVP + m_DOMAIN_Outline, 2))
allocate(Sorted_Am_Outline(2*Max_IVP + m_DOMAIN_Outline, 2))
allocate(Signed_Distance(2, Max_IVP))
allocate(tem_ALL_Coor(Num_Node + 2*Max_IVP + 4, 2))


m_Ap_Outline         = 0
m_Am_Outline         = 0
Num_In_and_Vex_Point = 0
Num_Inters_point     = 0
m_Cut_Segments       = 0
num_Cha_DOM_S_Flag   = 0
tem_count            = 0
tem_Num_Node         = Num_Node
m_new_DOMAIN_Outline = 0
num_node_1           = 0
num_node_2           = 0
num_node_3           = 0
cou_1                = 0

In_Ap = .False.
In_Am = .False.

changed_DOMAIN_Side_Flag(1:Max_IVP) = 0
new_DOMAIN_Outline(1:2*Max_IVP + m_DOMAIN_Outline,1:2)   =0
Uniq_new_DOMAIN_Outline(1:2*Max_IVP + m_DOMAIN_Outline,1:2)=0
Cut_Segments(1:Max_IVP,1:2)  =0
EndByEnd_Cut_Segments(1:Max_IVP,1:2)=0
Signed_Distance(1:2,1:Max_IVP)  =ZR
Ap_Outline(1:2*Max_IVP + m_DOMAIN_Outline,1:2)   =0
Am_Outline(1:2*Max_IVP + m_DOMAIN_Outline,1:2)   =0
Sorted_Ap_Outline(1:2*Max_IVP + m_DOMAIN_Outline,1:2)   =0
Sorted_Am_Outline(1:2*Max_IVP + m_DOMAIN_Outline,1:2)   =0
Inters_and_Vex_point(1:Max_IVP,1:2)   =ZR

kesi(1:Num_Gauss_Points)=ZR
yita(1:Num_Gauss_Points)=ZR
weight(1:Num_Gauss_Points)=ZR


Yes_Cutthrough =.False.
Yes_Ap_has_GP = .False.
Yes_Am_has_GP = .False.

tem_ALL_Coor(1:Num_Node + 2*Max_IVP + 4,1:2) = ZR
tem_ALL_Coor(1:Num_Node,1:2) = Coor

do i_Seg = 1,c_Cr_Point-1
    point_seg_1 = c_Crack(i_Seg,:)
    point_seg_2 = c_Crack(i_Seg+1,:)
    do i_Side = 1,m_DOMAIN_Outline
        point_side_1 = [Coor(DOMAIN_Outline(i_Side,1),1), Coor(DOMAIN_Outline(i_Side,1),2)]
        point_side_2 = [Coor(DOMAIN_Outline(i_Side,2),1), Coor(DOMAIN_Outline(i_Side,2),2)]    
        call Tool_Intersection(point_seg_1,point_seg_2,point_side_1,point_side_2,X,Y,Yes_Cross) 

        if(Yes_Cross.eqv..True.)then

            tem_count = tem_count + 1
            Num_In_and_Vex_Point = Num_In_and_Vex_Point +1
            Inters_and_Vex_point(Num_In_and_Vex_Point,1:2)=[X,Y]
            Num_Inters_point = Num_Inters_point +1
            Inters_point(Num_Inters_point,:) = [X,Y]
            tem_Num_Node = tem_Num_Node+1
            tem_ALL_Coor(tem_Num_Node,:) = [X,Y]
            m_new_DOMAIN_Outline = m_new_DOMAIN_Outline +1
            new_DOMAIN_Outline(m_new_DOMAIN_Outline,:) =  [DOMAIN_Outline(i_Side,1),tem_Num_Node]
            m_new_DOMAIN_Outline = m_new_DOMAIN_Outline +1
            new_DOMAIN_Outline(m_new_DOMAIN_Outline,:) = [DOMAIN_Outline(i_Side,2),tem_Num_Node]      

            num_Cha_DOM_S_Flag = num_Cha_DOM_S_Flag +1 
            changed_DOMAIN_Side_Flag(num_Cha_DOM_S_Flag)=i_Side                 
        end if

    end do
    if (tem_count.ne.0) then
        Yes_point_seg_1_in = .False.
        Yes_point_seg_2_in = .False.
        if (i_Seg.ne. 1) then
            EndByEnd_DOMAIN_Outline(1:m_DOMAIN_Outline) = DOMAIN_Outline(:,1)
            EndByEnd_DOMAIN_Outline(m_DOMAIN_Outline+1) = DOMAIN_Outline(1,1)
            Call Tool_Yes_In_Poly(point_seg_1(1),point_seg_1(2),Coor(EndByEnd_DOMAIN_Outline,1), &
            Coor(EndByEnd_DOMAIN_Outline,2),size(EndByEnd_DOMAIN_Outline), Yes_point_seg_1_in)
        end if
        if (i_Seg .ne. (c_Cr_Point-1)) then
            EndByEnd_DOMAIN_Outline(1:m_DOMAIN_Outline) = DOMAIN_Outline(:,1)
            EndByEnd_DOMAIN_Outline(m_DOMAIN_Outline+1) = DOMAIN_Outline(1,1)
            Call Tool_Yes_In_Poly(point_seg_2(1),point_seg_2(2), Coor(EndByEnd_DOMAIN_Outline,1), Coor(EndByEnd_DOMAIN_Outline,2), &
            size(EndByEnd_DOMAIN_Outline), Yes_point_seg_2_in)
        end if
        if(Yes_point_seg_1_in)then  
            call Vector_belongs_Matrix_Is_Dou(Max_IVP,2, Inters_and_Vex_point, [point_seg_1(1),point_seg_1(2)], Location,Yes)
            if    (Yes.eqv..False.) then
                Num_In_and_Vex_Point = Num_In_and_Vex_Point +1
                Inters_and_Vex_point(Num_In_and_Vex_Point,1:2) = [point_seg_1(1),point_seg_1(2)]
                Ver_point = [point_seg_1(1),point_seg_1(2)]
            end if
        end if
        if(Yes_point_seg_2_in)then 
            call Vector_belongs_Matrix_Is_Dou(Max_IVP,2, Inters_and_Vex_point, [point_seg_2(1),point_seg_2(2)], Location,Yes)
            if(Yes .eqv..False.) then
                Num_In_and_Vex_Point = Num_In_and_Vex_Point +1
                Inters_and_Vex_point(Num_In_and_Vex_Point,1:2) = [point_seg_2(1),point_seg_2(2)]
                Ver_point = [point_seg_2(1),point_seg_2(2)]
            end if                  
        end if
    end if
end do


if ((Num_In_and_Vex_Point.eq.0)  .or. (Num_Inters_point  .eq. 1)  .or. (Num_Inters_point  .ge. 3)) then
    Yes_Cutthrough =.False.
    Yes_Ap_has_GP = .True.
    Yes_Am_has_GP = .True.
    return
end if

call Vector_Unique_Int(Max_IVP, num_Cha_DOM_S_Flag,changed_DOMAIN_Side_Flag, &
Uniq_changed_DOMAIN_S_F, Uniq_num_Cha_DOM_S_Flag)

do i_Side =1,m_DOMAIN_Outline
    if (any(i_Side .eq. changed_DOMAIN_Side_Flag(1:Uniq_num_Cha_DOM_S_Flag)).eqv..False.) then
        m_new_DOMAIN_Outline = m_new_DOMAIN_Outline +1
        new_DOMAIN_Outline(m_new_DOMAIN_Outline,:) = DOMAIN_Outline(i_Side,:)
    end if
end do

call Matrix_Unique_Row_Int(2*Max_IVP + m_DOMAIN_Outline, 2, m_new_DOMAIN_Outline, new_DOMAIN_Outline, &
Uniq_new_DOMAIN_Outline,Uniq_m_new_DOMAIN_Outline, Uni_Mat_Count)


select case(Num_In_and_Vex_Point)
case(2)
    temp_vector_1(1:2) =  Inters_and_Vex_point(1,1:2)
    temp_vector_2(1:2) =  Inters_and_Vex_point(2,1:2)
    call Vector_belongs_Matrix_Is_Dou(Num_Node + 2*Max_IVP + 4, 2, tem_ALL_Coor, temp_vector_1, num_node_1,Yes) 
    call Vector_belongs_Matrix_Is_Dou(Num_Node + 2*Max_IVP + 4, 2, tem_ALL_Coor, temp_vector_2, num_node_2,Yes)   
    m_Cut_Segments = m_Cut_Segments +1
    Cut_Segments(m_Cut_Segments,:) = [num_node_1,num_node_2]
case(3)
    tem_Num_Node = tem_Num_Node + 1
    tem_ALL_Coor(tem_Num_Node,:) = Ver_point
    temp_vector_1(1:2) =  Inters_and_Vex_point(1,1:2)
    temp_vector_2(1:2) =  Inters_and_Vex_point(2,1:2)
    temp_vector_3(1:2) =  Inters_and_Vex_point(3,1:2)

    call Vector_belongs_Matrix_Is_Dou(Num_Node + 2*Max_IVP + 4, 2, tem_ALL_Coor, temp_vector_1, num_node_1,Yes) 
    call Vector_belongs_Matrix_Is_Dou(Num_Node + 2*Max_IVP + 4, 2, tem_ALL_Coor, temp_vector_2, num_node_2,Yes) 
    call Vector_belongs_Matrix_Is_Dou(Num_Node + 2*Max_IVP + 4, 2, tem_ALL_Coor, temp_vector_3, num_node_3,Yes) 
    PP = [num_node_1,num_node_2,num_node_3]
    call Vector_Unique_Int(3,3,PP, Uniqued_PP,n_Uniqued_PP) 
    call Vector_Sort_Int(n_Uniqued_PP,Uniqued_PP(1:n_Uniqued_PP))  
    m_Cut_Segments = m_Cut_Segments +1
    Cut_Segments(m_Cut_Segments,1:2)=[Uniqued_PP(1),tem_Num_Node]
    m_Cut_Segments = m_Cut_Segments +1
    Cut_Segments(m_Cut_Segments,1:2)=[Uniqued_PP(2),tem_Num_Node]
case(4)
    Yes_Cutthrough= .True.
    Yes_Ap_has_GP =.True.
    Yes_Am_has_GP =.True.
    return
case default
    Yes_Cutthrough=.False.
    Yes_Ap_has_GP=.False.
    Yes_Am_has_GP=.False.
    return
end select


if(((Num_In_and_Vex_Point.eq.2) .or. (Num_In_and_Vex_Point.eq.3)) .and. (m_Cut_Segments.ne.0))     then
    Yes_Cutthrough = .True.
    call Tool_Sort_by_End_to_End(Max_IVP,m_Cut_Segments,Cut_Segments, EndByEnd_Cut_Segments,cou)

    tol = sqrt(Ave_Elem_Area)*1.0D-8

    if(Num_In_and_Vex_Point.eq.2)then
        do i=1,m_new_DOMAIN_Outline
            do j=1,2
                c_point = Uniq_new_DOMAIN_Outline(i,j)
                c_point_Coor = tem_ALL_Coor(c_point,1:2)
                do k=1,Num_Inters_point-1
                    point_seg_1 = Inters_point(k,1:2)
                    point_seg_2 = Inters_point(k+1,1:2)
                    Line_AB(1,1:2)= point_seg_1
                    Line_AB(2,1:2)= point_seg_2 
                    call Cal_Signed_Distance(Line_AB,c_point_Coor, Signed_Distance(j,k))
                    if (abs(Signed_Distance(j,k)) .le. tol) then
                        Signed_Distance(j,k) = ZR
                    end if
                end do
            end do
            Ap_or_Am = maxval(Signed_Distance(:,1:Num_Inters_point-1))
            if (Ap_or_Am .gt. tol) then
                m_Ap_Outline = m_Ap_Outline + 1
                Ap_Outline(m_Ap_Outline,:)= Uniq_new_DOMAIN_Outline(i,:)
            else
                m_Am_Outline = m_Am_Outline + 1
                Am_Outline(m_Am_Outline,:)=  Uniq_new_DOMAIN_Outline(i,:)
            end if
        enddo
    endif
    if(Num_In_and_Vex_Point.eq.3)then
        Yes_Special =.False.
        do i=1,m_new_DOMAIN_Outline
            c_point_1 = Uniq_new_DOMAIN_Outline(i,1)
            c_point_Coor_1 = tem_ALL_Coor(c_point_1,1:2)
            c_point_2 = Uniq_new_DOMAIN_Outline(i,2)
            c_point_Coor_2 = tem_ALL_Coor(c_point_2,1:2)
            c_Mid_Point(1:2)=HLF*(c_point_Coor_1 + c_point_Coor_2)
            Line_AB(1,1:2)= Inters_point(1,1:2)
            Line_AB(2,1:2)= Inters_point(2,1:2)
            call Tool_Yes_On_Line(c_Mid_Point(1),c_Mid_Point(2),Inters_point(1,1:2),Inters_point(2,1:2),c_Yes_ON)
            if(c_Yes_ON)then
                Yes_Special =.True.
                exit
            endif
        end do
        if(Yes_Special)then
            xpol(1) = tem_ALL_Coor(PP(1),1)
            xpol(2) = tem_ALL_Coor(PP(2),1)
            xpol(3) = tem_ALL_Coor(PP(3),1)
            ypol(1) = tem_ALL_Coor(PP(1),2)
            ypol(2) = tem_ALL_Coor(PP(2),2)
            ypol(3) = tem_ALL_Coor(PP(3),2)
            call Tool_Area_Polygon(xpol,ypol,3,VBA_area)
            if(abs(VBA_area) >= 1.0D-3*Ave_Elem_Area)then
                Yes_Ap_has_GP =.True.
                Yes_Am_has_GP =.True.
                return
            endif
        else
            do i=1,m_new_DOMAIN_Outline
                c_point_1 = Uniq_new_DOMAIN_Outline(i,1)
                c_point_Coor_1 = tem_ALL_Coor(c_point_1,1:2)
                c_point_2 = Uniq_new_DOMAIN_Outline(i,2)
                c_point_Coor_2 = tem_ALL_Coor(c_point_2,1:2)
                c_Mid_Point(1:2)=HLF*(c_point_Coor_1 + c_point_Coor_2)
                Line_AB(1,1:2)= Inters_point(1,1:2)
                Line_AB(2,1:2)= Inters_point(2,1:2)
                call Cal_Signed_Distance(Line_AB,c_Mid_Point(1:2), c_Signed_Distance)

                if (c_Signed_Distance > ZR) then
                    m_Ap_Outline = m_Ap_Outline + 1
                    Ap_Outline(m_Ap_Outline,:)= Uniq_new_DOMAIN_Outline(i,1:2)
                elseif(c_Signed_Distance < ZR) then
                    m_Am_Outline = m_Am_Outline + 1
                    Am_Outline(m_Am_Outline,:)= Uniq_new_DOMAIN_Outline(i,1:2)
                end if
            enddo
        endif
    endif
    do i=1,m_Cut_Segments
        m_Ap_Outline = m_Ap_Outline + 1
        Ap_Outline(m_Ap_Outline,:) = EndByEnd_Cut_Segments(i,:)
        m_Am_Outline = m_Am_Outline + 1
        Am_Outline(m_Am_Outline,:) = EndByEnd_Cut_Segments(i,:)
    end do  



    call Tool_Sort_by_End_to_End(2*Max_IVP + m_DOMAIN_Outline, m_Ap_Outline, Ap_Outline, Sorted_Ap_Outline, cou_1)
    call Tool_Sort_by_End_to_End(2*Max_IVP + m_DOMAIN_Outline, m_Am_Outline, Am_Outline, Sorted_Am_Outline, cou_1)
    allocate(Coor_Ap_Closed(m_Ap_Outline+1,2))
    allocate(Coor_Am_Closed(m_Am_Outline+1,2))
    Coor_Ap_Closed(1:m_Ap_Outline+1,1:2)=ZR
    Coor_Am_Closed(1:m_Am_Outline+1,1:2)=ZR
    Coor_Ap_Closed(1:m_Ap_Outline,:) = tem_ALL_Coor(Sorted_Ap_Outline(1:m_Ap_Outline,1),:)
    Coor_Ap_Closed(m_Ap_Outline+1,:) = tem_ALL_Coor(Sorted_Ap_Outline(1,1),:)

    Coor_Am_Closed(1:m_Am_Outline,:) = tem_ALL_Coor(Sorted_Am_Outline(1:m_Am_Outline,1),:)   
    Coor_Am_Closed(m_Am_Outline+1,:) = tem_ALL_Coor(Sorted_Am_Outline(1,1),:) 

    if (Key_Integral_Sol  == 2)then
        call Cal_Gauss_Points_QUAD(Num_Gauss_Points, kesi,yita,weight)
    elseif (Key_Integral_Sol  == 3)then
        call Cal_Gauss_Points_QUAD_for_SUBQUAD(Num_Sub_Quads, kesi,yita,weight)
        Num_Gauss_Points = Num_Sub_Quads*4
    endif

    do i_E = 1,n_Domain_El
        NODES_iE = [Elem_Node(Domain_El(i_E),1), Elem_Node(Domain_El(i_E),2), &
        Elem_Node(Domain_El(i_E),3), Elem_Node(Domain_El(i_E),4)]
        X_NODES = Coor(NODES_iE,1)                           
        Y_NODES = Coor(NODES_iE,2)
        do i_G = 1,Num_Gauss_Points
            call Cal_Coor_by_KesiYita(kesi(i_G),yita(i_G), X_NODES,Y_NODES, x,y)

            call Tool_Yes_In_Poly(x,y, Coor_Ap_Closed(:,1), Coor_Ap_Closed(:,2), m_Ap_Outline+1,In_Ap)
            if (In_Ap) then
                Yes_Ap_has_GP = .True.
                goto 201
            end if
        end do
    end do        
    201     continue
    do i_E = 1,n_Domain_El
        NODES_iE = [Elem_Node(Domain_El(i_E),1), Elem_Node(Domain_El(i_E),2), Elem_Node(Domain_El(i_E),3), &
        Elem_Node(Domain_El(i_E),4)]
        X_NODES = Coor(NODES_iE,1)                           
        Y_NODES = Coor(NODES_iE,2)
        do i_G = 1,Num_Gauss_Points
            call Cal_Coor_by_KesiYita(kesi(i_G),yita(i_G), X_NODES,Y_NODES, x,y)

            call Tool_Yes_In_Poly(x,y, Coor_Am_Closed(:,1), Coor_Am_Closed(:,2), m_Am_Outline+1,In_Am)
            if (In_Am) then
                Yes_Am_has_GP = .True.
                goto 202
            end if
        end do
    end do        
    202     continue          
end if

return
end subroutine Cal_Ap_and_Am                         
