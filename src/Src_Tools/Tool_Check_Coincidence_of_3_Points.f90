 
function Tool_Check_Coincidence_of_3_Points(p1, p2, p3, tol) result(coincident)

implicit none
real, dimension(3), intent(in) :: p1, p2, p3
real, intent(in) :: tol
logical :: coincident

coincident = all(abs(p1 - p2) <= tol) .and. all(abs(p1 - p3) <= tol)

end function Tool_Check_Coincidence_of_3_Points