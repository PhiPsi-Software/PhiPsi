 
function Tool_Check_Coincidence_of_2_Points(p1, p2, tol) result(coincident)

implicit none
real, dimension(3), intent(in) :: p1, p2
real, intent(in) :: tol
logical :: coincident

coincident = all(abs(p1 - p2) <= tol)

end function Tool_Check_Coincidence_of_2_Points