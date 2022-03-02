
! just regular trapz
subroutine trapz(x, y, N_arr,r)
   !! Found this little routine online
   !! Calculates the integral of an array y with respect to x using the trapezoid
   !! approximation. Note that the mesh spacing of x does not have to be uniform.
   integer, intent(in)                         :: N_arr
   real(kind=8), dimension(N_arr), intent(in)  :: x      !! Variable x
   real(kind=8), dimension(N_arr), intent(in)  :: y      !! Function y(x)
   real(kind=8)                                :: r      !! Integral ∫y(x)·dx

   ! Integrate using the trapezoidal rule
   r = SUM((y(1+1:N_arr-0) + y(1+0:N_arr-1))*(x(1+1:N_arr-0) - x(1+0:N_arr-1)))/2
end subroutine
