
! just regular trapz
pure function trapz(x, y) result(r)
   !! Found this little routine online
   !! Calculates the integral of an array y with respect to x using the trapezoid
   !! approximation. Note that the mesh spacing of x does not have to be uniform.
   real(kind=8), intent(in)  :: x(:)         !! Variable x
   real(kind=8), intent(in)  :: y(size(x))   !! Function y(x)
   real(kind=8)              :: r            !! Integral ∫y(x)·dx

   ! Integrate using the trapezoidal rule
   associate(n => size(x))
     r = SUM((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
   end associate
end function