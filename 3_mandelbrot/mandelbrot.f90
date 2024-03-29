program mandelbrot
   use omp_lib
  implicit none
  integer, parameter :: sp = kind(1.0)
  integer, parameter :: dp = kind(1.0d0)
  integer :: i, j, iter, numoutside
  integer, parameter :: npoints = 2000, maxiter = 2000
  real (kind=dp) :: area, error, start, finish
  complex (kind=dp) :: c , z

  start = omp_get_wtime()
  numoutside = 0

  !$omp parallel do default(none) private(i, j, c, z, iter) reduction(+:numoutside)
  do i = 0, npoints-1
     do j= 0,npoints-1
        c = cmplx(-2.0+(2.5*i)/npoints + 1.0d-07,(1.125*j)/npoints + 1.0d-07)
        z = c
        iter = 0
        do while (iter < maxiter)
           z = z*z + c
           iter = iter + 1
           if (real(z)*real(z)+aimag(z)*aimag(z) > 4) then
              numoutside = numoutside + 1
!              print *, omp_get_thread_num(), omp_get_ancestor_thread_num(omp_get_level()-1)
              exit
           endif
        end do
     end do
  end do
  !$omp end parallel do

  finish = omp_get_wtime()

  area = 2.0*2.5*1.125 * real(npoints*npoints-numoutside)/real(npoints*npoints)
  error = area/real(npoints)
  print *, "Area of Mandelbrot set = ",area," +/- ",error
  print *, "Time  = ", finish-start, " seconds"
end
