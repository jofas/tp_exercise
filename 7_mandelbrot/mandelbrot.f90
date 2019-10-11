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

  !$omp parallel default(shared)

  !!$omp do private(i, j, c, z, iter) reduction(+:numoutside)

  !$omp master
  do i = 0, npoints-1
     do j= 0,npoints-1
        !$omp task default(none) private(c, z, iter) firstprivate(i, j) shared(numoutside)
        c = cmplx(-2.0+(2.5*i)/npoints + 1.0d-07,(1.125*j)/npoints + 1.0d-07)
        z = c
        iter = 0
        do while (iter < maxiter)
           z = z*z + c
           iter = iter + 1
           if (real(z)*real(z)+aimag(z)*aimag(z) > 4) then
             !$omp atomic
             numoutside = numoutside + 1
              exit
           endif
        end do
        !$omp end task
     end do
  end do
  !$omp end master

  !!$omp end do

  !$omp end parallel

  finish = omp_get_wtime()

  area = 2.0*2.5*1.125 * real(npoints*npoints-numoutside)/real(npoints*npoints)
  error = area/real(npoints)
  print *, "Area of Mandelbrot set = ",area," +/- ",error
  print *, "Time  = ", finish-start, " seconds"
end
