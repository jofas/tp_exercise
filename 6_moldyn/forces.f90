subroutine forces(npart, x, f, vir, epot, side, rcoff)
  use mod_setprecision
  use omp_lib
  implicit none
  integer :: npart
  real (kind=dprec), intent(in)    :: x(npart, 3)
  real (kind=dprec), intent(inout) :: f(npart, 3)
  real (kind=dprec), intent(inout) :: vir, epot
  real (kind=dprec), intent(in) :: side, rcoff

  integer :: i, j
  real (kind=dprec) :: rcoffs, sideh, xi, yi, zi, fxi, fyi, fzi, xx, yy, zz
  real (kind=dprec):: rd, rrd, rrd2, rrd3, rrd4, rrd6, rrd7, r148
  real (kind=dprec):: forcex, forcey, forcez

  sideh = 0.5d0*side
  rcoffs = rcoff*rcoff

  !$omp do &
  !$omp private(xi, yi, zi, fxi, fyi, fzi, xx, yy, zz, rd, &
  !$omp         rrd, rrd2, rrd3, rrd4, rrd6, rrd7, r148, &
  !$omp         forcex, forcey, forcez) &
  !$omp reduction(+:epot, f, vir)
  do i = 1,npart
     xi = x(i,1)
     yi = x(i,2)
     zi = x(i,3)
     fxi = 0.0d0
     fyi = 0.0d0
     fzi = 0.0d0
     do j = i+1,npart
        xx = xi - x(j,1)
        yy = yi - x(j,2)
        zz = zi - x(j,3)
        if (xx < -sideh) xx = xx + side
        if (xx >  sideh) xx = xx - side
        if (yy < -sideh) yy = yy + side
        if (yy >  sideh) yy = yy - side
        if (zz < -sideh) zz = zz + side
        if (zz >  sideh) zz = zz - side
        rd = xx*xx + yy*yy + zz*zz
        if (rd <= rcoffs) then
           rrd = 1.0d0/rd
           rrd2 = rrd*rrd
           rrd3 = rrd2*rrd
           rrd4 = rrd2*rrd2
           rrd6 = rrd2*rrd4
           rrd7 = rrd6*rrd
           r148 = rrd7 - 0.5d0*rrd4

           epot = epot + (rrd6 - rrd3)
           vir = vir - rd*r148

           forcex = xx * r148
           forcey = yy * r148
           forcez = zz * r148
           fxi = fxi + forcex
           fyi = fyi + forcey
           fzi = fzi + forcez
           f(j,1) = f(j,1) - forcex
           f(j,2) = f(j,2) - forcey
           f(j,3) = f(j,3) - forcez
        endif
     end do
     f(i,1) = f(i,1) + fxi
     f(i,2) = f(i,2) + fyi
     f(i,3) = f(i,3) + fzi
  end do
  !$omp end do
end subroutine forces
