program hello_world
  implicit none

  !$omp parallel
  print *, "Hello World"
  !$omp end parallel
end
