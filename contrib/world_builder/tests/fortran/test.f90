program test
use WorldBuilder

IMPLICIT NONE

  print *, "hello world!"
  ! Declare the types which will be needed.
  
  !REAL*8 :: temperature,x=120e3,y=500e3,z=0,depth=0,gravity = 10
  !INTEGER :: composition_number = 3
  !REAL*8 :: composition
  !character(len=256) :: path
  !INTEGER*4 :: k = 1
  !character(len=256) :: file_name != MY_FLAG//"/data/continental_plate.wb"//C_NULL_CHAR
  !logical(1) :: has_output_dir = .false.
  !character(len=256) :: output_dir = "../../../doc/"//C_NULL_CHAR

  !call getarg( k, file_name )
!  file_name = trim(file_name//C_NULL_CHAR
  ! Show how to call the functions.
!  CALL create_world(cworld, trim(file_name)//C_NULL_CHAR, has_output_dir, output_dir)

 ! write(*, *) '2d temperature:'
  !CALL temperature_2d(cworld,x,z,depth,gravity,temperature)
 ! write(*, *) 'temperature in fortran = ', temperature

  !write(*, *) '3d temperature:'
  !CALL temperature_3d(cworld,x,y,z,depth,gravity,temperature)
  !write(*, *) 'temperature in fortran = ', temperature

  !  write(*, *) '2d composition:'
  !CALL composition_2d(cworld,x,z,depth,composition_number,composition)
  !write(*, *) 'composition in fortran = ', composition

  !write(*, *) '3d composition:'
  !CALL composition_3d(cworld,x,y,z,depth,composition_number,composition)
  !write(*, *) 'composition in fortran = ', composition

  !CALL release_world(cworld)
END program
