program test
use WorldBuilder

IMPLICIT NONE

  ! Declare the types which will be needed.
  REAL*8 :: temperature,x=120e3,y=500e3,z=0,depth=0
  INTEGER :: composition_number = 3
  INTEGER*8 :: random_number_seed = 1.0 !! use a random number seed larger than zero
  REAL*8 :: composition
  character(len=256) :: file_name = "../../tests/data/continental_plate.wb"//C_NULL_CHAR
  logical(1) :: has_output_dir = .false.
  character(len=256) :: output_dir = "../../doc/"//C_NULL_CHAR

  ! Show how to call the functions.
  CALL create_world(cworld, file_name, has_output_dir, output_dir, random_number_seed)

  write(*, *) '2d temperature:'
  CALL temperature_2d(cworld,x,z,depth,temperature)
  write(*, *) 'temperature in fortran = ', temperature

  write(*, *) '3d temperature:'
  CALL temperature_3d(cworld,x,y,z,depth,temperature)
  write(*, *) 'temperature in fortran = ', temperature

    write(*, *) '2d composition:'
  CALL composition_2d(cworld,x,z,depth,composition_number,composition)
  write(*, *) 'composition in fortran = ', composition

  write(*, *) '3d composition:'
  CALL composition_3d(cworld,x,y,z,depth,composition_number,composition)
  write(*, *) 'composition in fortran = ', composition

  CALL release_world(cworld)
END program
