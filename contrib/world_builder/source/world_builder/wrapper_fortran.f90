!!
!!  Copyright (C) 2018-2024 by the authors of the World Builder code.
!!
!!  This file is part of the World Builder.
!!
!!  This program is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU Lesser General Public License as published
!!  by the Free Software Foundation, either version 2 of the License, or
!!  (at your option) any later version.
!!
!!  This program is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU Lesser General Public License for more details.
!!
!!  You should have received a copy of the GNU Lesser General Public License
!!  along with this program.  If not, see <https://www.gnu.org/licenses/>.
!!

MODULE WorldBuilder
USE, INTRINSIC :: ISO_C_BINDING!, ONLY: C_PTR
  IMPLICIT NONE

  !> This contains the interface to the world builder for fortran.
  INTERFACE
  !> Create an interface with the create world C function.
  !! This function creates an object of the world builder and returns a pointer
  !! to it. This pointer can then be used to call the temperature and composition
  !!functions. When done call the release world function to destroy the object.
  !! Please not that the basic type for the random number seed is a unsigned long,
  !! so it is advisable to not use negative numbers if you want to reproduce the 
  !! same result as other codes.
   SUBROUTINE create_world(cworld, file_name, has_output_dir, output_dir, random_number_seed) BIND(C, NAME='create_world') 
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_BOOL, C_CHAR, C_LONG
      IMPLICIT NONE
      ! This argument is a pointer passed by reference.
      TYPE(C_PTR), INTENT(OUT) :: cworld
      character(KIND=C_CHAR,len=1),  intent(in)  :: file_name
      logical(KIND=C_BOOL), intent(in) :: has_output_dir
      character(KIND=C_CHAR,len=1),  intent(in)  :: output_dir
      INTEGER(C_LONG), intent(in), value ::random_number_seed
    END SUBROUTINE create_world

    !> Create an interface with the 2d temperature C function of the World builder.
    !! This function return the temperature at a specific location given x, z and depth.
    SUBROUTINE temperature_2d(cworld, x, z, depth, temperature) BIND(C, NAME='temperature_2d')
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      ! This argument is a pointer passed by value.
      TYPE(C_PTR), INTENT(IN), VALUE :: cworld
      REAL(C_DOUBLE), intent(in), value :: x
      REAL(C_DOUBLE), intent(in), value :: z
      REAL(C_DOUBLE), intent(in), value :: depth
      REAL(C_DOUBLE), intent(out) :: temperature
    END SUBROUTINE temperature_2d

    !> Create an interface with the 3d temperature function of the World builder.
    !! This function return the temperature at a specific location given x, y, z and depth.
        SUBROUTINE temperature_3d(cworld, x, y, z, depth, temperature) BIND(C, NAME='temperature_3d')
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      ! This argument is a pointer passed by value.
      TYPE(C_PTR), INTENT(IN), VALUE :: cworld
      REAL(C_DOUBLE), intent(in), value :: x
      REAL(C_DOUBLE), intent(in), value :: y
      REAL(C_DOUBLE), intent(in), value :: z
      REAL(C_DOUBLE), intent(in), value :: depth
      REAL(C_DOUBLE), intent(out) :: temperature
    END SUBROUTINE temperature_3d

    !> Create an interface with the 2d composition function of the World builder.
    !! This function return the composition at a specific location given x, z, depth and
    !! composition number.
        SUBROUTINE composition_2d(cworld, x, z, depth, composition_number, composition) BIND(C, NAME='composition_2d')
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      ! This argument is a pointer passed by value.
      TYPE(C_PTR), INTENT(IN), VALUE :: cworld
      REAL(C_DOUBLE), intent(in), value :: x
      REAL(C_DOUBLE), intent(in), value :: z
      REAL(C_DOUBLE), intent(in), value :: depth
      INTEGER(C_INT), intent(in), value :: composition_number
      REAL(C_DOUBLE), intent(out) :: composition
    END SUBROUTINE composition_2d

    !> Create an interface with the 3d composition function of the World builder.
    !! This function return the composition at a specific location given x, y, z, depth and
    !! composition number.
      SUBROUTINE composition_3d(cworld, x, y, z, depth, composition_number, composition) BIND(C, NAME='composition_3d')
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      ! This argument is a pointer passed by value.
      TYPE(C_PTR), INTENT(IN), VALUE :: cworld
      REAL(C_DOUBLE), intent(in), value :: x
      REAL(C_DOUBLE), intent(in), value :: y
      REAL(C_DOUBLE), intent(in), value :: z
      REAL(C_DOUBLE), intent(in), value :: depth
      INTEGER(C_INT), intent(in), value :: composition_number
      REAL(C_DOUBLE), intent(out) :: composition
    END SUBROUTINE composition_3d

    !> Create an interface with the release world function.
    !! This is the destructor for the world builder class. Call this function when done
    !! with the world builder.
    SUBROUTINE release_world(cworld) BIND(C, NAME='release_world')
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR
      IMPLICIT NONE
      ! This argument is a pointer passed by value.
      TYPE(C_PTR), INTENT(IN), VALUE :: cworld
    END SUBROUTINE release_world
  END INTERFACE

  !> The C pointer to the World Builder world. It is generated by the create_world function.
  TYPE(C_PTR) :: cworld
  END MODULE WorldBuilder
