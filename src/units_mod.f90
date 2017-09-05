!! Copyright (C) 2013-2017 Paulo V. C. Medeiros
!!
!! This file is part of BandUP: 
!! Band Unfolding code for Plane-wave based calculations.
!!
!! BandUP is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! BandUP is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with BandUP.  If not, see <http://www.gnu.org/licenses/>.

!==============================================================================
! MODULE: math 
!
!> @author
!> Paulo V C Medeiros, Linköping University
!
! DESCRIPTION: 
!> Defines some mathematical functions and subroutines used by BandUP.
!==============================================================================

module units
use constants_and_types
implicit none
PRIVATE
PUBLIC :: to_ev, to_angstrom

CONTAINS

function to_angstrom(value, units) result(rtn)
implicit none
real(kind=dp) :: rtn
real(kind=dp), intent(in) :: value
character(len=*), intent(in) :: units
real(kind=dp) :: conv_factor
character(len=1) :: first_letter_units

    first_letter_units = trim(adjustl(units))
    select case(first_letter_units)
        case default 
            conv_factor = 1.0_dp
        case('b', 'B')
            conv_factor = bohr
    end select
    rtn = conv_factor * value

end function to_angstrom


function to_ev(value, units) result(rtn)
implicit none
real(kind=dp) :: rtn
real(kind=dp), intent(in) :: value
character(len=*), intent(in) :: units
real(kind=dp) :: conv_factor
real(kind=dp), parameter :: Ry =  13.60569172, Hartree = 27.2113834
character(len=1) :: first_letter_units

    first_letter_units = trim(adjustl(units))
    select case(first_letter_units)
        case default 
            conv_factor = 1.0_dp
        case('r', 'R')
            conv_factor = Ry
        case('h', 'H')
            conv_factor = Hartree
    end select
    rtn = conv_factor * value

end function to_ev

end module units
