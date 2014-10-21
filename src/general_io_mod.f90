!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
!!
!! This file is part of BandUP: Band Unfolding code for Plane-wave based calculations.
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

!===============================================================================
! MODULE: general_io
!
!> @author
!> Paulo V C Medeiros, Linköping University
!
! DESCRIPTION: 
!> Provides routines and constants to perform general I/O.
!===============================================================================

module general_io
use constants_and_types
use math
implicit none
SAVE
PRIVATE
PUBLIC :: available_io_unit, file_extension, filename_without_extension, get_file_size_in_bytes, & 
          str_len, package_version, file_header_BandUP, file_header_BandUP_short, &
          file_for_pc_reduced_to_prim_cell, file_for_SC_reduced_to_prim_cell

integer, parameter :: str_len=256
character(len=30), parameter :: package_version="2.5.0 (BETA), 2014-10-22"
character(len=str_len), parameter :: file_header_BandUP="# File created by BandUP - Band Unfolding code for Plane-wave based calculations, &
                                                        V"//trim(adjustl(package_version)), &
                                     file_header_BandUP_short="# File created by BandUP, V"//trim(adjustl(package_version)), &
                                     file_for_pc_reduced_to_prim_cell="BandUP_suggestion_of_pc_for_your_reference_unit_cell.POSCAR", &
                                     file_for_SC_reduced_to_prim_cell="BandUP_suggestion_of_smaller_SC_based_on_your_input_SC.POSCAR"

!! Functions and subroutines
CONTAINS 


function available_io_unit(min_unit,max_unit) result(unit_num)
! Returns a number unit_num which can be safely used in statements like
! open(file=unit_num)
implicit none
integer, intent(in),optional :: min_unit,max_unit
logical :: unit_is_opened
integer ::unit_num,int_aux,min_unit_num,max_unit_num

  if(present(min_unit).and.present(max_unit))then
    min_unit_num=min_unit
    max_unit_num=max_unit
    if(min_unit_num>max_unit_num)then
      int_aux=max_unit_num
      max_unit_num=min_unit_num
      min_unit_num=int_aux
    else
      if(min_unit_num==max_unit_num)then
        min_unit_num=10
        max_unit_num=90
      endif
    endif
  else
    min_unit_num=10
    max_unit_num=90
  endif

  unit_num=min_unit_num
  inquire(unit=unit_num,opened=unit_is_opened)
  do while((unit_is_opened).and.(unit_num<max_unit_num))
    unit_num=unit_num+1
    inquire(unit=unit_num,opened=unit_is_opened)
  enddo

end function available_io_unit

function file_extension(filename, extention_length) result(ext)
implicit none
character(len=:), allocatable :: ext
character(len=*), intent(in) :: filename
integer, intent(in), optional :: extention_length
integer :: ext_length, dot_position

    ext_length = 3
    if(present(extention_length))then
        ext_length = extention_length
    endif
    allocate(character(ext_length) :: ext)

    ext = ''
    if(len(filename) >= ext_length + 1)then
        dot_position = scan(filename, '.', back=.TRUE.)
        if(dot_position > 0)then
            ext = trim(adjustl(filename(dot_position + 1:len(filename))))
        endif
    endif

end function file_extension


function filename_without_extension(filename) result(rtn)
implicit none
character(len=*), intent(in) :: filename
character(len=len(filename)) :: rtn
integer :: dot_position

    if(trim(adjustl(filename))=='')then
        rtn = ''
        return
    endif

    dot_position = scan(filename, '.', back=.TRUE.)
    if(dot_position > 0)then
       rtn = filename(1:dot_position-1)
    else
        rtn = trim(adjustl(filename))
    endif

end function filename_without_extension


subroutine get_file_size_in_bytes(file_size, file)
implicit none
integer(8), intent(out) :: file_size
character(len=*), intent(in) :: file
logical :: file_exists

    file_size = 0
    inquire(file=file, exist=file_exists)
    if(file_exists) inquire(file=file, size=file_size)
    return

end subroutine get_file_size_in_bytes


end module general_io
