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

module lists_and_seqs
use constants_and_types
implicit none
PRIVATE
PUBLIC :: real_seq, integer_seq, append, list_index, kpts_line

interface append
  module procedure append_integer_list, append_character_list, &
                   append_complex_list, append_matrix_index_list
end interface append

interface list_index
    module procedure list_index_integer, list_index_matrix_indices
end interface list_index

CONTAINS

subroutine append_integer_list(item, list)
!! Copyright (C) 2014 Paulo V. C. Medeiros
!! Not fully tested, but works for BandUP
implicit none
integer, intent(in) :: item
integer, dimension(:), allocatable, intent(inout) :: list
integer, dimension(:), allocatable :: new_list

if(.not.allocated(list))then
    allocate(list(1:1))
    list(1) = item
else
    allocate(new_list(1:size(list)+1))
    new_list(1:size(list)) = list
    new_list(size(list)+1) = item
    deallocate(list)
    allocate(list(1:size(new_list)))
    list = new_list
    deallocate(new_list)
endif

end subroutine append_integer_list

subroutine append_matrix_index_list(item, list)
!! Copyright (C) 2017 Paulo V. C. Medeiros
!! Not fully tested, but works for BandUP
implicit none
integer, dimension(1:2), intent(in) :: item
type(MatrixIndices), dimension(:), allocatable, intent(inout) :: list
type(MatrixIndices), dimension(:), allocatable :: new_list

if(.not.allocated(list))then
    allocate(list(1:1))
    list(1)%indices = item
else
    allocate(new_list(1:size(list)+1))
    new_list(1:size(list)) = list
    new_list(size(list)+1)%indices = item
    deallocate(list)
    allocate(list(1:size(new_list)))
    list = new_list
    deallocate(new_list)
endif

end subroutine append_matrix_index_list

subroutine append_complex_list(item, list)
!! Copyright (C) 2017 Paulo V. C. Medeiros
!! Not fully tested, but works for BandUP
implicit none
complex(kind=kind_cplx_coeffs), intent(in) :: item
complex(kind=kind_cplx_coeffs), dimension(:), allocatable, intent(inout):: list
complex(kind=kind_cplx_coeffs), dimension(:), allocatable :: new_list

if(.not.allocated(list))then
    allocate(list(1:1))
    list(1) = item
else
    allocate(new_list(1:size(list)+1))
    new_list(1:size(list)) = list
    new_list(size(list)+1) = item
    deallocate(list)
    allocate(list(1:size(new_list)))
    list = new_list
    deallocate(new_list)
endif

end subroutine append_complex_list

subroutine append_character_list(item, list)
!! Copyright (C) 2014 Paulo V. C. Medeiros
!! Not fully tested, but works for BandUP
implicit none
character(len=*), intent(in) :: item
character(len=*), dimension(:), allocatable, intent(inout) :: list
character(len=len(list)), dimension(:), allocatable :: new_list
integer :: size_new_list

if(.not.allocated(list))then
    allocate(list(1:1))
    list(1) = item
else
    allocate(new_list(1:size(list)+1))
    new_list(1:size(list)) = list
    new_list(size(list)+1) = item
    deallocate(list)
    size_new_list = size(new_list)
    allocate(list(1:size_new_list))
    list = new_list
    deallocate(new_list)
endif

end subroutine append_character_list


function list_index_integer(item, list) result(rtn) 
! Index of first occurence of item in list
implicit none
integer :: rtn
integer, intent(in) :: item
integer, dimension(:), intent(in) :: list
integer :: i
logical :: found
    
    rtn = 0
    if(lbound(list,1)>0 .and. ubound(list,1)>=lbound(list,1))then
        found = .FALSE.
        do i=1,size(list)
            if(list(i) == item)then
                found = .TRUE.
                exit
            endif
        enddo
        if(found) rtn = i
    endif

end function list_index_integer

function list_index_matrix_indices(item, list) result(rtn) 
! Index of first occurence of item in list
implicit none
integer :: rtn
integer, dimension(1:2), intent(in) :: item
type(MatrixIndices), dimension(:), intent(in) :: list
integer :: i
logical :: found
    
    rtn = 0
    if(lbound(list,1)>0 .and. ubound(list,1)>=lbound(list,1))then
        found = .FALSE.
        do i=1,size(list)
            if(all(list(i)%indices==item))then
                found = .TRUE.
                exit
            endif
        enddo
        if(found) rtn = i
    endif

end function list_index_matrix_indices


subroutine real_seq(first_term,last_term,increment,return_list)
implicit none
real(kind=dp), dimension(:), allocatable, intent(out) :: return_list
real(kind=dp), intent(in) :: first_term,last_term,increment
integer :: i, n_terms

    n_terms = nint((last_term - first_term)/increment + 1.0_dp)
    allocate(return_list(1:n_terms))
    do i = 1, n_terms
         return_list(i) = first_term + real(i - 1, kind=dp)*increment
    enddo

end subroutine real_seq

subroutine integer_seq(first_term, last_term, increment, return_list)
implicit none
integer, dimension(:), allocatable, intent(out) :: return_list
integer, intent(in) :: first_term, last_term
integer, intent(in), optional :: increment
integer :: inc, i, n_terms

    inc = 1
    if(present(increment)) inc = increment
    n_terms = nint((last_term - first_term)/inc + 1.0_dp)
    allocate(return_list(1:n_terms))
    do i = 1, n_terms
         return_list(i) = first_term + (i - 1)*inc
    enddo

end subroutine integer_seq

function kpts_line(kstart,kend,nk) result(line)
!! Copyright (C) 2014 Paulo V. C. Medeiros
!! Not fully tested, but works for BandUP
implicit none
type(vec3d), dimension(:), allocatable :: line
real(kind=dp), dimension(1:3), intent(in) :: kstart, kend
integer, intent(in) :: nk
integer :: alloc_stat,ikpt
real(kind=dp), dimension(1:3) :: line_vec

    deallocate(line, stat=alloc_stat)
    allocate(line(1:nk))
    if(nk>1)then
        line_vec = (kend - kstart)/real(nk - 1, kind=dp)
    else
        line_vec = 0.0_dp
    endif
    do ikpt=1,nk
        line(ikpt)%coord(:) = kstart + real(ikpt - 1,kind=dp)*line_vec(:)
    enddo
    return

end function kpts_line

end module lists_and_seqs
