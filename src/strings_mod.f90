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
! MODULE: strings
!
!> @author
!> Paulo V C Medeiros, Linköping University
!
! DESCRIPTION: 
!> Provides basic methods to operate with strings.
!==============================================================================

module strings
use constants_and_types
implicit none
PRIVATE
PUBLIC :: lower_case, upper_case, split, compact

!==============================================================================
!> Splits the input string and puts the result in the output array.
!> Assumes spaces as separators. 
!> If the output array is of type *type*, then the interface will assume that 
!! all elements in the input string can be converted to the type *type*
!!
!! @param[in] string Input string
!! @param[out] splitted_string Output array
!==============================================================================
interface split
  module procedure split_string_into_strings, split_string_into_integers, &
                   split_string_into_floats
end interface split

contains

!==============================================================================
!> Returns a lower-case version of the input string.
!! Based on the extended ASCII table (see http://www.ascii-code.com)
!!
!! @param[in] input_string String of arbitrary length.
!==============================================================================
function lower_case(input_string) result(rtn)
implicit none
character(len=*), intent(in) :: input_string
character(len=len(input_string)) :: rtn
integer, parameter :: upper_to_lower=ichar('a')-ichar('A')
integer :: i

   do i=1, len(input_string)
      select case (ichar(input_string(i:i)))
         case(65:90, 192:214, 216:222)
           rtn(i:i)=achar(ichar(input_string(i:i))+upper_to_lower)
         case default
           rtn(i:i)=input_string(i:i)
      end select
   enddo

end function lower_case

!==============================================================================
!> Returns an upper-case version of the input string.
!! Based on the extended ASCII table (see http://www.ascii-code.com)
!!
!! @param[in] input_string String of arbitrary length.
!==============================================================================
function upper_case(input_string) result(rtn)
implicit none
character(len=*), intent(in) :: input_string
character(len=len(input_string)) :: rtn
integer, parameter :: lower_to_upper=ichar('A')-ichar('a')
integer :: i

   do i=1, len(input_string)
      select case (ichar(input_string(i:i)))
         case(97:122, 224:246, 248:254)
           rtn(i:i)=achar(ichar(input_string(i:i))+lower_to_upper)
         case default
           rtn(i:i)=input_string(i:i)
      end select
   enddo

end function upper_case

!==============================================================================
!> Returns the number of occurancies of a given substring in a given string.
!!
!! @param[in] string Input string
!! @param[in] substring Substring to be searched
!==============================================================================
recursive function count_occurencies_substring(string, substring) &
result(occ_count)
!! Copyright (C) 2014 Paulo V. C. Medeiros
implicit none
character(len=*), intent(in) :: substring, string
integer :: occ_count, new_start, pos_substr_in_str

    occ_count = 0
    pos_substr_in_str = index(string, substring)
    if(pos_substr_in_str > 0)then
        occ_count = occ_count + 1
        new_start = pos_substr_in_str + len(substring)
        if(new_start <= len(string))then
            occ_count = occ_count + &
                        count_occurencies_substring(&
                            string(new_start:), substring &
                        )
        endif
    endif 

end function count_occurencies_substring

!==============================================================================
!> Removes leading spaces, trailing blanks and double spaces from the input
!! string.
!!
!! @param[in, out] io_string String of arbitrary length.
!==============================================================================
subroutine compact(io_string)
implicit none
character(len=*), intent(inout) :: io_string
integer :: i, num_spaces_excluded, end_of_io_string, selected_position, &
           new_end_of_io_string

    io_string = trim(adjustl(io_string))
    if(len_trim(io_string) < 2)then
        return
    endif
    io_string = trim(adjustl(io_string))
    num_spaces_excluded = 0
    end_of_io_string = len_trim(io_string)
    do i=1,end_of_io_string-1
        selected_position = i - num_spaces_excluded
        new_end_of_io_string = end_of_io_string-num_spaces_excluded
        if((io_string(selected_position:selected_position)==' ') .and. & 
           (io_string(selected_position + 1:selected_position + 1)==' '))then
            io_string(selected_position:new_end_of_io_string - 1) =  &
                io_string(selected_position + 1:new_end_of_io_string)
            num_spaces_excluded = num_spaces_excluded + 1
        endif
    enddo
    io_string(new_end_of_io_string + 1:end_of_io_string) = ''

end subroutine compact

!==============================================================================
!> Splits the input string and puts the result in an array of strings.
!> Assumes spaces as separators.
!!
!! @param[in] string Input string
!! @param[out] splitted_string Output array of strings
!==============================================================================
subroutine split_string_into_strings(string, splitted_string)
!! Copyright (C) 2014 Paulo V. C. Medeiros
implicit none
character(len=*), intent(in) :: string
character(len=*), dimension(:), allocatable, intent(out) :: splitted_string
character(len=len(string)) :: aux_string
integer :: n_components, alloc_stat, icomp

    aux_string = string
    call compact(aux_string)
    deallocate(splitted_string,stat=alloc_stat)
    if(trim(adjustl(aux_string)) == '')then
        allocate(splitted_string(1:1))
        splitted_string = ''   
    else
        n_components = 1 + &
                       count_occurencies_substring(&
                           string=trim(adjustl(aux_string)), substring=' ' &
                       )
        allocate(splitted_string(1:n_components))
        read(aux_string,*)(splitted_string(icomp), icomp=1,n_components)
    endif

end subroutine split_string_into_strings

!==============================================================================
!> Splits the input string and puts the result in an array of integers.
!> Assumes spaces as separators. Assumes that the input string contains only
!! integer numbers. 
!!
!! @param[in] string Input string
!! @param[out] splitted_string Output array of integers
!==============================================================================
subroutine split_string_into_integers(string, splitted_string)
!! Copyright (C) 2014 Paulo V. C. Medeiros
implicit none
character(len=*), intent(in) :: string
integer, dimension(:), allocatable, intent(out) :: splitted_string
character(len=len(string)) :: aux_string
integer :: n_components, alloc_stat, icomp

    aux_string = string
    call compact(aux_string)
    deallocate(splitted_string,stat=alloc_stat)
    if(trim(adjustl(aux_string)) == '')then
        write(*,*)'WARNING (split_string_into_integers): Empty input string. &
                   Splitted string not allocated.'
        return
    else
        n_components = 1 + &
                       count_occurencies_substring(&
                           string=trim(adjustl(aux_string)), substring=' ' &
                       )
        allocate(splitted_string(1:n_components))
        read(aux_string,*)(splitted_string(icomp), icomp=1,n_components)
    endif

end subroutine split_string_into_integers

!==============================================================================
!> Splits the input string and puts the result in an array of floats.
!> Assumes spaces as separators. Assumes that the input string contains only
!! real numbers. 
!!
!! @param[in] string Input string
!! @param[out] splitted_string Output array of floats
!==============================================================================
subroutine split_string_into_floats(string, splitted_string)
!! Copyright (C) 2014 Paulo V. C. Medeiros
implicit none
character(len=*), intent(in) :: string
real(kind=dp), dimension(:), allocatable, intent(out) :: splitted_string
character(len=len(string)) :: aux_string
integer :: n_components, alloc_stat, icomp

    aux_string = string
    call compact(aux_string)
    deallocate(splitted_string,stat=alloc_stat)
    if(trim(adjustl(aux_string)) == '')then
        write(*,*)'WARNING (split_string_into_floats): Empty input string. &
                   Splitted string not allocated.'
        return
    else
        n_components = 1 + &
                       count_occurencies_substring(&
                           string=trim(adjustl(aux_string)), substring=' ' &
                       )
        allocate(splitted_string(1:n_components))
        read(aux_string,*)(splitted_string(icomp), icomp=1,n_components)
    endif

end subroutine split_string_into_floats


end module strings
