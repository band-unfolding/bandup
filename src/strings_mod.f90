!    Copyright (C) 2013 Paulo V. C. Medeiros
!
!    This file is part of BandUP: Band Unfolding code for Plane-wave based calculations.
!
!    BandUP is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    BandUP is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with BandUP.  If not, see <http://www.gnu.org/licenses/>.
module strings
implicit none
PRIVATE
PUBLIC :: lower_case, upper_case, split_string 
contains
!*******************************************************************************************************
function lower_case(input_string) result(rtn)
! Based on the extended ASCII table (see http://www.ascii-code.com)
implicit none
character(len=*) :: input_string
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
!*******************************************************************************************************
!*******************************************************************************************************
function upper_case(input_string) result(rtn)
! Based on the extended ASCII table (see http://www.ascii-code.com)
implicit none
character(len=*) :: input_string
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
!*******************************************************************************************************
!*******************************************************************************************************
subroutine compact(io_string)
implicit none
character(len=*), intent(inout) :: io_string
integer i,num_spaces_excluded,end_of_io_string,selected_position,new_end_of_io_string


io_string=adjustl(io_string)
num_spaces_excluded=0
end_of_io_string=len_trim(io_string)
do i=1,end_of_io_string-1
  selected_position=i-num_spaces_excluded
  new_end_of_io_string=end_of_io_string-num_spaces_excluded
  if((io_string(selected_position:selected_position)==' ').and.(io_string(selected_position+1:selected_position+1)==' '))then
    io_string(selected_position:new_end_of_io_string-1)=io_string(selected_position+1:new_end_of_io_string)
    num_spaces_excluded=num_spaces_excluded+1
  endif
enddo
io_string(new_end_of_io_string+1:end_of_io_string)=' '

end subroutine compact
!*******************************************************************************************************
!*******************************************************************************************************
subroutine split_string(working_string,splitted_string)
implicit none
character(len=*) :: working_string
character(len=*), dimension(:), intent(out) :: splitted_string
integer :: end_of_working_string,position_on_splitted_string,position_on_working_string


call compact(working_string)
end_of_working_string=len_trim(working_string)
splitted_string=''
position_on_splitted_string=1
do position_on_working_string=1,end_of_working_string
  if(working_string(position_on_working_string:position_on_working_string)/=' ')then
    splitted_string(position_on_splitted_string)=trim(splitted_string(position_on_splitted_string))//working_string(position_on_working_string:position_on_working_string)
  else
    position_on_splitted_string=position_on_splitted_string+1
  endif
enddo


end subroutine split_string
!*******************************************************************************************************

end module strings
