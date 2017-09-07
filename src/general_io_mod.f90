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
! MODULE: general_io
!
!> @author
!> Paulo V C Medeiros, Linköping University
!
! DESCRIPTION: 
!> Provides routines and constants to perform general I/O.
!==============================================================================

module general_io
use constants_and_types
implicit none
SAVE
PRIVATE
PUBLIC :: available_io_unit, file_extension, filename_without_extension, &
          get_file_size_in_bytes, str_len, file_header_BandUP, &
          file_for_pc_reduced_to_prim_cell, file_for_SC_reduced_to_prim_cell, &
          compiler_version, compilation_time, get_git_info_compiled_files, &
          timestamp

character(len=str_len), parameter :: &
    file_for_pc_reduced_to_prim_cell=&
        "BandUP_suggestion_of_pc_for_your_reference_unit_cell.POSCAR", &
    file_for_SC_reduced_to_prim_cell=&
        "BandUP_suggestion_of_smaller_SC_based_on_your_input_SC.POSCAR"

!! Functions and subroutines
CONTAINS 

function file_header_BandUP() result(header)
implicit none
character(len=:), allocatable :: header

    header = "# File created by BandUP (" // &
             trim(adjustl(package_version)) // ') at ' // &
             timestamp()

end function file_header_BandUP


function compiler_version() result(compiler_info_string)
implicit none
character(len=127) :: compiler_info_string
integer :: compiler_major, compiler_minor, compiler_patch, &
           compiler_release, compiler_build

#ifdef __INTEL_COMPILER
    write(compiler_info_string, '(A,I0)') 'Intel Fortran (ifort) V', &
                                          __INTEL_COMPILER
#elif defined __GFORTRAN__
    compiler_major = __GNUC__
    compiler_minor = __GNUC_MINOR__
    compiler_patch = __GNUC_PATCHLEVEL__
    write(compiler_info_string, '(A,2(I0,A),I0)') 'GNU Fortran (gfortran) V', &
        compiler_major,'.',compiler_minor,'.',compiler_patch
#elif defined NAGFOR
    compiler_release = 0
    compiler_build = 0
#   if defined __NAG_COMPILER_RELEASE
        compiler_release = __NAG_COMPILER_RELEASE
#   endif
#   if defined __NAG_COMPILER_BUILD
        compiler_build = __NAG_COMPILER_BUILD
#   endif
    write(compiler_info_string, '(2(A,I0))') 'NAG Fortran Release ', &
                                             compiler_release, &
                                             " Build ", &
                                             compiler_build
#else
    compiler_info_string='Unknown Compiler'
#endif

end function compiler_version

function compilation_time() result(timestamp)
implicit none
character(len=30) :: timestamp

#if defined (__INTEL_COMPILER) && defined (__TIMESTAMP__)
    write(timestamp, '(A)') __TIMESTAMP__
#elif defined (__GFORTRAN__) && defined (__DATE__) && defined (__TIME__)
    write(timestamp, '(A,X,A)') __DATE__, __TIME__
#elif defined __BUILD_START_TIME__
    write(timestamp, '(A)') __BUILD_START_TIME__
#else
    timestamp = 'unknown date and time'
#endif

end function compilation_time

function get_git_info_compiled_files(requested_info) result(req_info_val)
implicit none
character(len=*), intent(in) :: requested_info
character(len=127) :: req_info_val

    select case (trim(adjustl(requested_info)))
        case ('hash_latest_commit')
#           if defined (__COMMIT_HASH__)
                write(req_info_val, '(A)') __COMMIT_HASH__
#           else
                req_info_val = '(unknown hash from latest git commit)'
#           endif
        case ('branch_name')
#           if defined (__USED_BRANCH__)
                write(req_info_val, '(A)') __USED_BRANCH__
#           else
                req_info_val = 'unknown branch'
#           endif
        case default
            req_info_val = ''
    end select

end function get_git_info_compiled_files


function str_month(imonth) result(rtn)
implicit none
character(len=:), allocatable :: rtn
integer, intent(in) :: imonth
character(len=3), dimension(12) :: month_int2str

    month_int2str = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', &
                     'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    select case (imonth)
        case(1:12)
            rtn = month_int2str(imonth)
        case default
            rtn = 'ERR'
    end select

end function str_month

function timestamp() result(str_time)
implicit none
character(len=:), allocatable :: str_time
character(len=8) :: date
character(len=10) :: time
character(len=5) :: zone
integer, dimension(1:8) :: values
integer :: year, day, hour, minute
character(len=127) :: fmt_str, temp_str

    call date_and_time(date, time, zone, values)
    year = values(1)
    day = values(3)
    hour = values(5)
    minute = values(6)
    fmt_str = "(I0,':',I0.2,X,'UTC',A,X,'on',X,A,X,I0.2,',',X,I4)"
    write(temp_str, trim(fmt_str)) &
        hour, minute, zone, str_month(values(2)), day, year
    str_time = trim(adjustl(temp_str))

end function timestamp


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
integer(kind=long_int_kind), intent(out) :: file_size
character(len=*), intent(in) :: file
logical :: file_exists

    file_size = 0
    inquire(file=file, exist=file_exists)
    if(file_exists) inquire(file=file, size=file_size)
    return

end subroutine get_file_size_in_bytes


end module general_io
