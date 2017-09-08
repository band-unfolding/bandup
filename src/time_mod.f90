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

module time
use constants_and_types
!$ use omp_lib
implicit none
PRIVATE
PUBLIC :: time_now, initialize, formatted_time

CONTAINS

function time_now() result(rtn)
implicit none
real(kind=dp) :: rtn
integer :: count, count_rate

!$ rtn = omp_get_wtime()
!$ return
   CALL system_clock(count, count_rate)
   rtn = count/count_rate

end function time_now


subroutine initialize(times)
implicit none
type(timekeeping), intent(out) :: times

    times%start=time_now()
    times%read_wf=0.0_dp
    times%calc_spec_weights=0.0_dp
    times%calc_SF=0.0_dp
    times%calc_dN=0.0_dp
    times%calc_rho=0.0_dp
    times%calc_pauli_vec=0.0_dp
    times%calc_pauli_vec_projs=0.0_dp
    times%write_dN_files=0.0_dp
    times%write_unf_dens_op_files=0.0_dp
    times%end = time_now()

end subroutine initialize


function formatted_time(t_in_sec) result(rtn)
implicit none
character(len=127) :: rtn
real(kind=dp), intent(in) :: t_in_sec
real(kind=dp) :: unaccounted_t, secs
integer :: days, hours, minutes
character(len=127) :: aux_char

    unaccounted_t = t_in_sec
    days = int(unaccounted_t/86400.0_dp)

    unaccounted_t = unaccounted_t - real(days, kind=dp) * 86400.0_dp
    hours = int(unaccounted_t/3600.0_dp)

    unaccounted_t = unaccounted_t - real(hours, kind=dp) * 3600.0_dp

    minutes = int(unaccounted_t/60.0_dp)
    unaccounted_t = unaccounted_t - real(minutes, kind=dp) * 60.0_dp

    secs = unaccounted_t
   
    rtn = ''
    if(days > 0)then
        write(aux_char,'(I0,A1)') days, 'D '
        write(rtn,'(A)') trim(adjustl(aux_char))
    endif
    if(hours > 0)then
        write(aux_char,'(I0,A1)') hours, 'h '
        write(rtn,'(2A)')rtn(1:len_trim(rtn) + 1), trim(adjustl(aux_char))
    endif     
    if(minutes > 0)then
        write(aux_char,'(I0,A1)') minutes, 'm'
        write(rtn,'(2A)')rtn(1:len_trim(rtn) + 1), trim(adjustl(aux_char))
    endif     
    if(secs > epsilon(1.0_dp))then
        write(aux_char,'(f0.2,A1)') secs, 's'
        write(rtn,'(2A)')rtn(1:len_trim(rtn) + 1), trim(adjustl(aux_char))
    endif
    
    rtn = trim(adjustl(rtn)) 

end function formatted_time

end module time
