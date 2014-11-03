!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros 
!! =============================================================================================================
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

!! This module has been adapted from the Fortran 90 version of the WaveTrans code
!! =============================================================================================================
!! WaveTrans, Copyright (C) 2012 R. M. Feenstra and M. Widom <http://www.andrew.cmu.edu/user/feenstra/wavetrans>
!! =============================================================================================================
!! Permission to distribute this module as a part of BandUP and under the terms of the GNU GPL licence (or similar)
!! has been granted by WaveTrans' Copyright holders. Thanks!

!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module read_vasp_wavecar
! ===================
!  Contains a function "read_wavefunction", which returns information about any
!  selected k-point, skipping the others. This is very important for the
!  performance of the BandUP code (or any other code that would require similar functionality).
!  See the definition of the function for more information on the arguments.
!  I've only tested this module on linux and using the ifort compiler.
!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$************************* List of major changes from the WaveTrans code ****************************************
!!$*
!!$*       1: Turned the code into a module (read_wavecar), which can be used in the  program.
!!$*          1.a: The routine "read_wavecar" returns information about any selected k-point, skipping the others.
!!$*               In order to do that, I have changed a bit the way the positioning of the 
!!$*               register is handled (which is also needed for openmp, see item 2)
!!$*          1.b: The GCOEFF.txt file is no longer created
!!$*          1.c: Started using "implicit none"
!!$*               This made it possible, for instance, to spot a typo (npreci) typo of the variable 'nprec'
!!$*          1.d: Added the variable "iostat", which returns 0 if everything is OK and -1 if there is an error.
!!$*          1.e: Started using the constants "sp" and "dp" to controll the precision of the real variables 
!!$*               Replaced "real*8" by "real(kind=dp)" 
!!$*               Replaced "complex*8" by "complex(kind=sp)", ie., single precision real and imaginary parts 
!!$*               Replaced "complex*16" by "complex(kind=dp)", ie., double precision real and imaginary parts
!!$*          1.f: Replaced subroutine vcross(vtemp, v1, v2) by function cross(v1, v2)
!!$*          1.g: Added fuction norm(v) to compute the norm of a vector instead of doing it explicitly every time
!!$*          1.h: All the angles (phi12, ..., phi123) are now calculated using the function "angle"
!!$*          1.i: Added an option to renormalize the WF to unity w.r.t. the usual inner product
!!$*               This allows BandUP to get accurate spectral weights as the norm of the partial functions,
!!$*               and ensures that the calculated unfolding delta_Ns are always (non-negative) integers
!!$*               when using a perfect supercell, as it should be [see Phys. Rev. B 89, 041407(R) (2014)].
!!$*
!!$*       2: openmp is now used to parallelize (over bands) the reading of the WAVECAR file
!!$*          2.a: Introduced the function "available_io_unit"
!!$*               Each thread associates the WAVECAR file to a different IO unit in the "open" statements.
!!$*
!!$*       3: Added support for spinor wavefunctions
!!$*          3.a: Started using "access='stream'" instead of "access='direct'" to read the coefficients.     
!!$*               It made it easier reading the spinor components into separate sections of the coeff array 
!!$*          
!!$*       4: Some small rorganization of the code that are not so important, but are worth mentioning
!!$*          4.a: The size of the WAVECAR file is now passed to the user (optional)
!!$*          4.b: The constant "c" (vasp_2m_over_hbar_sqrd) is now defined by means of a "parameter" statement 
!!$*               instead of the original "data"
!!$*          4.c: Replaced uni=6 by unit=* in the write statements
!!$*          4.d: Pi is now declared as a parameter, and a parameter "twopi" has now been defined as well.
!!$*          4.e: The names of some variables have been changed (just because of my personal style of writing)
!!$*
!!$******* Below, comments in the header of the WaveTrans code (adapted) ****************************************       
!!$*       
!!$*   input the WAVECAR file in binary format from VASP 
!!$*   Compile with gfortran or ifort. Flag "-assume byterecl" is required
!!$*   for ifort.
!!$*
!!$*   the x,y,z components of each G value are given in terms of the
!!$*   ig1,ig2,ig3 values and the components of the recip. lattice vectors
!!$*   according to vector_G = sum(ign*vector_bn. n = 1, 2, 3)i, i.e,
!!$*   ig1*b1_x + ig2*b2_x + ig3*b3_x,
!!$*   ig1*b1_y + ig2*b2_y + ig3*b3_y, and
!!$*   ig1*b1_z + ig2*b2_z + ig3*b3_z, respectively
!!$*
!!$*   note that the energy eigenvalues are complex, as provided in the
!!$*   WAVECAR file, but the imaginary part is zero (at least for all cases
!!$*   investigated thus far)
!!$******* End of comments in the header of the WaveTrans code **************************************************       

module read_vasp_wavecar
use constants_and_types
use cla_wrappers
use math
use general_io, only: available_io_unit
!$ use omp_lib
implicit none
PRIVATE
PUBLIC :: read_wavecar, get_wavecar_record_lenght, vasp_2m_over_hbar_sqrd
!!$* Comments adapted from WaveTrans: The value of vasp_2m_over_hbar_sqrd has been adjusted in final decimal places to agree with VASP value
!!$*                                  The routine "read_wavecar" checks for discrepancy of any results between this and VASP values
real(kind=dp), parameter :: vasp_2m_over_hbar_sqrd = 0.262465831 ! c = 2m/hbar**2 in units of 1/eV Ang^2

CONTAINS


subroutine get_wavecar_record_lenght(nrecl, wf_file, iost)
implicit none
integer(kind=selected_int_kind(18)), intent(out) :: nrecl
character(len=*), intent(in) :: wf_file
integer, intent(out) :: iost
integer :: input_file_unit
real(kind=dp) :: xnrecl

    input_file_unit = available_io_unit()
    nrecl=24 ! Trial record length
    open(unit=input_file_unit, file=wf_file, access='direct', recl=nrecl, iostat=iost, status='old')
        read(unit=input_file_unit, rec=1, iostat=iost) xnrecl
        nrecl = nint(xnrecl)
    close(unit=input_file_unit)

end subroutine get_wavecar_record_lenght


subroutine get_total_n_SC_kpts_wavecar(nkpts, file)
implicit none
integer, intent(out) :: nkpts
character(len=*), intent(in) :: file
integer :: input_file_unit, iost
integer(kind=selected_int_kind(18)) :: nrecl
real(kind=dp) :: xnwk

    nkpts = 0
    call get_wavecar_record_lenght(nrecl, file, iost)
    if(iost==0)then
        input_file_unit = available_io_unit()
        open(unit=input_file_unit, file=file, access='direct', recl=nrecl, iostat=iost, status='old')
            read(unit=input_file_unit, rec=2, iostat=iost) xnwk
        close(input_file_unit)
        if(iost==0) nkpts = nint(xnwk)
    endif

end subroutine get_total_n_SC_kpts_wavecar


subroutine read_wavecar(wf, file, ikpt, read_coeffs, iostat)
implicit none
! Mandatory input and output variables
type(pw_wavefunction), intent(inout) :: wf
character(len=*), intent(in) :: file
integer, intent(in) :: ikpt
! Optional input and output variables
integer, optional, intent(out) :: iostat
logical, optional, intent(in) :: read_coeffs
! Parameters
integer, parameter :: size_of_a_coeff_in_bytes = sizeof(cmplx(1, kind=kind_cplx_coeffs))
real(kind=dp), parameter :: c = vasp_2m_over_hbar_sqrd ! 2m/hbar**2 in units of 1/eV Ang^2
! Local variables
real(kind=dp), dimension(1:3) :: frac_coords_selected_kpt 
complex(kind=dp), allocatable :: cener(:)  ! Re{cener(iband)} = E(K; iband)
real(kind=dp), dimension(:), allocatable :: occ
real(kind=dp), dimension(1:3) :: a1, a2, a3, b1, b2, b3, sumkg
real(kind=dp), dimension(1:3,1:3) :: b_matrix 
real(kind=dp) :: b1mag, b2mag, b3mag, ecut, gtot, etot, phi12, phi13, phi23, sinphi123, &
                 xnwk, xnband, xnrecl, xnspin, xnprec, xnplane
integer :: i_kpt, nplane, i, j, iost,  ig1, ig2, ig3, ig1p, ig2p, ig3p, iplane, ispin, &
           irec_start_first_spin, irec_start_chosen_spin, nwk, nprec, nband, &
           nb1maxA, nb2maxA, nb3maxA, nb1maxB, nb2maxB, nb3maxB, nb1maxC, nb2maxC, nb3maxC, &
           npmaxA, npmaxB, npmaxC, nb1max, nb2max, nb3max, npmax, iband, ncnt,  &
           input_file_unit, irec_before_loop,temp_unit, alloc_stat, i_wf_comp, &
           n_threads, ithread
integer(kind=selected_int_kind(18)) :: stream_pos, irec, nrecl
integer, dimension(:), allocatable :: io_unit_for_thread
complex(kind=kind_cplx_coeffs) :: inner_prod
logical :: read_coefficients, reset_spin

! Start
read_coefficients = .TRUE. ! Reading the coeffs by default
if(present(read_coeffs)) read_coefficients = read_coeffs
if(present(iostat)) iostat = -1

call get_wavecar_record_lenght(nrecl, file, iost)
if(iost /= 0)then
    write(*,'(A)')'ERROR (read_wavecar): Could not determine record legth for the file.'
    write(*,'(A)')'Stopping now.'     
    if(present(iostat)) iostat = iost
    stop
endif

input_file_unit = available_io_unit()
open(unit=input_file_unit, file=file, access='direct', recl=nrecl, iostat=iost, status='old')
    read(unit=input_file_unit, rec=1, iostat=iost) xnrecl, xnspin, xnprec
    if(iost /= 0)then
        write(*,'(2A)')'Error reading the header of the input wavefunction file ', trim(adjustl(file))
        write(*,'(A)')'Stopping now.'     
        if(present(iostat)) iostat = iost
        stop
    endif      
    nprec = nint(xnprec) ! Numerical precision flag
    if(nprec /= 45200 .and. kind_cplx_coeffs /= dp)then
        write(*,'(A)') 'ERROR (read_wavecar module): You seem to be using a double-precision WAVECAR (WAVECAR_double),'
        write(*,'(A)') '                             but you are using a version of BandUP compiled for single-precision.' 
        write(*,'(A)') '                             To generate a double-precision version of BandUP, please change the variable' 
        write(*,'(A)') '                             "kind_cplx_coeffs" from "sp" to "dp" in the file src/math_mod.f90, and recompile the code.'
        write(*,'(A)') 'Stopping now.'
        stop
    endif

    wf%n_spin = nint(xnspin) ! Number of spin channels
    reset_spin = (wf%i_spin < 1 .or. wf%i_spin > wf%n_spin)
    if(wf%i_spin < 1) wf%i_spin = 1
    if(wf%i_spin > wf%n_spin) wf%i_spin = wf%n_spin
    if(reset_spin) write(*,'(A,I1)')'WARNING (read_wavecar): Spin channel reset to ', wf%i_spin
    ispin = wf%i_spin
 
    read(unit=input_file_unit,rec=2, iostat=iost) xnwk,xnband,ecut,(a1(j),j=1,3),(a2(j),j=1,3),(a3(j),j=1,3)
    nwk = nint(xnwk) ! Number of k-points
    if(nwk==0 .or. iost/=0)then
        write(*,'(3A)')'ERROR (read_wavecar): Unexpected file format (file "',trim(adjustl(file)),'")'
        if(nwk==0) write(*,'(A)')"Ifort users: You can try adding/removing the flag '-assume byterecl' when compiling the code."
        write(*,'(A)')'Stopping now.'
        if(present(iostat)) iostat = iost
        stop
    endif
    nband=nint(xnband) ! Number of bands
    wf%n_bands = nband
    wf%Vcell= dabs(dot_product(a1, cross(a2,a3)))
    ! Vectors of the direct lattice.
    wf%A_matrix(1,:) = a1; wf%A_matrix(2,:) = a2; wf%A_matrix(3,:) = a3
    ! Vectors of the reciprocal lattice.
    b1 = twopi*cross(a2,a3)/wf%Vcell; b1mag = norm(b1); b_matrix(1,:) = b1
    b2 = twopi*cross(a3,a1)/wf%Vcell; b2mag = norm(b2); b_matrix(2,:) = b2
    b3 = twopi*cross(a1,a2)/wf%Vcell; b3mag = norm(b3); b_matrix(3,:) = b3
    wf%B_matrix = b_matrix

    wf%encut = ecut

    i_kpt = ikpt
    if((i_kpt > nwk).or.(i_kpt < 1))then
        write(*,'(2(A,I0))')"WARNING (read_wavecar): Invalid choice of i_kpt! &
                             i_kpt should lie between 1 and ", nwk, &
                             ". Reading info for k-point #1 instead of k-point #", &
                             i_kpt
        i_kpt = 1
    endif

    ! Positioning file register...
    irec_start_first_spin=2
    ! ... at the correct spin channel
    irec_start_chosen_spin = irec_start_first_spin + (ispin-1)*nwk*(1+nband) 
    ! ... at the correct k-point
    irec = irec_start_chosen_spin + (i_kpt-1)*(1+nband) + 1 

    allocate(cener(nband)) ! Band energies (as complex numbers)
    allocate(occ(nband)) ! Occupations 
    read(unit=input_file_unit,rec=irec) xnplane,(frac_coords_selected_kpt(i),i=1,3),(cener(iband),occ(iband),iband=1,nband)

close(unit=input_file_unit)

nplane=nint(xnplane)
wf%kpt_frac_coords = frac_coords_selected_kpt
deallocate(wf%band_energies, wf%band_occupations, stat=alloc_stat)
allocate(wf%band_energies(1:nband))
allocate(wf%band_occupations(1:nband))
wf%band_energies(:) = real(cener(:), kind=dp)
wf%band_occupations(:) = occ(:)
deallocate(cener)
deallocate(occ)

! Estimating the maximum integers (nb1,nb2,,nb3) used in the combination
! G(nb1,nb2,,nb3) = sum(nbi*bi; i=1,2,3)
! that defines a general vector G of the reciprocal lattice.
! nbimax == max(nbi), for i =1,2,3. 
phi12 = angle(b1, b2)
phi13 = angle(b1, b3)
phi23 = angle(b2, b3)

! phi123 is the angle between the vector b3 and the plane formed by the vectors b1 and b2
! Thus, phi123 = Pi/2 - angle(b3, cross(b1, b2)) => sin(phi123) = cos[angle(b3, cross(b1, b2))] 
sinphi123 = cos(angle(b3, cross(b1, b2)))
!     b1mag*abs(sin(phi12) gives the proj. of b1 in the direction perpendicular to b2
nb1maxA=(dsqrt(ecut*c)/(b1mag*abs(sin(phi12))))+1 
nb2maxA=(dsqrt(ecut*c)/(b2mag*abs(sin(phi12))))+1
nb3maxA=(dsqrt(ecut*c)/(b3mag*abs(sinphi123)))+1
npmaxA=nint(4.*pi*nb1maxA*nb2maxA*nb3maxA/3.)

! phi123 is now the angle between the vector b2 and the plane formed by the vectors b1 and b3
sinphi123 = cos(angle(b2, cross(b1,b3)))
nb1maxB=(dsqrt(ecut*c)/(b1mag*abs(sin(phi13))))+1
nb2maxB=(dsqrt(ecut*c)/(b2mag*abs(sinphi123)))+1
nb3maxB=(dsqrt(ecut*c)/(b3mag*abs(sin(phi13))))+1
npmaxB=nint(4.*pi*nb1maxB*nb2maxB*nb3maxB/3.)
     
! phi123 is now the angle between the vector b1 and the plane formed by the vectors b2 and b3
sinphi123 = cos(angle(b1, cross(b2,b3))) 
nb1maxC=(dsqrt(ecut*c)/(b1mag*abs(sinphi123)))+1
nb2maxC=(dsqrt(ecut*c)/(b2mag*abs(sin(phi23))))+1
nb3maxC=(dsqrt(ecut*c)/(b3mag*abs(sin(phi23))))+1 
npmaxC=nint(4.*pi*nb1maxC*nb2maxC*nb3maxC/3.)

nb1max=max0(nb1maxA,nb1maxB,nb1maxC)
nb2max=max0(nb2maxA,nb2maxB,nb2maxC)
nb3max=max0(nb3maxA,nb3maxB,nb3maxC)
npmax=min0(npmaxA,npmaxB,npmaxC)

! I'll now count the # of pw coeffs. 
! I do this before storing them and the accosiated RL vecs because, 
! if we're working with spinor wavefunctions, the count will be twice as much as nplane
ncnt=0
do ig3=0,2*nb3max
    ig3p=ig3
    ! Trick to make ig3p range from 0 to nb3max and then from -nb3max to -1
    if (ig3.gt.nb3max) ig3p=ig3-2*nb3max-1  
    do ig2=0,2*nb2max
        ig2p=ig2
        if (ig2.gt.nb2max) ig2p=ig2-2*nb2max-1 
        do ig1=0,2*nb1max
            ig1p=ig1
            if (ig1.gt.nb1max) ig1p=ig1-2*nb1max-1 
            do j=1,3
                sumkg(j)=(frac_coords_selected_kpt(1)+ig1p)*b1(j) + &
                         (frac_coords_selected_kpt(2)+ig2p)*b2(j) + &
                         (frac_coords_selected_kpt(3)+ig3p)*b3(j)
            enddo
            gtot=sqrt(sumkg(1)**2+sumkg(2)**2+sumkg(3)**2)
            etot=gtot**2/c
            if (etot.lt.ecut) then
                ncnt=ncnt+1
            end if
        enddo
    enddo
enddo

wf%n_spinor = 1
if(ncnt /= nplane)then
    if(nplane == 2*ncnt)then
        nplane = ncnt
        wf%n_spinor = 2
    else
        write(*,'(A,I0,A)')'ERROR reading coefficients for wave vector K(',i_kpt,'):' 
        write(*,'(A,I0)')  '    * Number of plane-waves expected for this wave vector: ', ncnt
        write(*,'(A,I0)')  '    * Number plane-waves found in the input file: ', nplane
        write(*,'(A)')     'Cannot continue. Stopping now.'
        stop
    endif
endif
wf%n_pw = nplane
wf%is_spinor = wf%n_spinor==2

! Allocating the matrix that will contain the coordinates (n1,n2,n3) of the vectors 
! G(n1,n2,n3)  = n1*b1 + n2*b2 + n3*b3 of the rec. latt. for which 
! E(K+G)<=ENCUT
deallocate(wf%G_cart, wf%G_frac, stat=alloc_stat)
allocate(wf%G_cart(1:ncnt), wf%G_frac(1:ncnt))
iplane=0
do ig3=0,2*nb3max
    ig3p=ig3
    ! Trick to make ig3p range from 0 to nb3max and then from -nb3max to -1
    if (ig3.gt.nb3max) ig3p=ig3-2*nb3max-1  
    do ig2=0,2*nb2max
        ig2p=ig2
        if (ig2.gt.nb2max) ig2p=ig2-2*nb2max-1 
        do ig1=0,2*nb1max
            ig1p=ig1
            if (ig1.gt.nb1max) ig1p=ig1-2*nb1max-1 
            do j=1,3
                sumkg(j)=(frac_coords_selected_kpt(1)+ig1p)*b1(j) + &
                         (frac_coords_selected_kpt(2)+ig2p)*b2(j) + &
                         (frac_coords_selected_kpt(3)+ig3p)*b3(j)
            enddo
            gtot=sqrt(sumkg(1)**2+sumkg(2)**2+sumkg(3)**2)
            etot=gtot**2/c
            if (etot.lt.ecut) then
                iplane = iplane + 1
                wf%G_cart(iplane)%coord(:) = ig1p*b1(:) + ig2p*b2(:) + ig3p*b3(:)
                wf%G_frac(iplane)%coord(1) = ig1p
                wf%G_frac(iplane)%coord(2) = ig2p
                wf%G_frac(iplane)%coord(3) = ig3p
            end if
        enddo
    enddo
enddo

if(read_coefficients)then
    deallocate(wf%pw_coeffs, stat=alloc_stat)
    allocate(wf%pw_coeffs(1:wf%n_spinor, nplane, nband))
    ! The variables "io_unit_for_thread", "irec_before_loop" and "temp_unit" have been added
    ! to make it possible to use openmp to parallelize (over bands) 
    ! the reading of the WAVECAR file.
    n_threads = 1
    !$ n_threads = omp_get_max_threads() 
    allocate(io_unit_for_thread(1:n_threads))
    io_unit_for_thread(1) = available_io_unit()
    if(n_threads>1)then
        do ithread=2, n_threads
            io_unit_for_thread(ithread) = available_io_unit(io_unit_for_thread(ithread-1)+1,10*n_threads)
        enddo
    endif

    irec_before_loop = irec
    !$omp parallel default(none) &
    !$    private(temp_unit, iost, irec, stream_pos, inner_prod) &
    !$    shared(io_unit_for_thread, file, irec_before_loop, nrecl, wf)
    temp_unit = io_unit_for_thread(1)
    !$ temp_unit = io_unit_for_thread(omp_get_thread_num() + 1)
    open(unit=temp_unit,file=file,access='stream',iostat=iost,status='old')
    !$omp do schedule(guided)
    do iband=1, wf%n_bands
        irec=irec_before_loop+iband
        stream_pos = 1 + (irec - 1) * nrecl
        read(unit=temp_unit, pos=stream_pos) (wf%pw_coeffs(i_wf_comp,:,iband), i_wf_comp=1, wf%n_spinor)
    enddo
    !$omp end do
    close(temp_unit)
    !$omp end parallel
    deallocate(io_unit_for_thread)
endif
if(present(iostat)) iostat = 0

end subroutine read_wavecar


end module read_vasp_wavecar
