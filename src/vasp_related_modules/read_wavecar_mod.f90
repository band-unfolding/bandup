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
!! has been granted by the Copyright holders. Thanks!

!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module read_wavecar
! ===================
!  Contains a function "read_wavecar", which returns information about any
!  selected k-point, skipping the others. This is very important for the
!  performance of the BandUP code (or any other code that would require similar functionality).
!  See the definition of the function for more information on the arguments.
!  I've only tested this module on linux and using the ifort compiler.
!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$************************* List of major changes from the WaveTrans code ****************************************
!!$*
!!$*       1: Turned the code into a module (read_wavecar), which can be used in the  program.
!!$*          1.a: The routine "read_from_wavefunc_file" returns information about any selected k-point, skipping the others.
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
!!$*
!!$*       2: openmp is now used to parallelize (over bands) the reading of the WAVECAR file
!!$*          2.a: Introduced the function "available_io_unit"
!!$*               Each thread associates the WAVECAR file to a different IO unit in the "open" statements.
!!$*
!!$*       3: Added support for spinor wavefunctions
!!$*          
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

module read_wavecar
use math
use general_io, only : available_io_unit, WF_file
!$ use omp_lib
implicit none
PRIVATE
PUBLIC :: read_from_wavefunc_file, vasp_2m_over_hbar_sqrd 
!!$* Comments adapted from WaveTrans: The value of vasp_2m_over_hbar_sqrd has been adjusted in final decimal places to agree with VASP value
!!$*                                  The routine "read_from_wavefunc_file" checks for discrepancy of any results between this and VASP values
real(kind=dp), parameter :: vasp_2m_over_hbar_sqrd = 0.262465831 ! c = 2m/hbar**2 in units of 1/eV Ang^2
CONTAINS


subroutine read_from_wavefunc_file(i_selected_kpt, kpt_frac_coords, n_spin_channels, spin_channel, n_plane_waves, &
                                   energies_bands, occupations, i_allg_in_file, coeff, &
                                   file_size_in_bytes, rec_latt_b_matrix, latt_vecs, total_nkpts, ENCUT, &
                                   elapsed_time, add_elapsed_time_to, iostat)
implicit none
!!$* Input and output variables
integer, optional, intent(in) :: i_selected_kpt ! Index of the selected SCBZ k-point
real(kind=dp), dimension(1:3), optional, intent(out) :: kpt_frac_coords ! Coords of the SCBZ k-point in the [b1, b2, b3] basis
integer, intent(out), optional :: n_plane_waves ! Number of SCBZ K-wave-vectors that differ from the selected kpt by vectors of the SC rec. latt.
integer(8), intent(out), optional :: file_size_in_bytes
integer, optional, intent(out) :: n_spin_channels
integer, optional, intent(in) :: spin_channel
real(kind=dp), dimension(:), allocatable, optional, intent(out) :: energies_bands, occupations  ! E = E(K; iband)
integer, dimension(:,:), allocatable, optional, intent(out) :: i_allg_in_file
complex(kind=sp), dimension(:,:,:), allocatable, optional, intent(out) :: coeff ! Coefficients, coeff(ikpt, iband)
! Matrices of the direct/reciprocal lattice vectors. b_matrix(i,:) = bi
real(kind=dp), dimension(1:3,1:3), intent(out), optional :: rec_latt_b_matrix, latt_vecs 
integer, optional, intent(out) :: total_nkpts ! Total number of k-points in the input file
real(kind=dp), optional, intent(out) :: ENCUT
real(kind=dp), optional, intent(inout) :: add_elapsed_time_to, elapsed_time
integer, optional, intent(out) :: iostat
!!$* Parameters
real(kind=dp), parameter :: c = vasp_2m_over_hbar_sqrd ! 2m/hbar**2 in units of 1/eV Ang^2
!!$* Variables used internally
real(kind=dp), dimension(1:3) :: frac_coords_selected_kpt 
integer :: nplane ! Number of SCBZ K-wave-vectors that differ from the selected kpt by vectors of the SC rec. latt.
complex(kind=dp), allocatable :: cener(:)  ! Re{cener(iband)} = E(K; iband)
integer, dimension(:,:), allocatable :: igall
real(kind=dp), dimension(:), allocatable :: occ
real(kind=dp), dimension(1:3) :: a1,a2,a3,b1,b2,b3,sumkg
real(kind=dp), dimension(1:3,1:3) :: b_matrix ! Matrix of the reciprocal lattice vectors. b_matrix(i,:) = bi
real(kind=dp) :: vcell, b1mag, b2mag, b3mag, ecut, gtot, etot, phi12, phi13, phi23, sinphi123, &
                 xnwk, xnband, xnrecl, xnspin, xnprec, xnplane, stime, ftime
integer :: i, j, iost, irec, ig1, ig2, ig3, ig1p, ig2p, ig3p, iplane, ispin,         &
           irec_start_first_spin, irec_start_chosen_spin, nrecl, nwk, nspin, nprec, nband,     &
           nb1maxA, nb2maxA, nb3maxA, nb1maxB, nb2maxB, nb3maxB, nb1maxC, nb2maxC, nb3maxC,    &
           npmaxA, npmaxB, npmaxC, nb1max, nb2max, nb3max, npmax, iband, ncnt,                 &
           input_file_unit, irec_before_loop,temp_unit
integer, dimension(:), allocatable :: available_io_units
complex(kind=sp), dimension(:), allocatable :: band_coeff
complex(kind=sp), dimension(:,:,:), allocatable :: aux_coeff ! Aux coefficients when using spinor wavefunctions
logical :: file_exists, assume_spinor_wavecar

! Start
stime = time()

if(present(iostat))then
    iostat = -1
endif

if(.not.present(i_selected_kpt))then
    if(present(n_plane_waves)   .or. &
       present(kpt_frac_coords) .or. &
       present(energies_bands)  .or. &
       present(occupations)     .or. &
       present(i_allg_in_file)  .or. &
       present(coeff)          )then
        write(*,'(A)')'ERROR (read_wavecar_mod): Not able to retrieve information from the wavefunctions file. Please specify a k-point index.'
        write(*,'(A)')'                          Stopping now.'
        stop
    endif
endif



input_file_unit = available_io_unit()
nrecl=24 ! Trial record length (recl = length of each record in a file connected for direct access)
open(unit=input_file_unit,file=WF_file,access='direct',recl=nrecl,iostat=iost,status='old')
    if (iost.ne.0) then
        write(*,'(2A)')'Error opening the input wavefunction file ', trim(adjustl(WF_file))
        inquire(file=WF_file, exist=file_exists)
        if(file_exists)then
            write(*,'(A)')'    * The file exists, but could not be read.'
        else
            write(*,'(A)')'    * It seems that the file is missing.'
        endif 
        write(*,'(A)')'Stopping now.'     
        stop
    else
        if(present(file_size_in_bytes)) inquire(unit=input_file_unit, size=file_size_in_bytes)
    endif       
    read(unit=input_file_unit,rec=1) xnrecl,xnspin,xnprec
close(unit=input_file_unit)

nrecl=nint(xnrecl) ! Correct record length
nprec=nint(xnprec) ! Numerical precision flag
if(nprec.eq.45210) then
    write(*,'(A)') 'Error: WAVECAR_double requires complex*16 (double precision). Stopping now.'
    stop
endif

nspin=nint(xnspin) ! Number of spin channels
if(present(n_spin_channels))then
    n_spin_channels = nspin
endif

ispin = 1
if(present(spin_channel))then
    if((spin_channel>=1).and.(spin_channel<=nspin))then
        ispin = spin_channel
    else
        write(*,'(A)')'WARNING (read_wavecar_mod): The wavefunction contains only 1 spin channel.'
        write(*,'(A)')'                            Reading information on the this spin.'
    endif
endif


input_file_unit = available_io_unit()
open(unit=input_file_unit,file=WF_file,access='direct',recl=nrecl,iostat=iost,status='old')
    if (iost.ne.0) then
        write(*,'(2A)')'Error opening the input wavefunction file ', trim(adjustl(WF_file))
        write(*,'(A)')'    * The file exists, but could not be read.'
        write(*,'(A)')'Stopping now.'
        close(unit=input_file_unit)
        stop
    endif       
    read(unit=input_file_unit,rec=2) xnwk,xnband,ecut,(a1(j),j=1,3),(a2(j),j=1,3),(a3(j),j=1,3)
    nwk=nint(xnwk) ! Number of k-points
    if(nwk==0)then
        write(*,'(A)')'ERROR (read_wavecar module): no KPTS could be read from the WAVECAR file.'
        write(*,'(A)')"For ifort users: Have you used the flag '-assume byterecl' when compiling the read_wavecar module? You should use it."
        write(*,'(A)')'Stopping now.'
        stop
    endif
    nband=nint(xnband) ! Number of bands
!   Compute properties of the reciprocal lattice
!   Volume of the cell
    Vcell= dabs(dot_product(a1, cross(a2,a3)))
!   Vectors of the reciprocal lattice.
    b1 = twopi*cross(a2,a3)/Vcell; b1mag = norm(b1); b_matrix(1,:) = b1
    b2 = twopi*cross(a3,a1)/Vcell; b2mag = norm(b2); b_matrix(2,:) = b2
    b3 = twopi*cross(a1,a2)/Vcell; b3mag = norm(b3); b_matrix(3,:) = b3

    if(present(latt_vecs))then
        latt_vecs(1,:) = a1
        latt_vecs(2,:) = a2
        latt_vecs(3,:) = a3
    endif
    if(present(ENCUT))then
        ENCUT = ecut
    endif
    if(present(total_nkpts))then
        total_nkpts = nwk
    endif
    if(present(rec_latt_b_matrix))then
        rec_latt_b_matrix = b_matrix
    endif
    if(present(n_plane_waves)   .or. &
       present(kpt_frac_coords) .or. &
       present(energies_bands)  .or. &
       present(occupations)     .or. &
       present(i_allg_in_file)  .or. &
       present(coeff)          )then
        if((i_selected_kpt > nwk).or.(i_selected_kpt < 1))then
            write(*,'(2(A,I0),A)') 'ERROR: Cannot parse i_selected_kpt = ',i_selected_kpt,'. i_selected_kpt should be in between 1 and ', &
                                    nwk,' (the number fo K-points in the wavefunctions file).'
            write(*,'(A)') 'Stopping now.'
            stop
        endif
    else
        close(unit=input_file_unit)
        if(present(iostat))then
           iostat = 0
        endif
        if(present(elapsed_time))then
            ftime = time()
            elapsed_time = ftime - stime
        endif
        if(present(add_elapsed_time_to))then
            ftime = time()
            add_elapsed_time_to = add_elapsed_time_to  + (ftime - stime)
        endif
        return
    endif

    irec_start_first_spin=2
    irec_start_chosen_spin = irec_start_first_spin + (ispin-1)*nwk*(1+nband) ! Positioning the register at the correct spin channel
    irec = irec_start_chosen_spin + (i_selected_kpt-1)*(1+nband) + 1 ! Positioning the register at the correct k-point

    allocate(occ(nband)) ! Occupations 
    allocate(cener(nband)) ! Band energies (as a complex number)
    read(unit=input_file_unit,rec=irec) xnplane,(frac_coords_selected_kpt(i),i=1,3),(cener(iband),occ(iband),iband=1,nband)

    nplane=nint(xnplane)
    if(present(n_plane_waves))then
        n_plane_waves = nplane
    endif
    if(present(kpt_frac_coords))then
        kpt_frac_coords = frac_coords_selected_kpt
    endif
    if(present(energies_bands))then
        allocate(energies_bands(1:nband))
        energies_bands(:) = real(cener(:), kind=dp)
    endif
    if(present(occupations))then
        allocate(occupations(1:nband))
        occupations(:) = occ(:)
    endif
 

    if((.not.present(i_allg_in_file)).and.(.not.present(coeff)))then
        close(unit=input_file_unit)
        if(present(iostat))then
           iostat = 0
        endif
        if(present(elapsed_time))then
            ftime = time()
            elapsed_time = ftime - stime
        endif
        if(present(add_elapsed_time_to))then
            ftime = time()
            add_elapsed_time_to = add_elapsed_time_to  + (ftime - stime)
        endif
        return
    endif
!   Estimating the maximum integers (nb1,nb2,,nb3) used in the combination
!   G(nb1,nb2,,nb3) = sum(nbi*bi; i=1,2,3)
!   that defines a general vector G of the reciprocal lattice.
!   nbimax == max(nbi), for i =1,2,3. 
!
    phi12 = angle(b1, b2)
    phi13 = angle(b1, b3)
    phi23 = angle(b2, b3)

!   phi123 is the angle between the vector b3 and the plane formed by the vectors b1 and b2
!   Thus, phi123 = Pi/2 - angle(b3, cross(b1, b2)) => sin(phi123) = cos[angle(b3, cross(b1, b2))] 
    sinphi123 = cos(angle(b3, cross(b1, b2)))
    nb1maxA=(dsqrt(ecut*c)/(b1mag*abs(sin(phi12))))+1 ! b1mag*abs(sin(phi12) gives the proj. of b1 in the direction perpendicular to b2
    nb2maxA=(dsqrt(ecut*c)/(b2mag*abs(sin(phi12))))+1
    nb3maxA=(dsqrt(ecut*c)/(b3mag*abs(sinphi123)))+1
    npmaxA=nint(4.*pi*nb1maxA*nb2maxA*nb3maxA/3.)

!   phi123 is now the angle between the vector b2 and the plane formed by the vectors b1 and b3
    sinphi123 = cos(angle(b2, cross(b1,b3)))
    nb1maxB=(dsqrt(ecut*c)/(b1mag*abs(sin(phi13))))+1
    nb2maxB=(dsqrt(ecut*c)/(b2mag*abs(sinphi123)))+1
    nb3maxB=(dsqrt(ecut*c)/(b3mag*abs(sin(phi13))))+1
    npmaxB=nint(4.*pi*nb1maxB*nb2maxB*nb3maxB/3.)
     
!   phi123 is now the angle between the vector b1 and the plane formed by the vectors b2 and b3
    sinphi123 = cos(angle(b1, cross(b2,b3))) 
    nb1maxC=(dsqrt(ecut*c)/(b1mag*abs(sinphi123)))+1
    nb2maxC=(dsqrt(ecut*c)/(b2mag*abs(sin(phi23))))+1
    nb3maxC=(dsqrt(ecut*c)/(b3mag*abs(sin(phi23))))+1 
    npmaxC=nint(4.*pi*nb1maxC*nb2maxC*nb3maxC/3.)

    nb1max=max0(nb1maxA,nb1maxB,nb1maxC)
    nb2max=max0(nb2maxA,nb2maxB,nb2maxC)
    nb3max=max0(nb3maxA,nb3maxB,nb3maxC)
    npmax=min0(npmaxA,npmaxB,npmaxC)


    ! Getting the indexes of the planewave coefficients
    ! I count them before storing them because, if we're working with spinor
    ! wavefunctions, the count will be twice as much as nplane
    ncnt=0
    do ig3=0,2*nb3max
        ig3p=ig3
        if (ig3.gt.nb3max) ig3p=ig3-2*nb3max-1  ! Trick to make ig3p range from 0 to nb3max and then from -nb3max to -1
        do ig2=0,2*nb2max
            ig2p=ig2
            if (ig2.gt.nb2max) ig2p=ig2-2*nb2max-1 
            do ig1=0,2*nb1max
                ig1p=ig1
                if (ig1.gt.nb1max) ig1p=ig1-2*nb1max-1 
                do j=1,3
                    sumkg(j)=(frac_coords_selected_kpt(1)+ig1p)*b1(j) + (frac_coords_selected_kpt(2)+ig2p)*b2(j) + (frac_coords_selected_kpt(3)+ig3p)*b3(j)
                enddo
                gtot=sqrt(sumkg(1)**2+sumkg(2)**2+sumkg(3)**2)
                etot=gtot**2/c
                if (etot.lt.ecut) then
                    ncnt=ncnt+1
                end if
            enddo
        enddo
    enddo

    assume_spinor_wavecar = .FALSE.
    if(ncnt /= nplane)then
        if(nplane == 2*ncnt)then
            assume_spinor_wavecar = .TRUE.
        else
            write(*,'(A,I0,A)')'ERROR reading coefficients for wave vector K(',i_selected_kpt,'):' 
            write(*,'(A,I0)')  '    * Number of plane-waves expected for this wave vector: ', ncnt
            write(*,'(A,I0)')  '    * Number plane-waves found in the input file: ', nplane
            write(*,'(A)')     'Cannot continue. Stopping now.'
            stop
        endif
    endif

!   Allocating the matrix that will contain the coordinates (n1,n2,n3) of the vectors 
!   G(n1,n2,n3)  = n1*b1 + n2*b2 + n3*b3 of the rec. latt. for which 
!   E(K+G)<=ENCUT
    allocate (igall(3,ncnt))
    iplane=0
    do ig3=0,2*nb3max
        ig3p=ig3
        if (ig3.gt.nb3max) ig3p=ig3-2*nb3max-1  ! Trick to make ig3p range from 0 to nb3max and then from -nb3max to -1
        do ig2=0,2*nb2max
            ig2p=ig2
            if (ig2.gt.nb2max) ig2p=ig2-2*nb2max-1 
            do ig1=0,2*nb1max
                ig1p=ig1
                if (ig1.gt.nb1max) ig1p=ig1-2*nb1max-1 
                do j=1,3
                    sumkg(j)=(frac_coords_selected_kpt(1)+ig1p)*b1(j) + (frac_coords_selected_kpt(2)+ig2p)*b2(j) + (frac_coords_selected_kpt(3)+ig3p)*b3(j)
                enddo
                gtot=sqrt(sumkg(1)**2+sumkg(2)**2+sumkg(3)**2)
                etot=gtot**2/c
                if (etot.lt.ecut) then
                    iplane = iplane + 1
                    igall(1, iplane)=ig1p
                    igall(2, iplane)=ig2p
                    igall(3, iplane)=ig3p
                end if
            enddo
        enddo
    enddo

  
    if(present(i_allg_in_file))then
        allocate (i_allg_in_file(3,ncnt))
        i_allg_in_file(:,:) = igall(:,1:ncnt)
        deallocate(igall)
    endif

    if(present(coeff))then
        allocate (coeff(1:1, nplane, nband)) ! I'll extend it later to a two-component object if working with a spinor WF

        ! The variables "available_io_units", "irec_before_loop" and "temp_unit" have been added
        ! to make it possible to use openmp to parallelize (over bands) 
        ! the reading of the WAVECAR file.
        allocate(available_io_units(1:nband))
        available_io_units(1) = available_io_unit()
        if(nband>1)then
            do iband=2,nband
                available_io_units(iband) = available_io_unit(available_io_units(iband-1)+1,10*nband)
            enddo
        endif
        irec_before_loop = irec
        !$omp parallel do default(none) schedule(guided) &
        !$    private(iband,irec,iplane,temp_unit,iost,band_coeff) &
        !$    shared(nband,irec_before_loop,available_io_units,WF_file,coeff,nplane,nrecl)
        do iband=1,nband
            irec=irec_before_loop+iband
            temp_unit = available_io_units(iband)
            open(unit=temp_unit,file=WF_file,access='direct',recl=nrecl,iostat=iost,status='old')
                read(unit=temp_unit,rec=irec) (coeff(1,iplane,iband), iplane=1,nplane)
            close(temp_unit)
        enddo
        deallocate(available_io_units)
    endif

close(unit=input_file_unit)


if(present(coeff) .and. assume_spinor_wavecar)then
    allocate(aux_coeff(1:2, 1:ncnt, 1:nband))
    aux_coeff(1,:,:) = coeff(1, 1:ncnt, :)
    aux_coeff(2,:,:) = coeff(1, ncnt+1:nplane, :)
    deallocate(coeff)
    allocate(coeff(1:2, 1:ncnt, 1:nband)) ! Now using a two-component matrix of coeffs
    coeff(:,:,:) = aux_coeff(:,:,:)
    deallocate(aux_coeff)
endif


if(present(elapsed_time))then
    ftime = time()
    elapsed_time = ftime - stime
endif
if(present(add_elapsed_time_to))then
    ftime = time()
    add_elapsed_time_to = add_elapsed_time_to  + (ftime - stime)
endif

if(present(iostat))then
   iostat = 0  ! The whole operation ended just fine.
endif

return

end subroutine read_from_wavefunc_file

end module read_wavecar
