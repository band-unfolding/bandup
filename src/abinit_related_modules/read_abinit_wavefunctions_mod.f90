!! Copyright (C) 2015 Paulo V. C. Medeiros
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
! MODULE: read_abinit_wavefunctions
!
!> @author
!> Paulo V C Medeiros, Linköping University
!
! DESCRIPTION: 
!> Provides routines to read info from ABINIT's wavefunctions (WFK files)
!! into BandUP's wavefunction object.
!==============================================================================
module read_abinit_wavefunctions
use general_io
use constants_and_types
use units, only: to_angstrom, to_ev
use math, only: cross, norm
implicit none
PRIVATE
PUBLIC :: read_abinit_wfk_file, read_abinit_wfk_header

CONTAINS


subroutine read_abinit_wfk_header(unit, n_sppol, n_kpt, A_matrix, kpts, iostat)
! Copyright (C) 2015 Paulo V. C. Medeiros
! Written based on the information found at
! www.abinit.org/documentation/helpfiles/for-v7.8/users/abinit_help.html#header
implicit none
! Mandatory input and output variables
integer,intent(in) :: unit
! Optional input and output variables
integer, intent(out), optional :: n_sppol, n_kpt, iostat
real(kind=dp), dimension(1:3,1:3), intent(out), optional :: A_matrix
type(vec3d), dimension(:), allocatable, intent(out), optional :: kpts
! Local variables
character(len=6) :: codvsn
character(len=132) :: title
integer :: headform, fform, bantot, date, intxc, ixc, natom, nkpt, nspden, &
           nspinor, nsppol, nsym, npsp, ntypat, occopt, pertcase, usepaw, &
           usewvl, pspso, pspdat, pspcod, pspxc, lmn_size, ipsp, ikpt, &
           alloc_stat
integer, dimension(1:3) :: ngfft
integer, dimension(:), allocatable :: istwfk, nband, npwarr, so_psp, symafm, &
                                      typat, nrhoijsel
integer, dimension(:,:), allocatable :: rhoijselect
integer, dimension(:,:,:), allocatable :: symrel
! ABINIT uses double precision to write the binary files
double precision :: ecut, ecutdg, ecutsm, ecut_eff, stmbias, tphysel, tsmear, &
                    znuclpsp, zionpsp
double precision, dimension(1:3) :: qptn
double precision, dimension(1:3, 1:3) :: rprimd
double precision, dimension(:), allocatable :: occ, znucltypat, wtk
double precision, dimension(:,:), allocatable :: kpt, tnons, xred, rhoij


    if(present(iostat)) iostat = -1

    read(unit) codvsn, headform, fform    
    read(unit) bantot, date, intxc, ixc, natom, ngfft(1:3), &
               nkpt, nspden, nspinor, nsppol, nsym, npsp, ntypat, occopt, &
               pertcase, usepaw, ecut, ecutdg, ecutsm, ecut_eff, qptn(1:3), &
               rprimd(1:3,1:3), stmbias, tphysel, tsmear, usewvl

    if(present(n_sppol)) n_sppol = nsppol
    if(present(n_kpt)) n_kpt = nkpt
    if(present(A_matrix)) A_matrix = rprimd


    allocate(istwfk(1:nkpt))
    allocate(nband(1:nkpt*nsppol))
    allocate(npwarr(1:nkpt))
    allocate(so_psp(1:npsp))
    allocate(symafm(1:nsym))
    allocate(symrel(1:3,1:3,1:nsym))
    allocate(typat(1:natom))
    allocate(nrhoijsel(1:nspden))
    allocate(rhoijselect(1:maxval(nrhoijsel),1:nspden))

    allocate(kpt(1:3,nkpt))
    allocate(occ(1:bantot))
    allocate(tnons(1:3,1:nsym))
    allocate(znucltypat(1:ntypat))
    allocate(wtk(1:nkpt))

    allocate(xred(1:3,1:natom))
    allocate(rhoij(1:maxval(nrhoijsel),nspden))

    read(unit) istwfk(1:nkpt), nband(1:nkpt*nsppol), npwarr(1:nkpt), &
               so_psp(1:npsp),symafm(1:nsym),symrel(1:3,1:3,1:nsym), &
               typat(1:natom), kpt(1:3,1:nkpt), occ(1:bantot), &
               tnons(1:3,1:nsym), znucltypat(1:ntypat), wtk(1:nkpt)


    if(present(kpts))then
        deallocate(kpts, stat=alloc_stat)
        allocate(kpts(1:nkpt))
        do ikpt=1, nkpt
            kpts(ikpt)%coord(:) = kpt(:,ikpt)
        enddo
    endif

    do ipsp=1,npsp
        read (unit) title, znuclpsp, zionpsp, pspso, pspdat, pspcod, pspxc, &
                    lmn_size
    end do
    read(unit) 
    if(usepaw==1)then
        read(unit)
        read(unit)
    end if

    if(present(iostat)) iostat = 0

end subroutine read_abinit_wfk_header


subroutine reconstruct_gamma_point_wf(wf)
! Copyright (C) 2015 Paulo V. C. Medeiros
! ABINIT uses the property that psi[k] - {psi[-k]}* at gamma sometimes
! to reduce the size of the wavefunction at gamma
! This routine reconstructs a full wavefunction
implicit none
type(pw_wavefunction), intent(inout) :: wf
integer :: npw, nspinor, nband, npw_gamma_aux, ipw, alloc_stat
type(pw_wavefunction) :: aux_wf

    if(all(wf%kpt_frac_coords(:) == 0))then
        npw = wf%n_pw
        nspinor = size(wf%pw_coeffs(:, :, :), dim=1)
        nband = size(wf%pw_coeffs(:, :, :), dim=3)

        npw_gamma_aux = 2*npw - 1
        aux_wf = wf
        ! Reconstructing coordinates of translations
        deallocate(aux_wf%G_frac, aux_wf%G_cart, stat=alloc_stat)
        allocate(aux_wf%G_frac(1:npw_gamma_aux),aux_wf%G_cart(1:npw_gamma_aux))
        do ipw=1, npw
            aux_wf%G_frac(ipw)%coord(:) = wf%G_frac(ipw)%coord(:)
            aux_wf%G_cart(ipw)%coord(:) = wf%G_cart(ipw)%coord(:)
            if(ipw==1) cycle
            aux_wf%G_frac(npw+ipw-1)%coord(:) = -wf%G_frac(ipw)%coord(:)
            aux_wf%G_cart(npw+ipw-1)%coord(:) = -wf%G_cart(ipw)%coord(:)
        enddo
        deallocate(wf%G_frac, wf%G_cart, stat=alloc_stat)
        allocate(wf%G_frac(1:npw_gamma_aux), wf%G_cart(1:npw_gamma_aux))
        wf%G_frac = aux_wf%G_frac
        wf%G_cart = aux_wf%G_cart

        ! Reconstructing plane wave coefficients
        deallocate(aux_wf%pw_coeffs, stat=alloc_stat)
        allocate(aux_wf%pw_coeffs(1:nspinor, 1:npw_gamma_aux, 1:nband))
        aux_wf%pw_coeffs(:, 1:npw, :) = wf%pw_coeffs(:, :, :)
        aux_wf%pw_coeffs(:, npw+1:, :) = conjg(wf%pw_coeffs(:, 2:, :))
        deallocate(wf%pw_coeffs)
        allocate(wf%pw_coeffs(1:nspinor, 1:npw_gamma_aux, 1:nband))
        wf%pw_coeffs = aux_wf%pw_coeffs

        ! Updating number of plane waves
        wf%n_pw = npw_gamma_aux
    endif

end subroutine reconstruct_gamma_point_wf


subroutine read_abinit_wfk_file(wf, file, ikpt, read_coeffs, iostat)
! Copyright (C) 2015 Paulo V. C. Medeiros
! Written based on the information found at
! www.abinit.org/documentation/helpfiles/for-v7.0/users/
! abinit_help.html#wavefctfile
implicit none
! Mandatory input and output variables
type(pw_wavefunction), intent(inout) :: wf
character(len=*), intent(in) :: file
integer, intent(in) :: ikpt
! Optional input and output variables
integer, optional, intent(out) :: iostat
logical, optional, intent(in) :: read_coeffs
! Local variables
type(vec3d), dimension(:), allocatable :: all_kpts
integer :: i_spin, n_kpt, i_kpt, chosen_kpt, nband, iband, npw, &
           ipw, nspinor, alloc_stat, io_unit
integer, dimension(:,:), allocatable :: kg ! plane wave reduced coordinates
real(kind=dp) :: b1mag, b2mag, b3mag, inner_prod
! Direct and rec. latt. vecs.
real(kind=dp), dimension(1:3) :: a1, a2, a3, b1, b2, b3 
! ABINIT uses double precision to write the binary file
real(kind=dp), dimension(:), allocatable :: eigen, occ
complex(kind=dp), dimension(:,:), allocatable :: cg
logical :: read_coefficients


read_coefficients = .TRUE. ! Reading the coeffs by default
if(present(read_coeffs)) read_coefficients = read_coeffs
if(present(iostat)) iostat = -1

io_unit = available_io_unit() 
open(unit=io_unit, file=file,form='unformatted', recl=1)
    call read_abinit_wfk_header(&
             unit=io_unit, &
             n_sppol=wf%n_spin, n_kpt=n_kpt, A_matrix=wf%A_matrix, &
             kpts=all_kpts &
         )

    ! Latt. vecs.
    wf%A_matrix = bohr * wf%A_matrix
    a1 = wf%A_matrix(1,:)
    a2 = wf%A_matrix(2,:)
    a3 = wf%A_matrix(3,:)
    wf%Vcell= dabs(dot_product(a1, cross(a2,a3)))
    ! Vectors of the reciprocal lattice.
    b1 = twopi*cross(a2,a3)/wf%Vcell; b1mag = norm(b1); wf%B_matrix(1,:) = b1
    b2 = twopi*cross(a3,a1)/wf%Vcell; b2mag = norm(b2); wf%B_matrix(2,:) = b2
    b3 = twopi*cross(a1,a2)/wf%Vcell; b3mag = norm(b3); wf%B_matrix(3,:) = b3

    chosen_kpt = ikpt
    if((ikpt > n_kpt).or.(ikpt < 1))then
        write(*,'(2(A,I0))')&
            "WARNING (read_abinit_wfk_file): Invalid choice of ikpt! &
            ikpt should lie between 1 and ", n_kpt, &
            ". Reading info for k-point #1 instead of k-point #", &
            ikpt
        chosen_kpt = 1
    endif

    wf%kpt_frac_coords(:)  = all_kpts(chosen_kpt)%coord(:)
    wf%kpt_cart_coords(:)  = wf%kpt_frac_coords(1) * b1(:) + &
                             wf%kpt_frac_coords(2) * b2(:) + &
                             wf%kpt_frac_coords(3) * b3(:)

    do i_spin=1, wf%n_spin
        do i_kpt=1, n_kpt
            read(io_unit) npw, nspinor, nband
            if((i_spin /= wf%i_spin) .or. (i_kpt /= chosen_kpt))then
                ! Skipping unwanted kpts and/or spin channel
                read(io_unit)
                read(io_unit)
                do iband=1, nband
                    read(io_unit)
                enddo
            else
                wf%n_pw = npw
                wf%n_spinor = nspinor
                wf%n_bands = nband
                wf%is_spinor = wf%n_spinor==2
                allocate(kg(1:3,1:npw))
                allocate(eigen(1:nband), occ(1:nband))

                read(io_unit) kg(:,:)
                read(io_unit) eigen(:), occ(:)

                deallocate(&
                    wf%band_energies, wf%band_occupations, stat=alloc_stat &
                )
                allocate(&
                    wf%band_energies(1:nband), wf%band_occupations(1:nband) &
                )
                wf%band_energies(:) = Hartree * eigen(:)
                wf%band_occupations(:) = occ(:)

                deallocate(wf%G_frac, stat=alloc_stat)
                deallocate(wf%G_cart, stat=alloc_stat)
                allocate(wf%G_frac(1:npw), wf%G_cart(1:npw))
                do ipw=1, npw
                    wf%G_frac(ipw)%coord(:) = kg(:,ipw)
                    wf%G_cart(ipw)%coord(:) = kg(1,ipw) * b1(:) + &
                                              kg(2,ipw) * b2(:) + &
                                              kg(3,ipw) * b3(:)
                enddo
                deallocate(kg, eigen, occ)

                if(read_coefficients)then
                    allocate(cg(1:npw*nspinor, 1:nband))
                    do iband=1,nband
                        ! wf coeffs for each band at this kpt
                        read(io_unit) (cg(ipw, iband), ipw=1, npw*nspinor) 
                         
                        ! Chasing huge coeffs
                        ! I don't know why they appear, but by taking them out
                        ! one can normally still get a good qualitative view
                        inner_prod=abs(dot_product(cg(:, iband), cg(:,iband)))
                        if(inner_prod > huge(1.0_kind_cplx_coeffs))then
                            write(*,'(A, I0, A, D12.4, A)') &
                                '    WARNING: Squared norm of SC wavefunction &
                                for band #', iband, &
                                ' is too big (> ', &
                                huge(1.0_kind_cplx_coeffs),').'
                            write(*, '(A)')&
                                '             The coeffs have been zeroed.'
                            write(*, '(A)')&
                                '             Please check your results.'
                            cg(:, iband) = 0.0_dp
                        endif
                    enddo

                    ! Passing coeffs to BandUP's pw_wavefunction object
                    deallocate(wf%pw_coeffs, stat=alloc_stat)
                    allocate(wf%pw_coeffs(1:nspinor, 1:npw, 1:nband))
                    wf%pw_coeffs = 0.0_dp
                    wf%pw_coeffs(1, :, :) = cg(1:npw,:)
                    if(wf%is_spinor) wf%pw_coeffs(2, :, :) = cg(npw+1:,:)
                    ! Reconstructing gamma-point wavefunctions
                    if(all(wf%kpt_frac_coords(:) == 0) .and. &
                       .not. wf%is_spinor)then
                        call reconstruct_gamma_point_wf(wf)
                    endif

                    deallocate(cg)
                else
                    ! Skipping coeffs
                    do iband=1, nband
                        read(io_unit)
                    enddo
                endif

            endif
        enddo
    enddo
close(io_unit)

if(present(iostat)) iostat = -1

end subroutine read_abinit_wfk_file

end module read_abinit_wavefunctions

