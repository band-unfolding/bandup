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
module read_vasp_files
use constants_and_types
use strings
use crystals, only: create_crystal
use general_io
implicit none
PRIVATE
PUBLIC :: get_crystal_from_file

CONTAINS

function vasp_version_of_file(file_name) result(version)
!! Copyright (C) 2014 Paulo V. C. Medeiros
implicit none
integer :: version
character(len=*), intent(in) :: file_name
character(len=360) :: read_string
character(len=80) :: description
integer :: i,j,unit_num, ios
real(kind=dp) :: scaling_factor
real(kind=dp), dimension(1:3,1:3) :: latt_vecs

    unit_num=available_io_unit()
    open(unit=unit_num,file=file_name,status='old',action='read',blank='null',pad='yes')
        read(unit_num,"(A)")description
        read(unit_num,*)scaling_factor
        do i=1,3
            read(unit_num,*)(latt_vecs(i,j),j=1,3)
        enddo
        read(unit_num,"(A)", iostat=ios)read_string
    close(unit_num)
    call compact(read_string)
    if(ios /= 0 .or. trim(adjustl(read_string)) == '')then
        version = 0 ! Means that the file only contains the header
    else
        if(verify(read_string,'0123456789 ')/=0)then
            version = 5
        else
            version = 4
        endif
    endif

end function vasp_version_of_file



subroutine get_crystal_from_file(crystal,input_file,stop_if_file_not_found,success)
!! Copyright (C) 2014 Paulo V. C. Medeiros
! The output positions will be given in cartesian coordinates
implicit none
type(crystal_3D), intent(out) :: crystal
character(len=*), intent(in) :: input_file
logical, intent(in), optional :: stop_if_file_not_found
logical, intent(out), optional :: success
integer :: vasp_version, unit_num, ivec, icoord, n_species, n_atoms, iatom, &
           i_species, ierr, i_atom_species
integer, dimension(:), allocatable :: n_atoms_species
character(len=9) :: coordinates_type
character(len=18) :: dynamics_type
character(len=127) :: description
character(len=360) :: read_string ! aux string to help read data
character(len=1), dimension(1:3) :: flags
character(len=3), dimension(:), allocatable :: atomic_symbols_species, atomic_symbols_basis_atoms
real(kind=dp) :: scaling_factor
real(kind=dp), dimension(1:3,1:3) :: latt_vecs
real(kind=dp), dimension(:,:), allocatable :: r, positions
logical, dimension(:,:), allocatable :: unconstrained_dof_basis_atoms
logical :: file_exists, stop_if_non_existing_file


if(present(success)) success = .FALSE.
stop_if_non_existing_file = .TRUE.
if(present(stop_if_file_not_found))then
    stop_if_non_existing_file = stop_if_file_not_found
endif
inquire(file=input_file, exist=file_exists)
if(.not. file_exists)then
    if(stop_if_non_existing_file)then
        write(*,'(3A)')'ERROR: File "',trim(adjustl(input_file)),'" not found. Stopping now.'
        stop
    else
        write(*,'(3A)')'WARNING: File "',trim(adjustl(input_file)),'" not found. The file has not been read.'
        return
    endif
endif

vasp_version = vasp_version_of_file(input_file)
unit_num = available_io_unit()
open(unit=unit_num,file=input_file,status='old',action='read',blank='null',pad='yes')
    read(unit_num,"(A)")description
    read(unit_num,*)scaling_factor
    do ivec=1,3
        read(unit_num,*) (latt_vecs(ivec,icoord), icoord=1,3)
        latt_vecs(ivec,:) = scaling_factor*latt_vecs(ivec,:)
    enddo
    if(vasp_version == 0)then
        ! Creating an artificial crystal
        call create_crystal(crystal, latt_vecs=latt_vecs)
    else
        if(vasp_version==5)then
            read(unit_num,"(A)")read_string
            call split(read_string, atomic_symbols_species)
        endif
        read(unit_num,"(A)")read_string
        call  split(read_string,n_atoms_species) 
        if(vasp_version==4)then
            allocate(atomic_symbols_species(1:size(n_atoms_species)))
            do i_species=1,size(n_atoms_species)
                write(read_string,*)i_species
                atomic_symbols_species(i_species) = trim(adjustl(read_string))
            enddo 
        endif
        read(unit_num,"(A)")read_string
        read_string = trim(adjustl(read_string))
        if(upper_case(read_string(1:1))=='S')then
            dynamics_type='Selective Dynamics'
            read(unit_num,"(A)")coordinates_type
        else
            dynamics_type=''
            read(read_string,"(A)")coordinates_type
        endif
        coordinates_type=trim(adjustl(coordinates_type))
        select case(upper_case(coordinates_type(1:1)))
            case('C','K')
                coordinates_type='Cartesian'
            case('D')
                coordinates_type='Direct'
            case default
                write(*,'(A)')'WARNING: Positions assumed to have been given in direct coordenates.'
                coordinates_type='Direct'
        end select

        n_species=size(n_atoms_species)
        n_atoms = sum(n_atoms_species)

        deallocate(atomic_symbols_basis_atoms, stat=ierr)
        allocate(atomic_symbols_basis_atoms(1:n_atoms))

        iatom = 0
        do i_species=1,n_species
            do i_atom_species=1,n_atoms_species(i_species)
                iatom = iatom + 1
                atomic_symbols_basis_atoms(iatom) = atomic_symbols_species(i_species)
            enddo
        enddo

        deallocate(r, stat=ierr)
        allocate(r(1:n_atoms,1:3))
        deallocate(positions, stat=ierr)
        allocate(positions(1:n_atoms,1:3),stat=ierr)
        deallocate(unconstrained_dof_basis_atoms,stat=ierr)
        allocate(unconstrained_dof_basis_atoms(1:n_atoms,1:3))
        unconstrained_dof_basis_atoms(:,:) = .TRUE. ! Atoms are considered to be unconstrained by default

        ! Reading atom coordinates from the file
        if(trim(adjustl(dynamics_type))=='Selective Dynamics')then
            do iatom=1,n_atoms
                read(unit_num,*)r(iatom,:), flags(:)
                do icoord=1,3
                    flags(icoord) = upper_case(flags(icoord))
                enddo
                unconstrained_dof_basis_atoms(iatom,:) = (flags(:) /= 'F')
            enddo
        else
            do iatom=1,n_atoms
                read(unit_num,*)r(iatom,:)
            enddo
        endif

        ! The output positions will be given in cartesian coordinates
        if(trim(adjustl(coordinates_type))=='Cartesian')then
            positions(:,:) = scaling_factor * r(:,:)
        else
            do iatom=1,n_atoms
                positions(iatom,:) = matmul(r(iatom,:),latt_vecs(:,:)) 
            enddo
        endif
        call create_crystal(crystal, description, latt_vecs, positions, &
                            atomic_symbols_basis_atoms, unconstrained_dof_basis_atoms)
    endif
close(unit_num)

if(present(success)) success = .TRUE.

end subroutine get_crystal_from_file



end module read_vasp_files
