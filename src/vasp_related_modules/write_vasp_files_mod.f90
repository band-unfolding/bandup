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
module write_vasp_files
use constants_and_types
use strings
use lists_and_seqs, only: append
use general_io
implicit none
PRIVATE
PUBLIC :: write_crystal_to_file

CONTAINS


subroutine write_crystal_to_file(crystal,output_file,scaling_factor,overwrite,success)
!! Copyright (C) 2014 Paulo V. C. Medeiros
! The output positions will be given in cartesian coordinates
implicit none
type(crystal_3D), intent(in) :: crystal
character(len=*), intent(in) :: output_file
real(kind=dp), intent(in), optional :: scaling_factor
logical, intent(in), optional :: overwrite
logical, intent(out), optional :: success
integer :: unit_num, ivec, icoord, n_species, n_atoms, iatom, i_species, status_io
integer, dimension(:), allocatable :: n_atoms_species, species_in_the_basis
character(len=1), dimension(1:3) :: flags
character(len=3), dimension(:), allocatable :: unique_atomic_symbols
character(len=127) :: aux_string
real(kind=dp) :: a0
logical :: file_exists, stop_if_file_exists


if(present(success)) success = .FALSE.
stop_if_file_exists = .FALSE.
if(present(overwrite))then
    stop_if_file_exists = .not. overwrite
endif
inquire(file=output_file, exist=file_exists)
if(file_exists .and. stop_if_file_exists)then
    write(*,'(3A)')'WARNING: File "',trim(adjustl(output_file)),'" already exists. &
                    The file wil not be written.'
    return
endif

a0 = 1.0d0
if(present(scaling_factor)) a0 = scaling_factor

n_atoms = size(crystal%coords_basis_atoms,dim=1)
call append(item=crystal%integer_types_basis_atoms(1),list=species_in_the_basis)
call append(item=crystal%atomic_symbols_basis_atoms(1),list=unique_atomic_symbols)
if(n_atoms>1)then
    do iatom=2,n_atoms
        if(all(species_in_the_basis /= crystal%integer_types_basis_atoms(iatom)))then
            call append(item=crystal%integer_types_basis_atoms(iatom), &
                        list=species_in_the_basis)
        endif
        if(all(unique_atomic_symbols /= crystal%atomic_symbols_basis_atoms(iatom)))then
            call append(item=crystal%atomic_symbols_basis_atoms(iatom), &
                        list=unique_atomic_symbols)
        endif
    enddo
endif
n_species = size(species_in_the_basis)
allocate(n_atoms_species(1:n_species))
do i_species=1,n_species
    n_atoms_species(i_species) = 0
    do iatom=1,n_atoms
        if(crystal%integer_types_basis_atoms(iatom) == species_in_the_basis(i_species))then
            n_atoms_species(i_species) = n_atoms_species(i_species) + 1
        endif
    enddo
enddo
n_species = count(n_atoms_species>0)

unit_num = available_io_unit()
open(unit=unit_num,file=output_file)
    write(unit_num,"(A,1X,A)")trim(adjustl(crystal%description)), &
                              file_header_BandUP()
    write(unit_num,'(1X,F13.9)')a0
    do ivec=1,3
        write(unit_num,'(3(1X,F13.9))') (crystal%latt_vecs(ivec,icoord)/a0, icoord=1,3)
    enddo
    if(any(len_trim(crystal%atomic_symbols_basis_atoms)>0))then
        write(aux_string,"(120A3,1X)",iostat=status_io) &
              (unique_atomic_symbols(i_species), ' ', &
               i_species=1, n_species)
        call compact(aux_string)
        write(unit_num,"(2X,A)")trim(adjustl(aux_string))
    endif

    write(aux_string,*)(n_atoms_species(i_species),i_species=1,n_species)
    call compact(aux_string)
    write(unit_num,"(2X,A)")trim(adjustl(aux_string))

    if(any(crystal%unconstrained_dof_basis_atoms) .eqv. .FALSE.) write(unit_num,'(A)')'Selective Dynamics'
    write(unit_num,'(A)') 'Cartesian' 
    ! Writing atom coordinates to the file
    do i_species=1, n_species
        if(any(crystal%unconstrained_dof_basis_atoms) .eqv. .FALSE.)then
            do iatom=1,n_atoms
                if(crystal%integer_types_basis_atoms(iatom) == species_in_the_basis(i_species))then
                    do icoord=1,3
                        if(crystal%unconstrained_dof_basis_atoms(iatom,icoord))then
                            flags(icoord) = 'T'
                        else
                            flags(icoord) = 'F'
                        endif
                    enddo
                    write(unit_num,'(3(1X,F13.9),3(2X,A))')crystal%coords_basis_atoms(iatom,:)/a0, flags 
                endif
            enddo
        else
            do iatom=1,n_atoms
                if(crystal%integer_types_basis_atoms(iatom) == species_in_the_basis(i_species))then
                    write(unit_num,'(3(1X,F13.9))')crystal%coords_basis_atoms(iatom,:)/a0
                endif
            enddo
        endif
    enddo
close(unit_num)

if(file_exists .and. .not. stop_if_file_exists)then
    write(*,'(3A)')'       WARNING: The file "',trim(adjustl(output_file)),'" has been overwritten.'
endif
if(present(success)) success = .TRUE.

end subroutine write_crystal_to_file



end module write_vasp_files
