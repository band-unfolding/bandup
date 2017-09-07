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
! MODULE: crystals 
!
!> @author
!> Paulo V C Medeiros, Linköping University
!
! DESCRIPTION: 
!> Defines some derived types and routines to work with crystals.
!==============================================================================

module crystals
use constants_and_types
use math, only: coords_cart_vec_in_new_basis, norm, same_vector, cross, &
                triple_product, inverse_of_3x3_matrix
use lists_and_seqs, only: append
implicit none
PRIVATE
PUBLIC :: get_rec_latt, vec_in_latt, reduce_point_to_bz, create_crystal, &
          check_if_pc_and_SC_are_commensurate

CONTAINS

subroutine check_if_pc_and_SC_are_commensurate(&
               commensurate, M, crystal_pc, crystal_SC,tol &
           )
!! Copyright (C) 2014 Paulo V. C. Medeiros
! Calculates the matrix M so that  
! A[i] = sum(M_{ij}*a[j]; j=1,2,3), where
! 'A' and 'a' are the real space latt vecs of the SC and pc, respectively. 
! The subroutine then uses M to check if the PC and SC are commensurate.
implicit none
logical, intent(out) :: commensurate
real(kind=dp), dimension(1:3,1:3), intent(out) :: M
type(crystal_3D), intent(in) :: crystal_pc, crystal_SC
real(kind=dp), dimension(1:3,1:3) :: b_matrix_pc, B_matrix_SC
real(kind=dp), intent(in), optional :: tol
integer, dimension(1:3,1:3) :: int_M
real(kind=dp), dimension(1:3,1:3) :: residue_M
real(kind=dp) :: max_residue

    b_matrix_pc = crystal_pc%rec_latt_vecs
    B_matrix_SC = crystal_SC%rec_latt_vecs
    commensurate = .TRUE.
    max_residue = default_tol_for_int_commens_test
    if(present(tol))then
        max_residue = dabs(tol)
    endif
    M = transpose(matmul(b_matrix_pc,inverse_of_3x3_matrix(B_matrix_SC)))
    int_M = nint(M)
    residue_M = dabs(M - int_M)
    if(any(residue_M > max_residue))then
        commensurate = .FALSE.
    endif

end subroutine check_if_pc_and_SC_are_commensurate

subroutine create_crystal(&
               crystal, description, latt_vecs, coords_basis_atoms, &
               atomic_symbols_basis_atoms, unconstrained_dof_basis_atoms &
           )
!! Copyright (C) 2014 Paulo V. C. Medeiros
implicit none
type(crystal_3D), intent(out) :: crystal
character(len=*), intent(in), optional :: description
real(kind=dp), dimension(1:3,1:3), intent(in) :: latt_vecs
! coords_basis_atoms = coords_basis_atoms(:,1:3)
real(kind=dp), dimension(:,:), intent(in), optional :: coords_basis_atoms 
character(len=3), dimension(:), optional :: atomic_symbols_basis_atoms
! dimension{unconstrained_dof_basis_atoms} = (:,1:3)
logical, dimension(:,:), optional :: unconstrained_dof_basis_atoms 
integer :: iatom, isymb
real(kind=dp), dimension(1:3) :: default_positions
logical, dimension(1:1,1:3) :: default_unconstrained_dof
real(kind=dp), dimension(:,:), allocatable :: atom_pos ! dimension(:,1:3)
character(len=3), dimension(:), allocatable :: symbols, unique_symbols
logical, dimension(:,:), allocatable :: dof_basis_atoms ! dimension(:,1:3)
type(crystal_3D), allocatable :: aux_crystal  
    ! Defaults of optional arguments 
    default_positions(:) = 0.0_dp * (latt_vecs(1,:) + &
                                     latt_vecs(2,:) + &
                                     latt_vecs(3,:))
    default_unconstrained_dof(1,:) = (/.TRUE.,.TRUE.,.TRUE./) 

    ! 'Cleaning' crystal
    allocate(aux_crystal)
    crystal = aux_crystal

    ! Handling optional arguments
    if(present(coords_basis_atoms))then
        allocate(atom_pos(1:size(coords_basis_atoms,dim=1), &
                          1:size(coords_basis_atoms,dim=2)))
        atom_pos = coords_basis_atoms
    else
        allocate(atom_pos(1:1,1:3))
        atom_pos(1,:) = default_positions
    endif
    if(present(unconstrained_dof_basis_atoms))then
        allocate(dof_basis_atoms(1:size(unconstrained_dof_basis_atoms,dim=1), &
                                 1:size(unconstrained_dof_basis_atoms,dim=2)))
        dof_basis_atoms = unconstrained_dof_basis_atoms
    else
        allocate(dof_basis_atoms(1:1,1:3))
        dof_basis_atoms = default_unconstrained_dof
    endif
    if(present(atomic_symbols_basis_atoms))then
        allocate(symbols(1:size(atomic_symbols_basis_atoms)))
        symbols = atomic_symbols_basis_atoms
    else
        allocate(symbols(1:size(atom_pos,dim=1)))
        symbols(1) = ''
    endif
    if(present(description))then
        crystal%description = trim(adjustl(description))
    else
        crystal%description = ''
    endif
    ! End of handling optional arguments
    
    ! Lattice
    crystal%latt_vecs = latt_vecs
    crystal%vol = abs(triple_product(&
                          latt_vecs(1,:), latt_vecs(2,:), latt_vecs(3,:)&
                      )&
                  )
    call get_rec_latt(latt_vecs, crystal%rec_latt_vecs, crystal%rec_latt_vol)
    ! Basis
    allocate(crystal%coords_basis_atoms(1:size(atom_pos,dim=1), &
                                        1:size(atom_pos,dim=2)))
    crystal%coords_basis_atoms = atom_pos
    allocate(crystal%fractional_coords_basis_atoms(1:size(atom_pos,dim=1), &
                                                   1:size(atom_pos,dim=2)))
    do iatom=1, size(atom_pos,dim=1)
        crystal%fractional_coords_basis_atoms(iatom,:) = &
            coords_cart_vec_in_new_basis(atom_pos(iatom,:), latt_vecs(:,:))
    enddo
    ! Atomic symbols
    allocate(crystal%atomic_symbols_basis_atoms(1:size(symbols)))
    crystal%atomic_symbols_basis_atoms = symbols
    call append(item=symbols(1), list=unique_symbols)
    if(size(symbols)>1)then
        do iatom=2, size(symbols)
            if(all(unique_symbols /= symbols(iatom)))then
                call append(item=symbols(iatom), list=unique_symbols)
            endif
        enddo
    endif
    allocate(crystal%integer_types_basis_atoms(1:size(symbols)))
    do iatom=1, size(symbols)
        do isymb=1, size(unique_symbols)
            if(symbols(iatom) == unique_symbols(isymb))then
                crystal%integer_types_basis_atoms(iatom) = isymb
                exit
            endif
        enddo
    enddo

    ! Degrees of freedom
    allocate(&
        crystal%unconstrained_dof_basis_atoms(1:size(dof_basis_atoms,dim=1), &
                                              1:size(dof_basis_atoms,dim=2)) &
    )
    crystal%unconstrained_dof_basis_atoms = dof_basis_atoms


end subroutine create_crystal


subroutine get_rec_latt(latt,rec_latt,rec_latt_vol)
implicit none
real(kind=dp), dimension(1:3,1:3), intent(in) :: latt
real(kind=dp), dimension(1:3,1:3), intent(out) :: rec_latt
real(kind=dp), intent(out), optional :: rec_latt_vol
real(kind=dp) :: vpc

    vpc = dabs(dot_product(latt(1,:), cross(latt(2,:),latt(3,:))))
    rec_latt(1,:) = twopi*cross(latt(2,:),latt(3,:))/vpc;
    rec_latt(2,:) = twopi*cross(latt(3,:),latt(1,:))/vpc;
    rec_latt(3,:) = twopi*cross(latt(1,:),latt(2,:))/vpc;

    if(present(rec_latt_vol))then
        rec_latt_vol = (8.0*(Pi**3))/vpc
    endif

end subroutine get_rec_latt

function vec_in_latt(vec, latt, tolerance) result(rtn)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
!! Not fully tested, but works for BandUP
!! Checks if the vector "vec" belongs to the Bravais lattice spanned by the
!! vectors latt(i,:), i=1,2,3.
implicit none
logical :: rtn
real(kind=dp), dimension(1:3), intent(in) :: vec
real(kind=dp), dimension(1:3,1:3), intent(in) :: latt
real(kind=dp), intent(in), optional :: tolerance
real(kind=dp), dimension(1:3) :: reduced_vec,frac_coords,reduced_frac_coords,g
real(kind=dp) :: tol
integer :: i,ig1,ig2,ig3,ig1p,ig2p,ig3p,nb1max,nb2max,nb3max

tol = default_tol_for_vec_equality
if(present(tolerance))then
     tol = dabs(tolerance)
endif

rtn = .FALSE.
frac_coords(:) = coords_cart_vec_in_new_basis(cart_vec=vec, new_basis=latt)
reduced_frac_coords(:) = mod(frac_coords(:),1.0_dp)
reduced_vec(:) = 0.0_dp
do i=1,3
    reduced_vec(:) = reduced_vec(:) + reduced_frac_coords(i)*latt(i,:)
enddo

nb1max = 1
nb2max = 1
nb3max = 1
do ig3=0,2*nb3max
    ig3p=ig3
    ! Trick to make the indices ig3p vary from 0 to nb3max and then 
    ! from -nb3max to -1
    if (ig3.gt.nb3max) ig3p=ig3-2*nb3max-1  
    do ig2=0,2*nb2max
        ig2p=ig2
        if (ig2.gt.nb2max) ig2p=ig2-2*nb2max-1
        do ig1=0,2*nb1max
            ig1p=ig1
            if (ig1.gt.nb1max) ig1p=ig1-2*nb1max-1
            g(:) = ig1p*latt(1,:) + ig2p*latt(2,:) + ig3p*latt(3,:)
            if(same_vector(reduced_vec,g,tol))then
                rtn = .TRUE.
                return
            endif
        enddo
    enddo
enddo
return

end function vec_in_latt


function point_is_in_bz(point,rec_latt,origin_point,tolerance,verbose) &
result(rtn)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
!! Checks if the point "point" belongs to the 1st Brillouin zone of a Bravais
!! latt for which the reciprocal lettice vectors are "rec_latt(i,:)", i=1,2,3
!! The opt variable "origin_point" is to be used if you want the mesh to be
!! centered in other point than (0,0,0).
!! Not fully tested, but works for BandUP
implicit none
logical :: rtn
real(kind=dp), dimension(1:3), intent(in) :: point
real(kind=dp), dimension(1:3,1:3), intent(in) :: rec_latt
real(kind=dp), dimension(1:3), optional, intent(in) :: origin_point
real(kind=dp), optional, intent(in) :: tolerance
logical, intent(in), optional :: verbose
real(kind=dp), dimension(1:3) :: origin, g
real(kind=dp) :: tol
integer :: ig1, ig2, ig3
logical :: print_stuff

rtn = .TRUE.
if(present(origin_point))then
    origin(:) = origin_point(:)
else
    origin(:) = 0.0_dp
endif

tol = abs(default_tol_for_vec_equality)
if(present(tolerance))then
    tol = abs(tolerance)
endif

print_stuff = .FALSE.
if(present(verbose))then
    print_stuff = verbose
endif

do ig3=-1,1
    do ig2=-1,1
        do ig1=-1,1
            if(ig1==0 .and. ig2==0 .and. ig3==0) cycle
            g(:) = origin(:) + real(ig1,kind=dp)*rec_latt(1,:) + &
                               real(ig2,kind=dp)*rec_latt(2,:) + & 
                               real(ig3,kind=dp)*rec_latt(3,:)
            if(norm(point(:) - origin(:)) > norm(point(:) - g(:)))then
                if(norm(point(:) - origin(:)) - norm(point(:) - g(:)) > &
                   min(1E-8_dp, 0.01_dp * default_tol_for_vec_equality))then
                    rtn = .FALSE.
                    return
                else
                    if(print_stuff)then
                        write(*,'(A)')&
                            'WARNING (point_is_in_bz): Tolerance applied &
                            when checking whether point lies inside 1BZ. ,&
                            This should be OK, but do check your results.'
                    endif
                endif
            endif

        enddo
    enddo
enddo

end function point_is_in_bz


recursive subroutine reduce_point_to_bz(point, rec_latt,point_reduced_to_bz, &
                                        max_index_lin_comb_RL_vecs)
!! Copyright (C) 2013, 2014 Paulo V. C. Medeiros
!! Takes a point "point" and returns another point "point_reduced_to_bz" which
!! is equivalent to "point", but lies inside the 1st BZ of the Bravais lattice
!! with reciprocal lattice vectors "rec_latt(i,:)", i=1,2,3.
!! Not fully tested, but works for BandUP
implicit none
real(kind=dp), dimension(1:3), intent(in) :: point
real(kind=dp), dimension(1:3,1:3), intent(in) :: rec_latt
! point_reduced_to_bz is returned in cartesian coords!
real(kind=dp), dimension(1:3), intent(out) :: point_reduced_to_bz 
integer, intent(in), optional :: max_index_lin_comb_RL_vecs
logical :: point_has_been_reduced
real(kind=dp), dimension(1:3) :: frac_coords_point, nearest_G, ref_G, &
                                 trial_nearest_G
integer :: ig1, ig2, ig3, igmax

    igmax = 0
    if(present(max_index_lin_comb_RL_vecs))then
        igmax = max_index_lin_comb_RL_vecs
    endif

    frac_coords_point = coords_cart_vec_in_new_basis(&
                            cart_vec=point, new_basis=rec_latt &
                        )
    ! Finding the rec. latt. point "nearest_G" that is 
    ! closest to the input point "point"
    ! NB.: ref_G will be in cartesian coords
    ref_G = matmul(real(nint(frac_coords_point), kind=dp),rec_latt) 
    nearest_G = ref_G
    do ig3=-igmax,igmax
        do ig2=-igmax,igmax
            do ig1=-igmax,igmax
                trial_nearest_G = ref_G + &
                                  ig1*rec_latt(1,:) + &
                                  ig2*rec_latt(2,:) + &
                                  ig3*rec_latt(3,:)
                if(norm(point - nearest_G) > norm(point - trial_nearest_G))then
                    nearest_G = trial_nearest_G
                endif
            enddo 
        enddo 
    enddo
    ! Reducing the point, now that we have the nearest rec. latt. point
    point_reduced_to_bz(:) = point(:) - nearest_G(:)

    ! Double-checking if the point has really been reduced 
    ! The routine will raise an error and stop otherwise
    point_has_been_reduced = &
        point_is_in_bz(point_reduced_to_bz,rec_latt,verbose=.TRUE.)
    if(.not. point_has_been_reduced)then
        if(igmax>=10)then !  Recursion limit
            write(*,'(A)')&
                'ERROR (reduce_point_to_bz): Failed to reduce &
                k-point to the 1st Brillouin zone.'
            write(*,'(A,I0,A)')&
                '                            A total of ',igmax+1,' layer(s) &
                of neighboring rec. latt. points &
                were checked for nearest neighbors.'
            write(*,'(A)')&
                '                            The nearest rec. latt. point was &
                not found.'
            write(*,'(3(A,f0.5),A)')&
                '                            * Cartesian coordinates of the &
                k-point:  (', point(1),', ',point(2),', ',point(3),').'
            write(*,'(3(A,f0.5),A)')&
                '                            * Fractional coordinates of the &
                k-point: (', frac_coords_point(1),', ',frac_coords_point(2), &
                ', ', frac_coords_point(3),').'
            write(*,'(A)')'Stopping now.'
            stop
        else
            ! If the previous scan fails, go for an extra layer
            if(igmax>=10)then !  Warn about use of too many layers of kpoints
                write(*,'(A)')&
                    'WARNING (reduce_point_to_bz): Quite large number of &
                    layers of neighboring rec. latt. points. &
                    This should be OK, but do check your results.'
            endif
            call reduce_point_to_bz(point,rec_latt,point_reduced_to_bz,igmax+1)
        endif
    endif

end subroutine reduce_point_to_bz

end module crystals
