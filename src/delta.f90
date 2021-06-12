! Copyright (C) 2020- Nakib Haider Protik <nakib.haider.protik@gmail.com>
! This file is part of elphbolt <https://github.com/nakib/elphbolt>.
!
! elphbolt is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! elphbolt is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with elphbolt. If not, see <http://www.gnu.org/licenses/>.

module delta
  !! Module containing the procedures related to delta function evaulation.

  use params, only: dp, k8
  use misc, only: exit_with_message, print_message, demux_vector, mux_vector, &
       binsearch, sort
  
  implicit none

  private
  public form_tetrahedra_3d, fill_tetrahedra_3d, delta_fn_tetra, &
       form_triangles, fill_triangles, delta_fn_triang
  
contains

  pure real(dp) function delta_fn_tetra(e, ik, ib, mesh, tetramap, tetracount, tetra_evals)
    !! Calculate delta function using the tetraheron method.
    !!
    !! e Sample energy
    !! ik Wave vector index
    !! ib Band index
    !! mesh Wave vector grid
    !! tetramap Wave vector to (tetrahedron, vertex) mapping
    !! tetracount Number of tetrahedra in which a wave vector belongs
    !! tetra_evals Tetrahedra populated with the eigenvalues

    real(dp), intent(in) :: e
    integer(k8), intent(in) :: ik, ib
    integer(k8), intent(in) :: mesh(3), tetramap(:,:,:), tetracount(:)
    real(dp), intent(in) :: tetra_evals(:,:,:)

    !Local variables
    integer(k8) :: iv, it, itk, num, numtetra
    logical :: c1, c2, c3
    real(dp) :: e1, e2, e3, e4, e1e, e2e, e3e, e4e, &
         e21, e31, e41, e32, e42, e43, tmp ! eji \equiv ej - ei

    tmp = 0.0_dp
    delta_fn_tetra = 0.0_dp

    !Total number of tetrahedra in the system
    numtetra = product(mesh)*6
    
    !Grab number of tetrahedra in which wave vector belongs
    num = tetracount(ik)

    do itk = 1, num !Run over tetrahedra
       it = tetramap(1, ik, itk) !Grab tetrahedron
       iv = tetramap(2, ik, itk) !Grab vertex

       !Grab vertex energies
       e1 = tetra_evals(it, ib, 1)
       e2 = tetra_evals(it, ib, 2)
       e3 = tetra_evals(it, ib, 3)
       e4 = tetra_evals(it, ib, 4)

       !Define the energy differences
       e1e = e1 - e
       e2e = e2 - e
       e3e = e3 - e
       e4e = e4 - e
       e21 = e2 - e1
       e31 = e3 - e1
       e41 = e4 - e1
       e32 = e3 - e2
       e42 = e4 - e2
       e43 = e4 - e3

       !Evaluate the three cases
       c1 = e1 <= e .and. e <= e2
       c2 = e2 <= e .and. e <= e3
       c3 = e3 <= e .and. e <= e4

       if(.not. (e < e1 .or. e > e4)) then
          !Evaluate the expressions for the three cases
          select case(iv)
          case(1)
             if(c1) then
                tmp = (e2e/e21 + e3e/e31 + e4e/e41)*(e1e**2)/e41/e31/e21

                if(e1 == e2) then
                   tmp = 0.0_dp
                end if
             else if(c2) then
                tmp = -0.5_dp*(e3e/(e31**2)*(e3e*e2e/e42/e32 + e4e*e1e/e41/e42 + e3e*e1e/e32/e41) &
                     + e4e/(e41**2)*(e4e*e1e/e42/e31 + e4e*e2e/e42/e32 + e3e*e1e/e31/e32))

                if(e2 == e3) then
                   tmp = -0.5_dp*(e4e*e1e/e41/e42 + e1e/e41 &
                        + e4e/(e41**2)*(e4e*e1e/e42/e31 + e4e/e42 + e1e/e31))
                end if
             else if(c3) then
                tmp = (e4e**3)/(e41**2)/e42/e43

                if(e3 == e4) then
                   tmp = (e4e**2)/(e41**2)/e42
                end if
             end if
          case(2)
             if(c1) then
                tmp = -(e1e**3)/(e21**2)/e31/e41

                if(e1 == e2) then
                   tmp = 0.0_dp
                end if
             else if(c2) then
                tmp = -0.5_dp*(e3e/(e32**2)*(e3e*e2e/e42/e31 + e4e*e2e/e42/e41 + e3e*e1e/e31/e41) &
                     + e4e/(e42**2)*(e3e*e2e/e32/e31 + e4e*e1e/e41/e31 + e4e*e2e/e32/e41))

                if(e2 == e3) then
                   tmp = -0.5_dp*(0.0 + e4e/e42/e41 + 0.0 &
                        + e4e/(e42**2)*(0.0 + e4e*e1e/e41/e31 + 1.0))
                end if
             else if(c3) then
                tmp = (e4e**3)/e41/(e42**2)/e43

                if(e3 == e4) then
                   tmp = 0.0_dp
                end if
             end if
          case(3)
             if(c1) then
                tmp = -(e1e**3)/e21/(e31**2)/e41

                if(e1 == e2) then
                   tmp = 0.0_dp
                end if
             else if(c2) then
                tmp = 0.5_dp*(e2e/(e32**2)*(e3e*e2e/e42/e31 + e4e*e2e/e42/e41 + e3e*e1e/e31/e41) &
                     + e1e/(e31**2)*(e3e*e2e/e42/e32 + e4e*e1e/e41/e42 + e3e*e1e/e32/e41))

                if(e2 == e3) then
                   tmp = 0.5_dp*(0.0 + e4e/e42/e41 + e1e/e31/e41 &
                        + e1e/(e31**2)*(0.0 + e4e*e1e/e41/e42 + e1e/e41))
                end if
             else if(c3) then
                tmp = (e4e**3)/e41/e42/(e43**2)

                if(e3 == e4) then
                   tmp = 0.0_dp
                end if
             end if
          case(4)
             if(c1) then
                tmp = -(e1e**3)/e21/e31/(e41**2)
                if(e1 == e2) then
                   tmp = 0.0_dp
                end if
             else if(c2) then
                tmp = 0.5_dp*(e2e/(e42**2)*(e3e*e2e/e32/e31 + e4e*e1e/e41/e31 + e4e*e2e/e32/e41) &
                     + e1e/(e41**2)*(e4e*e1e/e42/e31 + e4e*e2e/e42/e32 + e3e*e1e/e31/e32))

                if(e2 == e3) then
                   tmp = 0.5_dp*(0.0 &
                        + e1e/(e41**2)*(e4e*e1e/e42/e31 + e4e/e42 + e1e/e31))
                end if
             else if(c3) then
                tmp = -(e3e/e43 + e2e/e42 + e1e/e41)*(e4e**2)/e41/e42/e43

                if(e3 == e4) then
                   tmp = 0.0_dp
                end if
             end if
          end select

          if ((e1 == e2) .and. (e1 == e3) .and. (e1 == e4) .and. (e == e1)) then
             tmp = 0.25_dp
          end if

          delta_fn_tetra = delta_fn_tetra + tmp
       end if ! .not. (e <= e1 .or. e >= e4)
    end do !itk

    if(delta_fn_tetra < 1.0e-12_dp) delta_fn_tetra = 0.0_dp
    
    !Normalize with the total number of tetrahedra
    delta_fn_tetra = delta_fn_tetra/numtetra
  end function delta_fn_tetra
  
  subroutine form_tetrahedra_3d(nk, mesh, tetra, tetracount, tetramap, &
       blocks, indexlist)
    !! Form all the tetrahedra of a 3d FBZ mesh.
    !!
    !! nk Number of points in the list of FBZ wave vectors
    !! mesh Wave vector grid
    !! tetra List of the tetrahedra vertices
    !! tetracount Number of tetrahedra in which a wave vector belongs
    !! tetramap Wave vector to (tetrahedron, vertex) mapping
    !! blocks Is the FBZ wave vector list full or energy restricted?
    !! indexlist List of muxed indices of the FBZ wave vectors

    integer(k8), intent(in) :: nk, mesh(3)
    integer(k8), intent(out), allocatable :: tetra(:,:), tetracount(:), tetramap(:,:,:)
    logical, intent(in) :: blocks
    integer(k8), optional, intent(in) :: indexlist(:)

    !Local variables
    integer(k8) :: ik, i, j, k, ijk(3), ii, jj, kk, tk, tl, aux, count
    integer(k8) :: ip1, jp1, kp1, n1, n2, n3, tmp
    integer(k8) :: tetra_vertices_labels(6, 4)
    integer(k8) :: scvol_vertices(8, 3) ! subcell volume vertices

    n1 = mesh(1)
    n2 = mesh(2)
    n3 = mesh(3)

    !Label of the vertices of the tetrahedra for a given subcell
    tetra_vertices_labels = reshape((/ &
         1, 2, 3, 6, &
         1, 3, 5, 6, &
         3, 5, 6, 7, &
         3, 6, 7, 8, &
         3, 4, 6, 8, &
         2, 3, 4, 6 /), &
         shape(tetra_vertices_labels), order = (/2, 1/))

    !Allocate tetrahedra related variables
    allocate(tetra(6*nk, 4), tetracount(nk), tetramap(2, nk, 24))

    tetra(:,:) = 0
    tetracount(:) = 0
    tetramap(:,:,:) = 0
    count = 1 !tetrahedron counter
    
    do ik = 1, nk !Run over all wave vectors in FBZ
       if(blocks) then !For energy window restricted FBZ
          call demux_vector(indexlist(ik), ijk, mesh, 1_k8)
       else !For unrestristed FBZ
          call demux_vector(ik, ijk, mesh, 1_k8)
       end if
       i = ijk(1)
       j = ijk(2)
       k = ijk(3)

       !Apply periodic boundary condition
       if (i == n1) then
          ip1 = 1
       else
          ip1 = i + 1
       end if

       if (j == n2) then
          jp1 = 1
       else
          jp1 = j + 1
       end if

       if (k == n3) then
          kp1 = 1
       else
          kp1 = k + 1
       end if

       !For each subcell save the vertices
       scvol_vertices = reshape((/ &
            i,   j,   k,   &
            ip1, j,   k,   &
            i,   jp1, k,   &
            ip1, jp1, k,   &
            i,   j,   kp1, &
            ip1, j,   kp1, &
            i,   jp1, kp1, &
            ip1, jp1, kp1 /), &
            shape(scvol_vertices), order = (/2, 1/))

       do tk = 1, 6 !Run over 6 tetrahedra
          do tl = 1, 4 !Run over the labels of the vertices that
             !make up each tetrahedron
             aux = tetra_vertices_labels(tk, tl)
             ii = scvol_vertices(aux,1)
             jj = scvol_vertices(aux,2)
             kk = scvol_vertices(aux,3)
             aux = mux_vector((/ii, jj, kk/), mesh, 1_k8)
             tmp = aux !Guaranteed to be > 0
             if(blocks) then
                !Which point in indexlist does aux correspond to?
                call binsearch(indexlist, aux, tmp) !tmp < 0 if search fails.
             end if
             tetra(count, tl) = tmp

             if(tmp > 0) then
                !Save the mapping of a wave vector index to a (tetrahedron, vertex)
                tetracount(tmp) = tetracount(tmp) + 1
                tetramap(1, tmp, tetracount(tmp)) = count
                tetramap(2, tmp, tetracount(tmp)) = tl
             end if
          end do
          count = count + 1
       end do
    end do
  end subroutine form_tetrahedra_3d

  subroutine fill_tetrahedra_3d(tetra, evals, tetra_evals)
    !! Populate the (sorted along the vertices) eigenvalues on all the vertices of the tetrahedra
    !!
    !! tetra List of the tetrahedra vertices
    !! evals List of eigenvalues 
    !! tetra_evals Tetrahedra populated with the eigenvalues

    integer(k8), intent(in) :: tetra(:,:)
    real(dp), intent(in) :: evals(:,:)
    real(dp), allocatable, intent(out) :: tetra_evals(:,:,:)

    !Local variables
    integer(kind=4) :: iv, it, ib, numbands, aux, numtetra

    numtetra = size(tetra(:, 1))
    numbands = size(evals(1, :))

    allocate(tetra_evals(numtetra, numbands, 4))

    !Note: Eigenvalues outside the transport active window is taken to be zero.
    !      As such, close to the transport window boundary, this method is inaccurate.
    !      A large enough transport window must be chosen to obtain accurate transport coefficients.
    tetra_evals(:,:,:) = 0.0_dp
    
    do it = 1, numtetra !Run over tetrahedra
       do ib = 1, numbands !Run over bands
          do iv = 1, 4 !Run over vertices
             aux = tetra(it, iv)
             if(aux > 0) then !Only eigenvalues inside transport active region
                tetra_evals(it, ib, iv) = evals(aux, ib)
             end if
          end do
          call sort(tetra_evals(it, ib, :))
       end do
    end do
  end subroutine fill_tetrahedra_3d

  !Triangular method
  subroutine form_triangles(nk, mesh, triang, triangcount, triangmap, &
       blocks, indexlist)
    !! Form all the triangles of a 3d FBZ mesh for each z component.
    !!
    !! nk Number of points in the list of FBZ wave vectors
    !! mesh Wave vector grid
    !! triang List of the triangle vertices
    !! triangcount Number of triangles in which a wave vector belongs
    !! triangmap Wave vector to (triangle, vertex) mapping
    !! blocks Is the FBZ wave vector list full or energy restricted?
    !! indexlist List of muxed indices of the FBZ wave vectors

    integer(k8), intent(in) :: nk, mesh(3)
    integer(k8), intent(out), allocatable :: triang(:,:), triangcount(:), triangmap(:,:,:)
    logical, intent(in) :: blocks
    integer(k8), optional, intent(in) :: indexlist(:)

    !Local variables
    integer(k8) :: ik, i, j, k, ijk(3), ii, jj, kk, tk, tl, aux, count
    integer(k8) :: ip1, jp1, kp1, n1, n2, n3, tmp
    integer(k8) :: triang_vertices_labels(2, 3)
    integer(k8) :: scvol_vertices(4, 3) !subcell vertices

    n1 = mesh(1)
    n2 = mesh(2)
    n3 = mesh(3)

    !Label of the vertices of the triangles for a given subcell
    triang_vertices_labels = reshape((/ &
         1, 2, 3,&
         1, 3, 4 /), &
         shape(triang_vertices_labels), order = (/2, 1/))

    !Allocate triangles related variables
    allocate(triang(2*nk, 3), triangcount(nk), triangmap(2, nk, 6))

    triang(:,:) = 0
    triangcount(:) = 0
    triangmap(:,:,:) = 0
    count = 1 !triangles counter

    do ik = 1, nk !Run over all wave vectors in FBZ
       if(blocks) then !For energy window restricted FBZ
          call demux_vector(indexlist(ik), ijk, mesh, 1_k8)
       else !For unrestristed FBZ
          call demux_vector(ik, ijk, mesh, 1_k8)
       end if
       i = ijk(1)
       j = ijk(2)
       k = ijk(3)

       !Apply periodic boundary condition
       if (i == n1) then
          ip1 = 1
       else
          ip1 = i + 1
       end if

       if (j == n2) then
          jp1 = 1
       else
          jp1 = j + 1
       end if

       !For each subcell save the vertices
       scvol_vertices = reshape((/ &
            i,   j,   k,    &
            ip1, j,   k,    &
            i,   jp1, k,    &
            ip1, jp1, k /), &
            shape(scvol_vertices), order = (/2, 1/))

       !Run over the 2 triangles
       do tk = 1, 2
          !Run over the labels of the vertices that
          !make up each triangle
          do tl = 1, 3 
             aux = triang_vertices_labels(tk, tl)
             ii = scvol_vertices(aux, 1)
             jj = scvol_vertices(aux, 2)
             kk = scvol_vertices(aux, 3)
             aux = mux_vector((/ii, jj, kk/), mesh, 1_k8)
             tmp = aux !Guaranteed to be > 0
             if(blocks) then
                !Which point in indexlist does aux correspond to?
                call binsearch(indexlist, aux, tmp) !tmp < 0 if search fails.
             end if
             triang(count, tl) = tmp

             if(tmp > 0) then
                !Save the mapping of a wave vector index to a (triangle, vertex)
                triangcount(tmp) = triangcount(tmp) + 1
                triangmap(1, tmp, triangcount(tmp)) = count
                triangmap(2, tmp, triangcount(tmp)) = tl
             end if
          end do
          count = count + 1
       end do
    end do
  end subroutine form_triangles

  subroutine fill_triangles(triang, evals, triang_evals)
    !! Populate the (sorted along the vertices) eigenvalues on all the vertices of the triangles
    !!
    !! triang List of the triangle vertices
    !! evals List of eigenvalues 
    !! triang_evals Triangles populated with the eigenvalues

    integer(k8), intent(in) :: triang(:,:)
    real(dp), intent(in) :: evals(:,:)
    real(dp), allocatable, intent(out) :: triang_evals(:,:,:)

    !Local variables
    integer(kind=4) :: iv, it, ib, numbands, aux, numtriangs, numvertices

    numtriangs = size(triang(:, 1))
    numbands = size(evals(1, :))
    numvertices = 3

    allocate(triang_evals(numtriangs, numbands, numvertices))

    !Note: Eigenvalues outside the transport active window is taken to be zero.
    !      As such, close to the transport window boundary, this method is inaccurate.
    !      A large enough transport window must be chosen to obtain accurate transport coefficients.
    triang_evals(:,:,:) = 0.0_dp
    
    do it = 1, numtriangs !Run over triangles
       do ib = 1, numbands !Run over bands
          do iv = 1, numvertices !Run over vertices
             aux = triang(it, iv)
             if(aux > 0) then !Only eigenvalues inside transport active region
                triang_evals(it, ib, iv) = evals(aux, ib)
             end if
          end do
          call sort(triang_evals(it, ib, :))
       end do
    end do
  end subroutine fill_triangles

  !pure real(dp) function delta_fn_triang(e, ik, ib, mesh, triangmap, triangcount, triang_evals)
  real(dp) function delta_fn_triang(e, ik, ib, mesh, triangmap, triangcount, triang_evals)
    !! Calculate delta function using the triangle method a la
    !! Kurganskii et al. Phys. Stat. Sol.(b) 129, 293 (1985)
    !!
    !! e Sample energy
    !! ik Wave vector index
    !! ib Band index
    !! mesh Wave vector grid
    !! triangmap Wave vector to (triangle, vertex) mapping
    !! triangcount Number of triangles in which a wave vector belongs
    !! triang_evals Triangles populated with the eigenvalues

    real(dp), intent(in) :: e
    integer(k8), intent(in) :: ik, ib
    integer(k8), intent(in) :: mesh(3), triangmap(:,:,:), triangcount(:)
    real(dp), intent(in) :: triang_evals(:,:,:)

    !Local variables
    integer(k8) :: iv, it, itk, num, numtriangs
    logical :: c1, c2, c3, c4
    real(dp) :: e1, e2, e3, E12, E21, E13, E31, E23, E32, tmp

    tmp = 0.0_dp
    delta_fn_triang = 0.0_dp

    !Total number of triangles in the system
    numtriangs = product(mesh)*2
    
    !Grab number of triangles in which wave vector belongs
    num = triangcount(ik)
    
    do itk = 1, num !Run over triangles
       it = triangmap(1, ik, itk) !Grab triangle
       iv = triangmap(2, ik, itk) !Grab vertex

       !Grab vertex energies
       e1 = triang_evals(it, ib, 1)
       e2 = triang_evals(it, ib, 2)
       e3 = triang_evals(it, ib, 3)

       !Evaluate the four possible cases
       c1 = e <= e1
       c2 = e1 < e .and. e <= e2
       c3 = e2 < e .and. e <= e3
       c4 = e3 < e

       tmp = 0.0_dp

       if(c1 .or. c4) cycle
       
       !Define Eij
       ! Note that at this stage the quantities below might
       ! be ill defined due to degeneracies. But the conditionals
       ! that will follow will take this into account.
       E12 = (e - e2)/(e1 - e2)
       E21 = (e - e1)/(e2 - e1)
       E13 = (e - e3)/(e1 - e3)
       E31 = (e - e1)/(e3 - e1)
       E23 = (e - e3)/(e2 - e3)
       E32 = (e - e2)/(e3 - e2)
              
       select case(iv)
       case(1)
          if(c2) then
             tmp = E21*(E12 + E13)/(e3 - e1)
          else if(c3) then
             tmp = E23*E13/(e3 - e1)
          end if
       case(2)
          if(c2) then
             tmp = E21*E21/(e3 - e1)
          else if(c3) then
             tmp = E23*E23/(e3 - e1)
          end if
       case(3)
          if(c2) then
             tmp = E21*E31/(e3 - e1)
          else if(c3) then
             tmp = E23*(E31 + E32)/(e3 - e1)
          end if
       end select

       delta_fn_triang = delta_fn_triang + tmp
    end do !itk
    
    if(delta_fn_triang < 1.0e-12_dp) delta_fn_triang = 0.0_dp
    
    !Normalize with the total number of triangles
    delta_fn_triang = delta_fn_triang/numtriangs
  end function delta_fn_triang
end module delta
