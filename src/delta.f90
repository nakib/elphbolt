module delta
  !! Module containing the procedures related to delta function evaulation.

  use params, only: dp, k4
  use misc, only: exit_with_message, print_message, demux_vector, mux_vector, &
       binsearch, sort
  
  implicit none

  private
  public form_tetrahedra_3d, fill_tetrahedra_3d!, delta_fn
  
contains
  
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

    integer(k4), intent(in) :: nk, mesh(3)
    integer(k4), intent(out), allocatable :: tetra(:,:), tetracount(:), tetramap(:,:,:)
    logical, intent(in) :: blocks
    integer(k4), optional, intent(in) :: indexlist(:)

    !Local variables
    integer(k4) :: ik, i, j, k, ijk(3), ii, jj, kk, tk, tl, aux, count
    integer(k4) :: ip1, jp1, kp1, n1, n2, n3, tmp
    integer(k4), dimension(6,4) :: tetra_vertices_labels
    integer(k4), dimension(8,3) :: scvol_vertices ! subcell volume vertices

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
    
    do ik = 1, nk !run over all wave vectors in FBZ
       if(blocks) then !for energy window restricted FBZ
          call demux_vector(indexlist(ik) - 1, ijk, mesh, 1_k4)
       else !for unrestristed FBZ
          call demux_vector(ik, ijk, mesh, 1_k4)
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

       !For each subcell get the vertices
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

       do tk = 1, 6 ! run over 6 tetrahedra
          do tl = 1, 4 ! run over the labels of the vertices that
             !Make up each tetrahedron
             aux = tetra_vertices_labels(tk, tl)
             ii = scvol_vertices(aux,1)
             jj = scvol_vertices(aux,2)
             kk = scvol_vertices(aux,3)
             aux = mux_vector((/ii, jj, kk/), mesh, 1_k4)
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

    integer(k4), intent(in) :: tetra(:,:)
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

!!$  !populate coefficients tetra_c for given energy e and a single band combo, for a single q index
!!$  subroutine fillTetraWeight_e_iq(e,i,iq,weight,sp)
!!$    real(kind=8), intent(in) :: e
!!$    integer(kind=4), intent(in) :: i
!!$    integer(kind=4), intent(in) :: iq
!!$    integer(kind=4), intent(in) :: sp !0 for phonons, 1 for electrons
!!$
!!$    real(kind=8), intent(out) :: weight
!!$
!!$    integer(kind=4) :: iv,itet,itetq,num
!!$    logical :: c1,c2,c3
!!$    real(kind=8) :: e21,e31,e41,e32,e42,e43,tmp ! eji \equiv ej - ei
!!$
!!$    tmp = 0.0
!!$    weight = 0.0
!!$
!!$    if(sp .eq. 0) then
!!$       num = q2tetcount(iq)
!!$    else
!!$       num = q2tetcount_blocks(iq)
!!$    end if
!!$
!!$    do itetq = 1,num
!!$       if(sp .eq. 0) then
!!$          itet = q2tet(1,iq,itetq)
!!$          iv = q2tet(2,iq,itetq)
!!$       else
!!$          itet = q2tet_blocks(1,iq,itetq)
!!$          iv = q2tet_blocks(2,iq,itetq)
!!$       end if
!!$
!!$       if(sp .eq. 0) then !phonons
!!$          e1 = tetra_evals_ph(itet,i,1)
!!$          e2 = tetra_evals_ph(itet,i,2)
!!$          e3 = tetra_evals_ph(itet,i,3)
!!$          e4 = tetra_evals_ph(itet,i,4)
!!$       else
!!$          e1 = tetra_evals_el(itet,i,1)
!!$          e2 = tetra_evals_el(itet,i,2)
!!$          e3 = tetra_evals_el(itet,i,3)
!!$          e4 = tetra_evals_el(itet,i,4)
!!$       end if
!!$
!!$       e1e = e1-e
!!$       e2e = e2-e
!!$       e3e = e3-e
!!$       e4e = e4-e
!!$       e21 = e2-e1
!!$       e31 = e3-e1
!!$       e41 = e4-e1
!!$       e32 = e3-e2
!!$       e42 = e4-e2
!!$       e43 = e4-e3
!!$
!!$       !eval the three cases
!!$       c1 = e1 <= e .and. e <= e2
!!$       c2 = e2 <= e .and. e <= e3
!!$       c3 = e3 <= e .and. e <= e4
!!$
!!$       if(.not. (e < e1 .or. e > e4)) then
!!$
!!$          select case(iv)
!!$          case(1)
!!$             if(c1) then
!!$                tmp = (e2e/e21 + e3e/e31 + e4e/e41)*(e1e**2)/e41/e31/e21
!!$
!!$                if(e1 .eq. e2) then
!!$                   tmp = 0.0
!!$                end if
!!$
!!$             else if(c2) then
!!$                tmp = -0.5*(e3e/(e31**2)*(e3e*e2e/e42/e32 + e4e*e1e/e41/e42 + e3e*e1e/e32/e41) &
!!$                     + e4e/(e41**2)*(e4e*e1e/e42/e31 + e4e*e2e/e42/e32 + e3e*e1e/e31/e32))
!!$
!!$                if(e2 .eq. e3) then
!!$                   tmp = -0.5*(e4e*e1e/e41/e42 + e1e/e41 &
!!$                        + e4e/(e41**2)*(e4e*e1e/e42/e31 + e4e/e42 + e1e/e31))
!!$                end if
!!$
!!$             else if(c3) then
!!$                tmp = (e4e**3)/(e41**2)/e42/e43
!!$
!!$                if(e3 .eq. e4) then
!!$                   tmp = (e4e**2)/(e41**2)/e42
!!$                end if
!!$
!!$             end if
!!$          case(2)
!!$             if(c1) then
!!$                tmp = -(e1e**3)/(e21**2)/e31/e41
!!$
!!$                if(e1 .eq. e2) then
!!$                   tmp = 0.0
!!$                end if
!!$
!!$             else if(c2) then
!!$                tmp = -0.5*(e3e/(e32**2)*(e3e*e2e/e42/e31 + e4e*e2e/e42/e41 + e3e*e1e/e31/e41) &
!!$                     + e4e/(e42**2)*(e3e*e2e/e32/e31 + e4e*e1e/e41/e31 + e4e*e2e/e32/e41))
!!$
!!$                if(e2 .eq. e3) then
!!$                   tmp = -0.5*(0.0 + e4e/e42/e41 + 0.0 &
!!$                        + e4e/(e42**2)*(0.0 + e4e*e1e/e41/e31 + 1.0))
!!$                end if
!!$
!!$             else if(c3) then
!!$                tmp = (e4e**3)/e41/(e42**2)/e43
!!$
!!$                if(e3 .eq. e4) then
!!$                   tmp = 0.0
!!$                end if
!!$             end if
!!$          case(3)
!!$             if(c1) then
!!$                tmp = -(e1e**3)/e21/(e31**2)/e41
!!$
!!$                if(e1 .eq. e2) then
!!$                   tmp = 0.0
!!$                end if
!!$
!!$             else if(c2) then
!!$                tmp = 0.5*(e2e/(e32**2)*(e3e*e2e/e42/e31 + e4e*e2e/e42/e41 + e3e*e1e/e31/e41) &
!!$                     + e1e/(e31**2)*(e3e*e2e/e42/e32 + e4e*e1e/e41/e42 + e3e*e1e/e32/e41))
!!$
!!$                if(e2 .eq. e3) then
!!$                   tmp = 0.5*(0.0 + e4e/e42/e41 + e1e/e31/e41 &
!!$                        + e1e/(e31**2)*(0.0 + e4e*e1e/e41/e42 + e1e/e41))
!!$                end if
!!$
!!$             else if(c3) then
!!$                tmp = (e4e**3)/e41/e42/(e43**2)
!!$
!!$                if(e3 .eq. e4) then
!!$                   tmp = 0.0
!!$                end if
!!$             end if
!!$          case(4)
!!$             if(c1) then
!!$                tmp = -(e1e**3)/e21/e31/(e41**2)
!!$                if(e1 .eq. e2) then
!!$                   tmp = 0.0
!!$                end if
!!$             else if(c2) then
!!$                tmp = 0.5*(e2e/(e42**2)*(e3e*e2e/e32/e31 + e4e*e1e/e41/e31 + e4e*e2e/e32/e41) &
!!$                     + e1e/(e41**2)*(e4e*e1e/e42/e31 + e4e*e2e/e42/e32 + e3e*e1e/e31/e32))
!!$
!!$                if(e2 .eq. e3) then
!!$                   tmp = 0.5*(0.0 &
!!$                        + e1e/(e41**2)*(e4e*e1e/e42/e31 + e4e/e42 + e1e/e31))
!!$                end if
!!$
!!$             else if(c3) then
!!$                tmp = -(e3e/e43 + e2e/e42 + e1e/e41)*(e4e**2)/e41/e42/e43
!!$
!!$                if(e3 .eq. e4) then
!!$                   tmp = 0.0
!!$                end if
!!$             end if
!!$          end select
!!$
!!$          if ((e1 .eq. e2) .and. (e1 .eq. e3) .and. (e1 .eq. e4) .and. (e .eq. e1)) then
!!$             tmp = 0.25
!!$          end if
!!$
!!$          weight = weight + tmp
!!$
!!$       end if ! .not. (e <= e1 .or. e >= e4)
!!$    end do !itetq
!!$
!!$    if(weight < 1.d-12) weight = 0.d0
!!$
!!$    if(sp .eq. 0) then !phonons
!!$       weight = weight/ntetra
!!$    else
!!$       weight = weight/product(kfgrid)/6.d0
!!$    end if
!!$
!!$  end subroutine fillTetraWeight_e_iq
end module delta
