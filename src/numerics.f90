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

module numerics_module
  !! Module containing type and procedures related to the numerics.

  use params, only: dp, k8, twopi
  use misc, only: exit_with_message, subtitle

  implicit none

  private
  public numerics

  type numerics
     !! Data and procedures related to the numerics.

     integer(k8) :: qmesh(3)
     !! Phonon wave vector mesh.
     integer(k8) :: mesh_ref
     !! Electron mesh refinement factor compared to the phonon mesh.
     real(dp) :: fsthick
     !! Fermi surface thickness in (eV).
     character(len = 1024) :: cwd
     !character, allocatable :: cwd(:)
     !! Current working directory.
     character(len = 1024) ::datadumpdir
     !character, allocatable ::datadumpdir(:)
     !! Runtime data dump repository.
     character(len = 1024) :: g2dir
     !character, allocatable :: g2dir(:)
     !! Directory for e-ph vertex.
     character(len = 1024) :: Vdir
     !character, allocatable :: Vdir(:)
     !! Directory for ph-ph vertex.
     logical :: read_gq2
     !! Choose if earlier e-ph (IBZ q) vertices are to be used.
     logical :: read_gk2
     !! Choose if earlier e-ph (IBZ k) vertices are to be used.
     logical :: read_V
     !! Choose if earlier ph-ph vertices are to be used.
     logical :: tetrahedra
     !! Choose if the tetrahedron method for delta function evaluation will be used.
     logical :: phe
     !! Choose if ph-e interaction will be included.
     logical :: phiso
     !! Use phonon-isotope scattering?
     logical :: onlyphbte
     !! Choose if only phonon BTE will be solved.
     logical :: onlyebte
     !! Choose if electron BTE will be solved.
     logical :: elchimp
     !! Use electron-charged impurity scattering?
     logical :: drag
     !! Choose if the drag effect will be included.
     integer(k8) :: maxiter
     !! Maximum number of iterations in the BTE solver.
     real(dp) :: conv_thres
     !! BTE iteration convergence criterion.
     logical :: plot_along_path
     !! Plot Wannierized quantities along high symmetry wave vectors?
   contains

     procedure :: initialize=>read_input_and_setup
     
  end type numerics

contains

  subroutine read_input_and_setup(n)
    !! Read input file for information related to the numerics.

    class(numerics), intent(out) :: n

    !Local variables
    integer(k8) :: mesh_ref, qmesh(3), maxiter
    real(dp) :: fsthick, conv_thres
    character(len = 1024) :: datadumpdir
    logical :: read_gq2, read_gk2, read_V, tetrahedra, phe, phiso, onlyphbte, &
         onlyebte, elchimp, drag, plot_along_path

    namelist /numerics/ qmesh, mesh_ref, fsthick, datadumpdir, read_gq2, read_gk2, &
         read_V, tetrahedra, phe, phiso, onlyphbte, onlyebte, maxiter, conv_thres, drag, &
         elchimp, plot_along_path

    call subtitle("Reading numerics information...")
    
    !Open input file
    open(1, file = 'input.nml', status = 'old')

    !Read numerics information
    qmesh = (/1, 1, 1/)
    mesh_ref = 1
    fsthick = 0.0_dp
    datadumpdir = './'
    read_gq2 = .false.
    read_gk2 = .false.
    read_V = .false.
    tetrahedra = .true.
    phe = .false.
    phiso = .false.
    onlyphbte = .false.
    onlyebte = .false.
    elchimp = .false.
    drag = .true.
    plot_along_path = .false.
    maxiter = 50
    conv_thres = 1e-4_dp
    read(1, nml = numerics)
    if(any(qmesh <= 0) .or. mesh_ref < 1 .or. fsthick < 0) then
       call exit_with_message('Bad input(s) in numerics.')
    end if
    !if(.not. onlyphbte .or. .not. onlyebte .or. drag) then
    !   call exit_with_message('Bad input(s) in numerics.')
    !end if
    if(.not. tetrahedra) then
       call exit_with_message('At present only the tetrahedra method is supported.')
    end if
    n%qmesh = qmesh
    n%mesh_ref = mesh_ref
    n%fsthick = fsthick
    n%datadumpdir = trim(datadumpdir)
    n%read_gq2 = read_gq2
    n%read_gk2 = read_gk2
    n%read_V = read_V
    n%tetrahedra = tetrahedra
    n%phe = phe
    n%phiso = phiso
    n%onlyphbte = onlyphbte
    n%onlyebte = onlyebte
    n%elchimp = elchimp
    n%maxiter = maxiter
    n%conv_thres = conv_thres
    n%drag = drag
    n%plot_along_path = plot_along_path

    !Set BTE solution type
    if(n%onlyphbte) then
       n%onlyebte = .false.
       n%drag = .false.
    end if
    if(n%onlyebte) then
       n%onlyphbte = .false.
       n%drag = .false.
       n%phe = .false.
    end if
    if(n%drag) then
       n%onlyebte = .false.
       n%onlyphbte = .false.
       n%phe = .true.
    end if
    
    !Create data dump directory
    if(this_image() == 1) call system('mkdir -p ' // trim(adjustl(n%datadumpdir)))

    !Create matrix elements data directories
    n%g2dir = trim(adjustl(n%datadumpdir))//'g2'
    if(this_image() == 1) call system('mkdir -p ' // trim(adjustl(n%g2dir)))
    n%Vdir = trim(adjustl(n%datadumpdir))//'V'
    if(this_image() == 1) call system('mkdir -p ' // trim(adjustl(n%Vdir)))
    
    !Close input file
    close(1)

    !Set current work directory.
    call getcwd(n%cwd)
    n%cwd = trim(n%cwd)

    !Print out information.
    if(this_image() == 1) then
       write(*, "(A, (3I5,x))") "q-mesh = ", n%qmesh
       write(*, "(A, (3I5,x))") "k-mesh = ", n%mesh_ref*n%qmesh
       write(*, "(A, 1E16.8, A)") "Fermi window thickness = ", n%fsthick, " eV"
       write(*, "(A, A)") "Working directory = ", trim(n%cwd)
       write(*, "(A, A)") "Data dump directory = ", trim(n%datadumpdir)
       write(*, "(A, A)") "e-ph directory = ", trim(n%g2dir)
       write(*, "(A, A)") "ph-ph directory = ", trim(n%Vdir)
       write(*, "(A, L)") "Reuse e-ph matrix elements: ", n%read_gk2
       write(*, "(A, L)") "Reuse ph-e matrix elements: ", n%read_gq2
       write(*, "(A, L)") "Reuse ph-ph matrix elements: ", n%read_V
       write(*, "(A, L)") "Use tetrahedron method: ", n%tetrahedra
       write(*, "(A, L)") "Include ph-e interaction: ", n%phe
       write(*, "(A, L)") "Include ph-isotope interaction: ", n%phiso
       write(*, "(A, L)") "Include electron-charged impurity interaction: ", n%elchimp
       if(n%onlyphbte) write(*, "(A, L)") "Calculate only phonon BTE: ", n%onlyphbte
       if(n%onlyebte) write(*, "(A, L)") "Calculate only electron BTE: ", n%onlyebte
       write(*, "(A, L)") "Include drag: ", n%drag
       write(*, "(A, L)") "Plot quantities along path: ", n%plot_along_path
       write(*, "(A, I5)") "Maximum number of BTE iterations = ", n%maxiter
       write(*, "(A, 1E16.8)") "BTE convergence threshold = ", n%conv_thres
    end if
  end subroutine read_input_and_setup
end module numerics_module
