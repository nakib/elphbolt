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
     logical :: phbte
     !! Choose if phonon BTE will be solved.
     logical :: ebte
     !! Choose if electron BTE will be solved.
     logical :: drag
     !! Choose if the drag effect will be included.
     integer(k8) :: maxiter
     !! Maximum number of iterations in the BTE solver.
     real(dp) :: conv_thres
     !! BTE iteration convergence criterion.
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
    logical :: read_gq2, read_gk2, read_V, tetrahedra, phe, phbte, ebte, drag

    namelist /numerics/ qmesh, mesh_ref, fsthick, datadumpdir, read_gq2, read_gk2, &
         read_V, tetrahedra, phe, phbte, ebte, maxiter, conv_thres, drag

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
    tetrahedra = .false.
    phe = .true.
    phbte = .true.
    ebte = .true.
    drag = .true.
    maxiter = 50
    conv_thres = 0.0001_dp
    read(1, nml = numerics)
    if(any(qmesh <= 0) .or. mesh_ref < 1 .or. fsthick < 0 .or. .not.(phbte .or. ebte)) then
       call exit_with_message('Bad input(s) in numerics.')
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
    n%phbte = phbte
    n%ebte = ebte
    n%maxiter = maxiter
    n%conv_thres = conv_thres
    n%drag = drag
    
    !Create data dump directory
    if(this_image() == 1) call system('mkdir ' // trim(adjustl(n%datadumpdir)))

    !Create matrix elements data directories
    n%g2dir = trim(adjustl(n%datadumpdir))//'g2'
    if(this_image() == 1) call system('mkdir ' // trim(adjustl(n%g2dir)))
    n%Vdir = trim(adjustl(n%datadumpdir))//'V'
    if(this_image() == 1) call system('mkdir ' // trim(adjustl(n%Vdir)))
    
    !Close input file
    close(1)

    !Set current work directory.
    call getcwd(n%cwd)
    n%cwd = trim(n%cwd)
  end subroutine read_input_and_setup
end module numerics_module
