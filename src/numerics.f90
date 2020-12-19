module numerics_module
  !! Module containing type and procedures related to the numerics.

  use params, only: dp, k4, twopi
  use misc, only: exit_with_message

  implicit none

  private
  public numerics

  type numerics
     !! Data and procedures related to the numerics.

     integer(k4) :: qmesh(3)
     !! Phonon wave vector mesh.
     integer(k4) :: mesh_ref
     !! Electron mesh refinement factor compared to the phonon mesh.
     real(dp) :: fsthick
     !! Fermi surface thickness in (eV).
     character(len=500) :: cwd
     !character, allocatable :: cwd(:)
     !! Current working directory.
     character(len=500) ::datadumpdir
     !character, allocatable ::datadumpdir(:)
     !! Runtime data dump repository.
     character(len=500) :: g2dir
     !character, allocatable :: g2dir(:)
     !! Directory for e-ph vertex.
     character(len=500) :: Vdir
     !character, allocatable :: Vdir(:)
     !! Directory for ph-ph vertex.
     logical :: read_g2
     !! Choose if earlier e-ph vertex is to be used.
     logical :: read_V
     !! Choose if earlier p-ph vertex is to be used.
   
   contains

     procedure :: initialize=>read_input_and_setup
     
  end type numerics

contains

  subroutine read_input_and_setup(n)
    !! Read input file for information related to the numerics.

    class(numerics), intent(out) :: n

    !Local variables
    integer(k4) :: mesh_ref, qmesh(3)
    real(dp) :: fsthick
    character(len=500) :: datadumpdir
    logical :: read_g2, read_V

    namelist /numerics/ qmesh, mesh_ref, fsthick, datadumpdir, read_g2, read_V

    !Open input file
    open(1, file = 'input.nml', status = 'old')

    !Read numerics information
    n%qmesh = (/1, 1, 1/)
    n%mesh_ref = 1
    n%fsthick = 0.0_dp
    n%datadumpdir = './'
    n%read_g2 = .false.
    n%read_V = .false.
    read(1, nml = numerics)
    if(any(qmesh <= 0) .or. mesh_ref < 1 .or. fsthick < 0) then
       call exit_with_message('Bad input(s) in numerics.')
    end if
    n%qmesh = qmesh
    n%mesh_ref = mesh_ref
    n%fsthick = fsthick
    n%datadumpdir = trim(datadumpdir)
    
    ! Create data dump directory
    if(this_image() == 1) call system('mkdir ' // trim(adjustl(n%datadumpdir)))

    ! Create matrix elements data directories
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
