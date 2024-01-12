! Copyright 2020 elphbolt contributors.
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

  use params, only: r64, i64, twopi
  use misc, only: exit_with_message, subtitle
  use crystal_module, only: crystal

  implicit none
  
  private
  public numerics

  !external system, getcwd
  
  type numerics
     !! Data and procedures related to the numerics.

     integer(i64) :: qmesh(3)
     !! Phonon wave vector mesh.
     integer(i64) :: mesh_ref
     !! Electron mesh refinement factor compared to the phonon mesh.
     logical :: eco_mode
     !! Use coarse grained matrix elements?
     integer(i64) :: econess(3)
     !! How much coarse graining for the matrix elements?
     real(r64) :: fsthick
     !! Fermi surface thickness in eV.
     character(len = 1024) :: cwd
     !! Current working directory.
     character(len = 1024) ::datadumpdir
     !! Runtime data dump repository.
     character(len = 1024) ::datadumpdir_T
     !! Runtime temperature dependent data dump repository.
     character(len = 1024) ::datadumpdir_T_chempot
     !! Runtime temperature and chemical potential dependent data dump repository.
     character(len = 1024) :: g2dir
     !! Directory for e-ph vertex.
     character(len = 1024) :: Vdir
     !! Directory for ph-ph vertex.
     character(len = 1024) :: Wdir
     !! Directory for ph-ph transition rates.
     character(len = 1024) :: Xdir
     !! Directory for e-ph transition rates.
     character(len = 1024) :: Ydir
     !! Directory for ph-e transition rates.
     character(len = 1024) :: scdir
     !! Directory for the superconductivity temporary data.
     character(len = 1024) :: epsilondir
     !! Directory for the dielectric function data.
     logical :: read_gq2
     !! Choose if earlier e-ph (IBZ q) vertices are to be used.
     logical :: read_gk2
     !! Choose if earlier e-ph (IBZ k) vertices are to be used.
     logical :: read_V
     !! Choose if earlier ph-ph (IBZ q) vertices are to be used.
     logical :: read_W
     !! Choose if earlier ph-ph (IBZ q) transition probabilities are to be used.
     logical :: tetrahedra
     !! Choose if the tetrahedron method for 3d delta function evaluation will be used.
     logical :: phe
     !! Choose if ph-e interaction will be included.
     logical :: phiso
     !! Use phonon-isotope scattering?
     character(len = 6) :: phiso_1B_theory
     !! Choose phonon-isotope 1st Born scattering theory from "Tamura" or "DIB-1B".
     logical :: phiso_Tmat
     !! Calculate phonon-isotope scattering from general subsitutional defect T-matrix formalism?
     logical :: phsubs
     !! Use phonon-substitution scattering?
     logical :: phbound
     !! Use phonon-boundary scattering?
     logical :: fourph
     !! Use 4-ph scattering?
     integer(i64) :: fourph_mesh_ref
     !! Mesh refinement factor of phonon wavectors with respect to external 4-ph calculation 
     logical :: phthinfilm
     !! Use phonon-thin-film scattering?
     logical :: phthinfilm_ballistic
     !! Use ballistic limit of Fuchs-Sondheimer formula?
     logical :: phdef_Tmat
     !! Calculate phonon-defect scattering T-matrix?
     logical :: onlyphbte
     !! Choose if only phonon BTE will be solved.
     logical :: onlyebte
     !! Choose if electron BTE will be solved.
     logical :: elchimp
     !! Use electron-charged impurity scattering?
     logical :: elbound
     !! Use electron-boundary scattering?
     logical :: drag
     !! Choose if the drag effect will be included.
     logical :: use_Wannier_ifc2s
     !! Choose if Wannier 2nd order ifcs will be used for usual phonon calculations.
     integer(i64) :: maxiter
     !! Maximum number of iterations in the BTE/Migdal-Eliashberg equations solver.
     real(r64) :: conv_thres
     !! BTE/Migdal-Eliashberg euqations iteration convergence criterion.
     logical :: plot_along_path
     !! Plot Wannierized quantities along high symmetry wave vectors?
     integer(i64) :: runlevel
     !! Control for the type of calculation.
     real(r64) :: ph_en_min, ph_en_max
     !! Bounds of equidistant phonon energy mesh.
     integer(i64) :: ph_en_num
     !! Number of equidistant phonon energy mesh points.
     real(r64) :: el_en_min, el_en_max
     !! Bounds of equidistant electron energy mesh.
     integer(i64) :: el_en_num
     !! Number of equidistant electron energy mesh points.
     integer(i64) :: ph_mfp_npts
     !! Number of equidistant phonon mean-free-path mesh points.
     logical :: Bfield_on
     !! Is B-field on?
     real(r64) :: Bfield(3)
     !! B-field vector
   contains

     procedure :: initialize=>read_input_and_setup, create_chempot_dirs
     
  end type numerics

contains

  subroutine read_input_and_setup(self, crys)
    !! Read input file for information related to the numerics.
    !!
    !! self Numerics object
    !! crys Crytal object

    class(numerics), intent(out) :: self
    type(crystal), intent(in) :: crys
    
    !Local variables
    integer(i64) :: mesh_ref, qmesh(3), maxiter, runlevel, el_en_num, &
         ph_en_num, ph_mfp_npts, fourph_mesh_ref
    integer :: i
    real(r64) :: fsthick, conv_thres, ph_en_min, ph_en_max, el_en_min, el_en_max, Bfield(3)
    character(len = 1024) :: datadumpdir, tag
    character(len = 6) :: phiso_1B_theory
    character(len = 1) :: numcols
    logical :: read_gq2, read_gk2, read_V, read_W, tetrahedra, phe, phiso, phsubs, &
         phbound, phdef_Tmat, onlyphbte, onlyebte, elchimp, elbound, drag, plot_along_path, &
         phthinfilm, phthinfilm_ballistic, fourph, use_Wannier_ifc2s, phiso_Tmat, Bfield_on

    namelist /numerics/ qmesh, mesh_ref, fsthick, datadumpdir, read_gq2, read_gk2, &
         read_V, read_W, tetrahedra, phe, phiso, phsubs, onlyphbte, onlyebte, maxiter, &
         conv_thres, drag, elchimp, plot_along_path, runlevel, ph_en_min, ph_en_max, &
         ph_en_num, el_en_min, el_en_max, el_en_num, phbound, elbound, phdef_Tmat, &
         ph_mfp_npts, phthinfilm, phthinfilm_ballistic, fourph, fourph_mesh_ref, use_Wannier_ifc2s, &
         phiso_Tmat, phiso_1B_theory, Bfield_on, Bfield

    call subtitle("Reading numerics information...")
    
    !Open input file
    open(1, file = 'input.nml', status = 'old')

    !Read numerics information
    qmesh = [1, 1, 1]
    mesh_ref = 1
    fourph_mesh_ref = 1
    fsthick = 0.0_r64
    datadumpdir = './'
    read_gq2 = .false.
    read_gk2 = .false.
    read_V = .false.
    read_W = .false.
    tetrahedra = .false.
    phe = .false.
    phiso = .false.
    phiso_1B_theory = "DIB-1B"
    phiso_Tmat = .false.
    phsubs = .false.
    phbound = .false.
    fourph = .false.
    phthinfilm = .false.
    phthinfilm_ballistic = .true.
    phdef_Tmat = .false.
    onlyphbte = .false.
    onlyebte = .false.
    elchimp = .false.
    elbound = .false.
    drag = .true.
    use_Wannier_ifc2s = .false.
    plot_along_path = .false.
    Bfield_on = .false.
    Bfield = [0.0_r64, 0.0_r64, 0.0_r64]
    maxiter = 50
    conv_thres = 1e-4_r64
    runlevel = 1
    ph_en_min = 0.0_r64
    ph_en_max = 1.0_r64
    ph_en_num = 100
    el_en_min = -10.0_r64
    el_en_max = 10.0_r64
    el_en_num = 100
    ph_mfp_npts = 100
    read(1, nml = numerics)

    if(any(qmesh <= 0) .or. fourph_mesh_ref < 1 .or. mesh_ref < 1 .or. fsthick < 0) then
       call exit_with_message('Bad input(s) in numerics.')
    end if

    if(crys%twod .and. tetrahedra) then
       call exit_with_message('The tetrahedra method only works for 3d. Exiting.')
    end if

    if(crys%epsilon0 == 0.0 .and. elchimp) then
       call exit_with_message(&
            'Need to provide non-zero epsilon0 for e-ch. imp. interaction. Exiting.')
    end if

    if(phbound .and. phthinfilm) then
       call exit_with_message('ph-boundary and ph-thin-film scattering not allowed together. Exiting.')
    end if

    if(phiso) then
       if((phiso_1B_theory /= "Tamura") .and. (phiso_1B_theory /= "DIB-1B")) then
          call exit_with_message("phiso_1B_theory can be either 'Tamura' or 'DIB-1B'. Exiting.")
       end if
       if(phiso_1B_theory == "DIB-1B" .and. crys%VCA) then
          call exit_with_message("phiso_1B_theory can't be 'DIB-1B' if 'VCA' is true. Exiting.")
       end if
       if(phiso_1B_theory == "Tamura" .and. crys%DIB) then
          call exit_with_message("phiso_1B_theory can't be 'Tamura' if 'DIB' is true. Exiting.")
       end if
    end if

    !Only allow B-field along z for now.
    if(Bfield_on) then
       if(.not. (Bfield(1) == 0.0_r64 .and. Bfield(2) == 0.0_r64 .and. Bfield(3) /= 0.0_r64)) then
          call exit_with_message("B-field has to be of the form [0 0 B]. Exiting.")
       end if
!!$       if( .not. (Bfield(1) /= 0.0_r64 .and. Bfield(2) == 0.0_r64 .and. Bfield(3) == 0.0_r64) .or. &
!!$            (Bfield(1) == 0.0_r64 .and. Bfield(2) /= 0.0_r64 .and. Bfield(3 == 0.0_r64)) .or. &
!!$            (Bfield(1) == 0.0_r64 .and. Bfield(2) == 0.0_r64 .and. Bfield(3 /= 0.0_r64))) then
!!$          call exit_with_message("B-field has to be of the form [B 0 0], [0 B 0], or [0 0 B]. Exiting.")
!!$       end if
    end if

    !TODO
    !! [ ] Read eco mode info from input
    !! [ ] Check for valid choice of econess
    
    !eco mode DBG
    self%eco_mode = .true.
    self%econess = [2, 2, 2]
    !!

    self%Bfield_on = Bfield_on
    self%Bfield = Bfield
    
    self%qmesh = qmesh
    self%runlevel = runlevel
    !Runlevels:
    !1 BTE
    !2 BTE postproc
    !3 Superconductivity
    if(self%runlevel /= 3) then !Non-superconductivity mode
       self%mesh_ref = mesh_ref
       self%fourph_mesh_ref = fourph_mesh_ref
       self%read_gq2 = read_gq2
       self%read_V = read_V
       self%read_W = read_W
       self%phe = phe
       self%phiso = phiso
       self%phiso_1B_theory = phiso_1B_theory
       self%phiso_Tmat = phiso_Tmat
       self%phsubs = phsubs
       self%phbound = phbound
       self%fourph = fourph
       self%phthinfilm = phthinfilm
       self%phthinfilm_ballistic = phthinfilm_ballistic
       self%phdef_Tmat = phdef_Tmat
       self%onlyphbte = onlyphbte
       self%onlyebte = onlyebte
       self%elchimp = elchimp
       self%elbound = elbound
       self%drag = drag
    else
       self%mesh_ref = 1 !Enforce this for superconductivity mode
    end if
    self%use_Wannier_ifc2s = use_Wannier_ifc2s
    self%read_gk2 = read_gk2
    self%fsthick = fsthick
    self%datadumpdir = trim(datadumpdir)
    self%tetrahedra = tetrahedra
    self%maxiter = maxiter
    self%conv_thres = conv_thres
    self%plot_along_path = plot_along_path

    if(runlevel == 2) then
       self%ph_en_min = ph_en_min
       self%ph_en_max = ph_en_max
       self%ph_en_num = ph_en_num
       self%el_en_min = el_en_min
       self%el_en_max = el_en_max
       self%el_en_num = el_en_num
       self%ph_mfp_npts = ph_mfp_npts
    end if
    
    if(crys%twod .and. self%qmesh(3) /= 1) then
       call exit_with_message('For 2d systems, qmesh(3) must be equal to 1.')
    end if
    
    !Set BTE solution type
    if(self%onlyphbte) then
       self%onlyebte = .false.
       self%drag = .false.
    end if
    if(self%onlyebte) then
       self%onlyphbte = .false.
       self%phiso = .false.
       self%phiso_Tmat = .false.
       self%phsubs = .false.
       self%drag = .false.
       self%phe = .false.
    end if
    if(self%drag) then
       self%onlyebte = .false.
       self%onlyphbte = .false.
       self%phe = .true.
    end if

    !Check if T-matrix and tetrahedron method consistency
    if(self%phdef_Tmat .and. .not. self%tetrahedra) then
       call exit_with_message("Currently T-matrix method is only supported with tetrahedron method. Exiting.")
    end if

    if(self%phiso_Tmat .and. .not. self%phdef_Tmat) then
       call exit_with_message("For ph-iso scattering from T-matrix, need both phiso_Tmat and phdef_Tmat. Exiting.")
    end if
    
    !Create data dump directory
    if(this_image() == 1) call system('mkdir -p ' // trim(adjustl(self%datadumpdir)))

    !Create matrix elements data directories
    self%g2dir = trim(adjustl(self%datadumpdir))//'g2'
    if(this_image() == 1) call system('mkdir -p ' // trim(adjustl(self%g2dir)))
    self%Vdir = trim(adjustl(self%datadumpdir))//'V2'
    if(this_image() == 1) call system('mkdir -p ' // trim(adjustl(self%Vdir)))

    !Create superconductivity data dump directory
    self%scdir = trim(adjustl(self%datadumpdir))//'sc'
    if(this_image() == 1) call system('mkdir -p ' // trim(adjustl(self%scdir)))

    !Create T dependent subdirectory
    write(tag, "(E9.3)") crys%T
    self%datadumpdir_T = trim(adjustl(self%datadumpdir))//'T'//trim(adjustl(tag))
    if(this_image() == 1) call system('mkdir -p ' // trim(adjustl(self%datadumpdir_T)))

    !Create T-dependent ph-ph transition probability directory
    self%Wdir = trim(adjustl(self%datadumpdir_T))//'/W'
    if(this_image() == 1) call system('mkdir -p ' // trim(adjustl(self%Wdir)))
    
    !Close input file
    close(1)

    !Set current work directory.
    call getcwd(self%cwd)
    self%cwd = trim(self%cwd)

    !Print out information.
    if(this_image() == 1) then
       write(numcols, "(I0)") 3
       open(1, file = "ph.reclattvecs", status = "replace")
       do i = 1, 3
          write(1, "(" // trim(adjustl(numcols)) // "E20.10)") &
               crys%reclattvecs(:, i)
       end do
       write(1, *) self%qmesh
       close(1)

       open(1, file = "el.reclattvecs", status = "replace")
       do i = 1, 3
          write(1, "(" // trim(adjustl(numcols)) // "E20.10)") &
               crys%reclattvecs(:, i)
       end do
       if(crys%twod) then
          write(1, *) self%mesh_ref*self%qmesh(1), self%mesh_ref*self%qmesh(2), 1
       else
          write(1, *) self%mesh_ref*self%qmesh(1), self%mesh_ref*self%qmesh(2), &
               self%mesh_ref*self%qmesh(3)
       end if
       close(1)
       
       write(*, "(A, (3I5,x))") "q-mesh = ", self%qmesh
       if(crys%twod) then
          write(*, "(A, (3I5,x))") "k-mesh = ", self%mesh_ref*self%qmesh(1), self%mesh_ref*self%qmesh(2), 1
       else
          write(*, "(A, (3I5,x))") "k-mesh = ", self%mesh_ref*self%qmesh(1), self%mesh_ref*self%qmesh(2), &
               self%mesh_ref*self%qmesh(3)
       end if
       if(self%fourph) then
          if(crys%twod) then
             write(*, "(A, (3I5,x))") "external 4ph q-mesh = ", &
                  self%qmesh(1)/self%fourph_mesh_ref, self%qmesh(2)/self%fourph_mesh_ref, 1
          else
             write(*, "(A, (3I5,x))") "external 4ph q-mesh = ", &
                  self%qmesh(1)/self%fourph_mesh_ref, self%qmesh(2)/self%fourph_mesh_ref, &
                  self%qmesh(3)/self%fourph_mesh_ref
          end if
       end if
       write(*, "(A, 1E16.8, A)") "Fermi window thickness (each side of reference energy) = ", self%fsthick, " eV"
       write(*, "(A, A)") "Working directory = ", trim(self%cwd)
       write(*, "(A, A)") "Data dump directory = ", trim(self%datadumpdir)
       write(*, "(A, A)") "T-dependent data dump directory = ", trim(self%datadumpdir_T)
       write(*, "(A, A)") "e-ph directory = ", trim(self%g2dir)
       write(*, "(A, L)") "Use Wannier ifcs: ", self%use_Wannier_ifc2s
       if(self%runlevel /= 3) write(*, "(A, A)") "ph-ph directory = ", trim(self%Vdir)
       write(*, "(A, L)") "Reuse e-ph matrix elements: ", self%read_gk2
       if(self%runlevel /= 3) then
          write(*, "(A, L)") "Reuse ph-e matrix elements: ", self%read_gq2
          write(*, "(A, L)") "Reuse ph-ph matrix elements: ", self%read_V
          write(*, "(A, L)") "Reuse ph-ph transition probabilities: ", self%read_W
          write(*, "(A, L)") "Calculate ph-e interaction: ", self%phe
          write(*, "(A, L)") "Calculate 1st Born ph-isotope interaction: ", self%phiso
          if(self%phiso) &
               write(*, "(A, A)") "Selected ph-isotope 1st Born theory: ", self%phiso_1B_theory
          write(*, "(A, L)") "Calculate ph-isotope interaction via T-matrix: ", self%phiso_Tmat
          write(*, "(A, L)") "Include ph-substitution interaction: ", self%phsubs
          write(*, "(A, L)") "Include ph-boundary interaction: ", self%phbound
          write(*, "(A, L)") "Include 4-ph interaction: ", self%fourph
          write(*, "(A, L)") "Include ph-thin-film interaction: ", self%phthinfilm
          if(self%phthinfilm) write(*, "(A, L)") "    Use ballistic limit: ", self%phthinfilm_ballistic
          write(*, "(A, L)") "Include ph-defect interaction using the T-matrix: ", self%phdef_Tmat
          if(self%phbound) then
             write(*,"(A,(1E16.8,x),A)") 'Characteristic length for ph-boundary scattering =', &
                  crys%bound_length, 'mm'
          end if
          if(self%phthinfilm) then
             write(*,"(A,(1E16.8,x),A,A,A)") 'Height for ph-thin-film scattering =', &
                  crys%thinfilm_height, 'mm along the ', crys%thinfilm_normal, ' direction'
          end if
          write(*, "(A, L)") "Include el-charged impurity interaction: ", self%elchimp
          write(*, "(A, L)") "Include el-boundary interaction: ", self%elbound
          if(self%elbound) then
             write(*,"(A,(1E16.8,x),A)") 'Characteristic length for el-boundary scattering =', &
                  crys%bound_length, 'mm'
          end if
          if(self%onlyphbte) write(*, "(A, L)") "Calculate only phonon BTE: ", self%onlyphbte
          if(self%onlyebte) write(*, "(A, L)") "Calculate only electron BTE: ", self%onlyebte
          write(*, "(A, L)") "Include drag: ", self%drag
       end if
       write(*, "(A, L)") "Use tetrahedron method: ", self%tetrahedra
       write(*, "(A, L)") "Plot quantities along path: ", self%plot_along_path
       write(*, "(A, I5)") "Maximum number of BTE/Migdal-Eliashberg equations iterations = ", self%maxiter
       write(*, "(A, 1E16.8)") "BTE/Migdal-Eliashberg equations convergence threshold = ", self%conv_thres
       write(*, "(A, L)") "B-field is on: ", self%Bfield_on
       if(self%Bfield_on) write(*, "(A, 3E16.8)") "B-field  = ", self%Bfield
    end if
    sync all
  end subroutine read_input_and_setup

  subroutine create_chempot_dirs(self, chempot)
    !! Subroutine to create data dump directory tagged by the chemical potential
    !! and subdirectories within.
    
    class(numerics), intent(inout) :: self
    real(r64), intent(in) :: chempot

    !Local variables
    character(len = 1024) :: tag

    !Create chemical potential dependent data dump directory
    write(tag, "(E14.8)") chempot
    self%datadumpdir_T_chempot = trim(adjustl(self%datadumpdir_T)) // '/mu' // trim(adjustl(tag))

    !Create e-ph and ph-e transition probability data directories
    self%Xdir = trim(adjustl(self%datadumpdir_T_chempot)) // '/X'
    self%Ydir = trim(adjustl(self%datadumpdir_T_chempot)) // '/Y'

    !Create chemical potential and T-dependent dielectric data directories
    self%epsilondir = trim(adjustl(self%datadumpdir_T_chempot)) // '/epsilon'
    
    if(this_image() == 1) then
       call system('mkdir -p ' // trim(adjustl(self%datadumpdir_T_chempot)))
       call system('mkdir -p ' // trim(adjustl(self%Xdir)))
       call system('mkdir -p ' // trim(adjustl(self%Ydir)))
    end if
    sync all
  end subroutine create_chempot_dirs
end module numerics_module
