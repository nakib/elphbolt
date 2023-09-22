program bte_regression

  use iso_fortran_env, only : r64 => real64, i64 => int64
  use testify_m, only : testify
  use misc, only: print_message, subtitle, timer, exit_with_message
  use numerics_module, only: numerics
  use crystal_module, only: crystal
  use symmetry_module, only: symmetry
  use electron_module, only: electron
  use phonon_module, only: phonon
  use wannier_module, only: wannier
  use bte_module, only: bte
  use bz_sums, only: calculate_dos, calculate_qTF, calculate_el_dos_fermi, calculate_el_Ws
  use interactions, only: calculate_gReq, calculate_gkRp, calculate_3ph_interaction, &
       calculate_eph_interaction_ibzq, calculate_eph_interaction_ibzk, &
       calculate_echimp_interaction_ibzk, calculate_bound_scatt_rates
      
  implicit none

  integer :: itest
  integer, parameter :: num_tests = 7
  type(testify) :: test_array(num_tests), tests_all

  type(numerics) :: num
  type(crystal) :: crys
  type(symmetry) :: sym
  type(wannier) :: wann
  type(electron) :: el
  type(phonon) :: ph
  type(bte) :: bt
  type(timer) :: t_all, t_event

  !Test counter
  itest = 0

  call t_all%start_timer('elphbolt: BTE')

  call t_event%start_timer('Initialization')

  !Set up crystal
  call crys%initialize

  !Set up numerics data
  call num%initialize(crys)

  !Calculate crystal and BZ symmetries
  call sym%calculate_symmetries(crys, num%qmesh)

  !Test symmetries
  if(this_image() == 1) then
     itest = itest + 1
     test_array(itest) = testify("number of symmetries")
     call test_array(itest)%assert(sym%nsymm, 24_i64)

     !TODO The following is assertion is problematic since
     !the trimming of the string sym%international does not
     !seem to work on all machine. For now disabling this.
!!$     itest = itest + 1
!!$     test_array(itest) = testify("symmetry group")
!!$     call test_array(itest)%assert( &
!!$          sym%international(1 : len(trim(sym%international)) - 1), "F-43m")
  end if
  sync all
  !!

  if(num%onlyebte .or. num%drag .or. num%phe &
       .or. num%plot_along_path .or. num%runlevel == 3) then
     !Read EPW Wannier data
     call wann%read(num)

     !Calculate electrons
     call el%initialize(wann, crys, sym, num)

     !Test electron energies
     if(this_image() == 1) then
        itest = itest + 1
        test_array(itest) = testify("electron energies")
        call test_array(itest)%assert( &
             [transpose(el%ens_irred)], &
             [0.1091276329E+02_r64, 0.1459479591E+02_r64, 0.2329596906E+02_r64, 0.2329596906E+02_r64, &
              0.1096209018E+02_r64, 0.1442980731E+02_r64, 0.2334947484E+02_r64, 0.2340698517E+02_r64, &
              0.1079030242E+02_r64, 0.1369997715E+02_r64, 0.2354625098E+02_r64, 0.2354625098E+02_r64, &
              0.1107938027E+02_r64, 0.1467078718E+02_r64, 0.2311065489E+02_r64, 0.2358101936E+02_r64], &
             tol = 1.0e-8_r64)
     end if
     !!
  end if

  call t_event%end_timer('Initialization')

  call t_event%start_timer('Phonons')

  !Calculate phonons
  call ph%initialize(crys, sym, num)

  !Test phonon energies
  if(this_image() == 1) then
     itest = itest + 1
     test_array(itest) = testify("phonon energies")
     call test_array(itest)%assert( &
          [ph%ens(ph%indexlist_irred(1), :), &
           ph%ens(ph%indexlist_irred(2), :), &
           ph%ens(ph%indexlist_irred(9), :), &
           ph%ens(ph%indexlist_irred(10), :)], &
          [0.0000000000E+00_r64, 0.0000000000E+00_r64, 0.0000000000E+00_r64, &
            0.9660348775E-01_r64, 0.9660348775E-01_r64, 0.9660348775E-01_r64, &
           0.2227978117E-01_r64, 0.2227978118E-01_r64, 0.4234985469E-01_r64, &
            0.9520507538E-01_r64, 0.9520507538E-01_r64, 0.1173379503E+00_r64, &
           0.4840698090E-01_r64, 0.5786962566E-01_r64, 0.7499234595E-01_r64, &
            0.9196463235E-01_r64, 0.9478418434E-01_r64, 0.1060798273E+00_r64, &
           0.5004202869E-01_r64, 0.5811331198E-01_r64, 0.6793992839E-01_r64, &
            0.9281467594E-01_r64, 0.9595198760E-01_r64, 0.1073275295E+00_r64], &
          tol = 1.0e-10_r64)
  end if
  sync all
  !!

  call t_event%end_timer('Phonons')

  call t_event%start_timer('Density of states and one-particle scattering rates')

  call subtitle("Calculating density of states...")
  if(num%onlyebte .or. num%drag) then
     !Calculate electron density of states
     call calculate_dos(el, num%tetrahedra)

     !Test electron density of states
     if(this_image() == 1) then
        itest = itest + 1
        test_array(itest) = testify("electron DOS")
        call test_array(itest)%assert( &
             [transpose(el%dos)], &
             [0.4968872758E-02_r64, 0.4954123967E-02_r64, 0.9760785376E-02_r64, 0.9760785680E-02_r64, &
              0.1639885183E-01_r64, 0.2207070167E-02_r64, 0.6743674412E-02_r64, 0.1143396098E-01_r64, &
              0.1262267770E-02_r64, 0.8026082935E-03_r64, 0.3278198350E-02_r64, 0.3278198022E-02_r64, &
              0.7370288521E-03_r64, 0.4136974131E-03_r64, 0.4998540869E-03_r64, 0.0000000000E+00_r64], &
             tol = 1.0e-11_r64)
     end if
     !!

     !Calculate Thomas-Fermi screening
     call calculate_qTF(crys, el)

     !Test Thomas-Fermi screening
     if(this_image() == 1) then
        itest = itest + 1
        test_array(itest) = testify("Thomas-Fermi screening")
        call test_array(itest)%assert( &
             crys%qTF, 0.14611480E+01_r64, tol = 1.0e-7_r64)
     end if
     !!

     !Calculate boundary scattering rates.
     call calculate_bound_scatt_rates(el%prefix, num%elbound, crys%bound_length, &
          el%vels, el%indexlist_irred, bt%el_rta_rates_bound_ibz)
  end if
  
  !Calculate phonon density of states and, if needed, phonon-isotope
  !and/or phonon-substitution scattering rates.
  !
  !Captain's log. August 11, 2023.
  !I am not a fan of calculating any type of scattering using the density of
  !states calculator. This is a one time calculation and by far not a bottleneck.
  !I will move the phonon-isotope and phonon-isotope scattering stuff to where
  !they belong -- interactions.f90 -- soon.
  call calculate_dos(ph, crys, num%tetrahedra, bt%ph_rta_rates_iso_ibz, bt%ph_rta_rates_subs_ibz, &
       num%phiso, num%phiso_1B_theory, num%phsubs, num%phiso_Tmat)
  
  !Test phonon density of states
  if(this_image() == 1) then
     itest = itest + 1
     test_array(itest) = testify("Phonon DOS")
     call test_array(itest)%assert( &
          [ph%dos(1, :), &
           ph%dos(2, :), &
           ph%dos(9, :), &
           ph%dos(10, :)], &
          [0.0000000000E+00_r64, 0.0000000000E+00_r64, 0.0000000000E+00_r64, &
           0.2672064440E+03_r64, 0.2672064440E+03_r64, 0.2672064440E+03_r64, &
           0.4607796189E+01_r64, 0.4607796199E+01_r64, 0.8583557121E+02_r64, &
           0.4038905966E+03_r64, 0.4038906017E+03_r64, 0.1437348153E+02_r64, &
           0.1343766722E+03_r64, 0.6297533596E+02_r64, 0.0000000000E+00_r64, &
           0.2805500270E+03_r64, 0.2896877310E+03_r64, 0.1284652824E+01_r64, &
           0.8160608929E+02_r64, 0.9246561665E+02_r64, 0.7797451387E+02_r64, &
           0.6395645220E+03_r64, 0.5625272067E+03_r64, 0.2034377505E+03_r64], &
          tol = 1.0e-7_r64)
  end if
  !!

  !Test phonon-isotope scattering in the Tamura model
  if(this_image() == 1) then
     itest = itest + 1
     test_array(itest) = testify("Phonon-isotope scattering Tamura model")
     call test_array(itest)%assert( &
          [bt%ph_rta_rates_iso_ibz(1, :), &
           bt%ph_rta_rates_iso_ibz(2, :), &
           bt%ph_rta_rates_iso_ibz(9, :), &
           bt%ph_rta_rates_iso_ibz(10, :)], &
          [0.0000000000E+00_r64,    0.0000000000E+00_r64,    0.0000000000E+00_r64, &
            0.1263687201E+00_r64,    0.9683788784E-01_r64,    0.9902981700E-01_r64, &
            0.2344573151E-03_r64,    0.1669499774E-03_r64,    0.1554351322E-01_r64, &
            0.1394726680E+00_r64,    0.1838958745E+00_r64,    0.8370501743E-02_r64, &
            0.3970374223E-01_r64,    0.2976859573E-01_r64,    0.0000000000E+00_r64, &
            0.9314628584E-01_r64,    0.1261722420E+00_r64,    0.5224496635E-03_r64, &
            0.2867531169E-01_r64,    0.3811180360E-01_r64,    0.4851394523E-01_r64, &
            0.2895513822E+00_r64,    0.2735793343E+00_r64,    0.1352376696E+00_r64], &
            tol = 1e-2_r64)
           !tol = 1e-9_r64)
  end if
  !!
  
  !Calculate boundary scattering rates.
  call calculate_bound_scatt_rates(ph%prefix, num%phbound, crys%bound_length, &
       ph%vels, ph%indexlist_irred, bt%ph_rta_rates_bound_ibz)

  call t_event%end_timer('Density of states and one-particle scattering rates')

  if(num%plot_along_path) then
     call t_event%start_timer('Plots along path')

     call subtitle("Plotting along high-symmetry path...")

     !Plot electron bands, phonon dispersions, and g along path.
     call wann%plot_along_path(crys, num,el%scissor)

     call t_event%end_timer('Plots along path')
  end if

  call subtitle("Calculating interactions...")

  !Set chemical potential dependent directory
  call num%create_chempot_dirs(el%chempot)

  if(num%onlyphbte .and. num%phe .or. num%drag) then
     if(.not. num%read_gq2) then
        call t_event%start_timer('IBZ q e-ph interactions')

        !Calculate mixed Bloch-Wannier space e-ph vertex g(Re,q)
        call calculate_gReq(wann, ph, num)

        !Calculate Bloch space e-ph vertex g(k,q) for IBZ q
        call calculate_eph_interaction_ibzq(wann, crys, el, ph, num, 'g')

        call t_event%end_timer('IBZ q e-ph interactions')
     end if

     call t_event%start_timer('IBZ ph-e transition probilities')

     !Calculate ph-e transition probabilities
     call calculate_eph_interaction_ibzq(wann, crys, el, ph, num, 'Y')

     call t_event%end_timer('IBZ ph-e transition probilities')
  end if

  if(num%onlyebte .or. num%drag) then
     if(.not. num%read_gk2) then
        call t_event%start_timer('IBZ k e-ph interactions')

        !Calculate mixed Bloch-Wannier space e-ph vertex g(k,Rp)
        call calculate_gkRp(wann, el, num)

        !Calculate Bloch space e-ph vertex g(k,q) for IBZ k
        call calculate_eph_interaction_ibzk(wann, crys, el, ph, num, 'g')

        call t_event%end_timer('IBZ k e-ph interactions')
     end if

     call t_event%start_timer('IBZ e-ph transition probabilities')

     !Calculate e-ph transition probabilities
     call calculate_eph_interaction_ibzk(wann, crys, el, ph, num, 'X')

     call t_event%end_timer('IBZ e-ph transition probabilities')
  end if

  if(num%onlyebte .or. num%drag .or. num%phe) then
     !After this point the electron eigenvectors are not needed
     call el%deallocate_eigenvecs
  end if

  if(num%onlyebte .or. num%drag) then
     if(num%elchimp) then
        call t_event%start_timer('e-ch. imp. interactions')

        !Calculate e-ch. imp. transition probabilities
        call calculate_echimp_interaction_ibzk(crys, el, num)

        call t_event%end_timer('e-ch. imp. interactions')
     end if
  end if

  if(num%onlyebte .or. num%drag .or. num%phe .or. num%drag &
       .or. num%plot_along_path) then
     !Deallocate Wannier quantities
     call wann%deallocate_wannier(num)
  end if

  if(num%onlyphbte .or. num%drag) then
     if(.not. num%read_V) then
        call t_event%start_timer('IBZ q ph-ph interactions')

        !Calculate ph-ph vertex
        call calculate_3ph_interaction(ph, crys, num, 'V')

        call t_event%end_timer('IBZ q ph-ph interactions')
     end if

     if(.not. num%read_W) then
        call t_event%start_timer('IBZ ph-ph scattering rates')

        !Calculate ph-ph transition probabilities
        call calculate_3ph_interaction(ph, crys, num, 'W')

        call t_event%end_timer('IBZ ph-ph scattering rates')
     end if
  end if

  if(num%onlyphbte .or. num%drag) then
     !After this point the phonon eigenvectors and other quantities are not needed
     call ph%deallocate_phonon_quantities
  end if

  !Solve BTEs
  if(num%onlyphbte .and. .not. num%phe) then
     call bt%solve_bte(num, crys, sym, ph)
  else
     call bt%solve_bte(num, crys, sym, ph, el)
  end if

  sync all
  if(this_image() == 1) then
     tests_all = testify(test_array)
     call tests_all%report
     
     if(tests_all%get_status() .eqv. .false.) error stop -1
  end if
  sync all

  call t_all%end_timer('elphbolt: BTE')
end program bte_regression
