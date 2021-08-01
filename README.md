![](./logo/elphbolt_logo.png)

# elphbolt

`elphbolt` (short for electron-phonon Boltzmann transport) is a modern
Fortran (2018 standard) code for solving the coupled electron and phonon
Boltzmann transport equations (BTEs). It is a "free as in freedom" code
distributed under the GNU Public License (GPL) version 3. You can read
more about the philosophy of software freedom here:
<https://www.gnu.org/philosophy/free-sw.en.html>.

Using **ab initio** electron-phonon and phonon-phonon interactions and a
fully wave vector and electron band/phonon branch resolved formulation
of the BTEs, `elphbolt` can calculate the

  - phonon and electronic thermal conductivities;
  - electronic conductivity;
  - phonon and electronic contributions to the thermopower; and
  - effect of the mutual electron-phonon drag on the transport
    coefficients listed above.

Stylistically, it is designed to be simple, small, fast, and extensible.
Object oriented programming concepts are combined with the procedural
style, which allows fast development, while resulting in a rather
compact code. The symmetries of the crystal are fully exploited and the
transport active Fermi window is used to allow the sampling of extremely
fine wave vector meshes needed for an accurate solution of the BTEs.
Parallelism is achieved through modern Fortran's intrinsic `coarrays`
feature that is fully supported by recent versions of both the `gcc` and
`intel` compilers.

`elphbolt` currently interfaces with the `Quantum Espresso` suite for
the phonon quantities and the `EPW` code for the Wannier space
information.

The project is currently in the beta phase and being used to do
publication level calculations.

## Installation on HPC systems with `gcc`

### 1\. Get `spglib`

`elphbolt` uses the spglib (<https://spglib.github.io/spglib/>) library
for crystal symmetry analysis. Currently, versions 1.6.0 and 1.14.1 have
been tried and tested. A lot of HPCs provide this library as a module,
so check before building from source.

### 2\. Get `OpenCoarrays`

`OpenCoarrays` (<http://www.opencoarrays.org>) is an implementation of
the `coarrays` functionalities. Some HPC systems include it (for example
MareNostrum4 and, more recently, LaPalma of the Barcelona Supercomputing
Center (BSC)). In that case, load it. If not, you can build it from
source. I was able to build version 2.0.0 on the LaPalma system of BSC
by first loading `gcc`, `openmpi`, and `cmake`, then going into the
`OpenCoarrays-2.0.0` directory and saying `./install.sh`. The
executables `caf` and `cafrun` were then installed in
`~/OpenCoarrays-2.0.0/prerequisites/installations/bin/`. I added the
above directory to my `$PATH` in my `bashrc`.

On the Sirius cluster of Boston College, I was able to build version
2.9.2 by first loading `gnu_gcc/9.2.0`, `openmpi/4.0.5gcc9.2.0`, and
`cmake/3.18.4`, then going into the `OpenCoarrays-2.9.2` directory and
saying `./install.sh --with-fortran
/usr/public/gnu_gcc/9.2.0/bin/gfortran --with-cxx
/usr/public/gnu_gcc/9.2.0/bin/g++ --with-c
/usr/public/gnu_gcc/9.2.0/bin/gcc --with-mpi
/usr/public/openmpi/4.0.5gcc.9.2.0/ --with-cmake
/usr/public/cmake/3.18.4/bin/cmake`. The executables `caf` and `cafrun`
were then installed in
`~/OpenCoarrays-2.9.2/prerequisites/installations/opencoarrays/2.9.2/bin/`.
I added the above directory to my `$PATH` in my `bashrc`.

(You can change the build path with the additional argument
`--install-prefix <path to build directory>` to `install.sh`.)

### 3\. Create your `<title>.make` file

In the `elphbolt` src directory you will find a few `<title>.make`
files. Adapt them to your own architecture. Then, in the `Makefile`
include your personal `<title>.make` file at the top, commenting out any
exiting one.

Once you have done the above, simply say `make`. This should build the
executable `elphbolt.x` one directory above.

## Installation on HPC systems with `intel`

\[Working on it.\]

## Tests

A full example for cubic silicon is provided.

## Workflow

This is a transport code. And it comes after doing some DFT, DFPT, and
Wannier calculations. The following quantities are the inputs for a
calculation with `elphbolt`:

### Input file

The input file - `input.nml` - contains the information about the
crystal and the various parameters of the calculation. A full
description of all the input parameters will be provided soon in the
documentation. I will also provide example input files for the test
materials.

### Second order interatomic force constants

This comes out of the usual `ph.x` and `q2r.x` calculation from `Quantum
Espresso`. This file is needed to calculate phonon quantities and must
be named `espresso.ifc2`.

### Third order interatomic force constants

This file, which must be named `FORCE_CONSTANTS_3RD`, is needed to
calculate the 3-ph scattering rates. This is a required file if you seek
a solution of the decoupled phonon BTE or the coupled electron-phonon
BTEs.

This must be provided for a solution to the phonon BTE or the coupled
electron-phonon BTEs. See documentation for the code `thirdorder.py`
(<https://bitbucket.org/sousaw/thirdorder/src/master>) for how to
generate this file.

### Wannier space information

These are required if you want to solve a decoupled electron BTE,
include phonon-electron interaction in the decoupled phonon BTE, or
solve the coupled electron-phonon BTEs.

These include the files `rcells_k`, `rcells_q`, `rcells_g`, `wsdeg_k`,
`wsdeg_q`, and `wsdeg_g` which must be printed out of an `EPW`
calculation. I will provide a patched `EPW/src/ephwann_shuffle.f90` code
which will print these quantities out during `EPW`'s Bloch -\> Wannier
calculation step.

We will also need the files `epmatwp1` and `epwdata.fmt`, both of which
are outputted by `EPW` after the Bloch -\> Wannier calculation step. The
first contains the Wannier space electron-phonon matrix elements and the
second contains the Wannier space dynamical matrix and Hamiltonian.

### High symmetry electron and phonon wave vector path and initial electron wave vector

You need to provide a wave vector path file named `highsympath.txt` (to
be used as both the electron and phonon wave vectors) and an initial
electron wave vector file named `initialk.txt` if you want the electron
bands, phonon dispersions, and electron-phonon matrix elements
calculated along the path.

## Description of `input.nml`

There are 5 Namelists in the `input.nml` file: `allocations`,
`crystal_info`, `electrons`, `numerics`, and `wannier`. Users of the
`ShengBTE` code will find the format of this file familiar. Below the
keys for each Namelist are described.

### `allocations`

| key           | Type    | Default | Description                     |
| ------------- | ------- | ------- | ------------------------------- |
| `numelements` | Integer | 0       | Number of types of basis atoms. |
| `numatoms`    | Integer | 0       | Number of basis atoms.          |

### `crystal_info`

| key               | Type                                  | Default   | Description                                                                                                                                                                                                                                |
| ----------------- | ------------------------------------- | --------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `name`            | String                                | "Crystal" | Name of material.                                                                                                                                                                                                                          |
| `elements`        | String array of size `numelements`    | 'X'       | Elements in the basis.                                                                                                                                                                                                                     |
| `atomtypes`       | Integer array of size `numatoms`      | 0         | Integer tagging unique elements in the basis.                                                                                                                                                                                              |
| `masses`          | Real array of size `numelements`      | \-1.0     | Masses of the basis atoms in amu. If masses are not provided, set `autoisotopes` to .True..                                                                                                                                                |
| `autoisotopes`    | Logical                               | .True.    | Use isotopic mix for masses?                                                                                                                                                                                                               |
| `lattvecs`        | 3 x 3 real array                      | 0.0       | Lattice vectors in Cartesian coordinates in units of nm. If `twod` is .True., the crystal must be positioned on the x-y plane and the third lattice vector must be of the form (0 0 layer thickness).                                      |
| `basis`           | 3 x `numatoms` real array             | 0.0       | Atomic basis vectors in crystal coordinates (i.e. fraction of `lattvecs`).                                                                                                                                                                 |
| `polar`           | Logical                               | .False.   | Is the system polar?                                                                                                                                                                                                                       |
| `born`            | 3 x 3 x `numatoms` rank-3 real tensor | 0.0       | Born effective charge tensor (from phonon calculation).                                                                                                                                                                                    |
| `epsilon`         | 3 x 3 rank-2 real tensor              | 0.0       | High-frequency dielectric tensor (from phonon calculation).                                                                                                                                                                                |
| `read_epsiloninf` | Real                                  | .False.   | Read high-frequency dielectric constant from input?                                                                                                                                                                                        |
| `epsiloninf`      | Real                                  | 0.0       | High-frequency scalar dielectric constant. If `read_epsiloninf` is .True. (.False.), this is read from the input (set equal to the trace-average of `epsilon`). Currently this quantity is not used in any calculation.                    |
| `epsilon0`        | Real                                  | 0.0       | Static scalar dielectric constant. Used for screening electron-charged impurity interaction, if included. Look up `elchimp` under the Namelist `numerics`. For the default value of `epsilon0`, the electron-charged interaction blows up. |
| `T`               | Real                                  | \-1.0\_dp | Crystal temperature in K.                                                                                                                                                                                                                  |
| `twod`            | Logical                               | .False.   | Is the system (quasi)-2-dimensional? See description of `lattvecs` also.                                                                                                                                                                   |
| `subs_masses`     | Real array of size `numelements`      | 0.0       | Masses of substitution atoms in amu. This is needed if `phsubs` is .True. See table of keys for Namelist `numerics`.                                                                                                                       |
| `subs_conc`       | Real array of size `numelements`      | 0.0       | Concentration of the substitutional atoms in cm<sup>-3</sup> (or cm<sup>-2</sup> if `twod` is .True.). This is needed if `phsubs` is .True. See table of keys for Namelist `numerics`.                                                     |

### `electrons`

| key                | Type                         | Default        | Description                                                                                                                                                               |
| ------------------ | ---------------------------- | -------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `spindeg`          | Integer                      | 2              | Spin degeneracy of the bands.                                                                                                                                             |
| `enref`            | Real                         | \-999999.99999 | Electron referenc energy in eV. This is the center of the transport active window. Also see description for `fsthick`. See table of keys for Namelist 'numerics'.         |
| `chempot`          | Real                         | \-999999.99999 | Chemical potential in eV.                                                                                                                                                 |
| `metallic`         | Logical                      | .False.        | Is the system metallic?                                                                                                                                                   |
| `numbands`         | Integer                      | 0              | Total number of electronic Wannier bands.                                                                                                                                 |
| `indlowband`       | Integer                      | 0              | Lowest transport band index.                                                                                                                                              |
| `indhighband`      | Integer                      | 0              | Highest transport band index.                                                                                                                                             |
| `indlowconduction` | Integer                      | 0              | Lowest conduction band index. For `metallic` .True., this or `indhighvalence` must be provided.                                                                           |
| `indhighvalence`   | Integer                      | 0              | Highest valence band index. For `metallic` .True., this or `indlowconduction` must be provided.                                                                           |
| `dopingtype`       | Character                    | 'x'            | Type of doping ('n' or 'p'). This is needed for `runlevel` 0 only. See table of keys for Namelist 'numerics'.                                                             |
| `numconc`          | Integer                      | 100            | Number of carrier concentration points. This is needed for `runlevel` 0 only. See table of keys for Namelist 'numerics'.                                                  |
| `conclist`         | Real array of size `numconc` | 0.0            | List carrier concentrations in cm<sup>-3</sup> (or cm<sup>-2</sup> if `twod` is .True.). This is needed for `runlevel` 0 only. See table of keys for Namelist 'numerics'. |
| `numT`             | Integer                      | 100            | Number of temperature points. This is needed for `runlevel` 0 only. See table of keys for Namelist 'numerics'.                                                            |
| `Tlist`            | Real array of size `numT`    | 100            | List of temperatures in K. This is needed for `runlevel` 0 only. See table of keys for Namelist 'numerics'.                                                               |
| `Zn`               | Real                         | 0.0            | Ionization number if donor impurities. This is needed only when `elchimp` is .True. and `metallic` is .False. See table of keys for Namelist 'numerics'.                  |
| `Zp`               | Real                         | 0.0            | Ionization number if acceptor impurities. This is needed only when `elchimp` is .True. and `metallic` is .False. See table of keys for Namelist 'numerics'.               |

### `numerics`

| key               | Type                    | Default | Description                                                                                                                                                                                                                                                                                                                                                   |
| ----------------- | ----------------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `qmesh`           | Integer array of size 3 | 1 1 1   | Phonon wave vector mesh (q).                                                                                                                                                                                                                                                                                                                                  |
| `mesh_ref`        | Integer                 | 1       | Electron wave vector mesh (k) refinement factor with respect to the phonon mesh.                                                                                                                                                                                                                                                                              |
| `fsthick`         | Real                    | 0.0     | Fermi surface thickness in eV.                                                                                                                                                                                                                                                                                                                                |
| `datadumpdir`     | String                  | "./"    | Runtime data dump directory.                                                                                                                                                                                                                                                                                                                                  |
| `read_gq2`        | Logical                 | .False. | Read electron-phonon (irreducible wedge q) vertices from disk?                                                                                                                                                                                                                                                                                                |
| `read_gk2`        | Logical                 | .False. | Read electron-phonon (irreducible wedge k) verticesfrom disk?                                                                                                                                                                                                                                                                                                 |
| `read_V`          | Logical                 | .False. | Read phonon-phonon (irreducible wedge q) vertices from disk?                                                                                                                                                                                                                                                                                                  |
| `read_W`          | Logical                 | .False. | Read phonon-phonon (irreducible wedge q) transition probabilities from disk?                                                                                                                                                                                                                                                                                  |
| `tetrahedra`      | Logical                 | .False. | Use the analytic tetrahedron method intead of the triangular method for 3d delta function evaluation?                                                                                                                                                                                                                                                         |
| `phe`             | Logical                 | .False. | Include phonon-electron interaction in phonon BTE?                                                                                                                                                                                                                                                                                                            |
| `phiso`           | Logical                 | .False. | Include phonon-isotope interaction in phonon BTE?                                                                                                                                                                                                                                                                                                             |
| `phsubs`          | Logical                 | .False. | Include phonon-substitution interaction in phonon BTE? If .True., look up `subs_masses` and `subs_conc` under the Namelist `crystal_info`.                                                                                                                                                                                                                    |
| `onlyphbte`       | Logical                 | .False. | Calculate phonon BTE without electron drag?                                                                                                                                                                                                                                                                                                                   |
| `elchimp`         | Logical                 | .False. | Include electron-charged impurity scattering in electron BTE? If .True., look up `epsilon0` under Namelist `crystal_info` and `Zn` and `Zp` under Namelist `electrons`.                                                                                                                                                                                       |
| `onlyebte`        | Logical                 | .False. | Calculate electron BTE without phonon drag?                                                                                                                                                                                                                                                                                                                   |
| `drag`            | Logical                 | .True.  | Include electron and phonon drag term in the phonon and electron BTE, respectively.                                                                                                                                                                                                                                                                           |
| `maxiter`         | Intger                  | 50      | Maximum number of iteration steps for the BTE(s).                                                                                                                                                                                                                                                                                                             |
| `conv_thres`      | Real                    | 1e-4    | Convergence threshold for the BTE(s).                                                                                                                                                                                                                                                                                                                         |
| `runlevel`        | Integer                 | 1       | Control for the type of calculation. 0: Calculate table of chemical potentials for a given doping type, temperature range, and carrier concentrations. Look up `dopingtype`, `numconc`, `conclist`, `numT`, and `Tlist` under Namelist `electrons`. 1: Transport calculation(s). 2: Post-processing results to calculate the spectral transport coefficients. |
| `plot_along_path` | Logical                 | .False. | Plot Wannier interpolated quantities along high symmetry wave vectors?                                                                                                                                                                                                                                                                                        |
| `ph_en_min`       | Real                    | 0.0     | Lower bound of equidistant phonon energy mesh in eV. Only needed for `runlevel` 2.                                                                                                                                                                                                                                                                            |
| `ph_en_max`       | Real                    | 1.0     | Upper bound of equidistant phonon energy mesh in eV. Only needed for `runlevel` 2.                                                                                                                                                                                                                                                                            |
| `ph_en_num`       | Integer                 | 100     | Number of equidistant phonon energy mesh points. Only needed for `runlevel` 2.                                                                                                                                                                                                                                                                                |
| `el_en_min`       | Real                    | \-10.0  | Lower bound of equidistant electron energy mesh in eV. Only needed for `runlevel` 2.                                                                                                                                                                                                                                                                          |
| `el_en_max`       | Real                    | 10.0    | Upper bound of equidistant electron energy mesh in eV. Only needed for `runlevel` 2.                                                                                                                                                                                                                                                                          |
| `el_en_num`       | Integer                 | 100     | Number of equidistant electron energy mesh points. Only needed for `runlevel` 2.                                                                                                                                                                                                                                                                              |

### `wannier`

| key            | Type                    | Default | Description                                                                                                                                               |
| -------------- | ----------------------- | ------- | --------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `coarse_qmesh` | Integer array of size 3 | 0 0 0   | Coarse phonon wave vector mesh employed in the Wannier calculation. This must match the q-mesh in the Quantum Espresso second order force constants file. |

## Description of output files

The code produces a large amount of data. Here, we provide a description
of the various types output files.

Below I(F)BZ = irreducible (full) Brillouin zone; RTA = relaxation time
approximation; ch. imp. = charged impurities; `numbranches` = number of
phonon branches.

### Zero temperature data

| File name                  | Directory         | Units                                        | Description                                                                                                                                    |
| -------------------------- | ----------------- | -------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| `gk2.istate*`              | `datadumpdir/g2/` | eV<sup>2</sup>                               | Squared e-ph (1-phonon) vertices for every IBZ electron state. Binary.                                                                         |
| `gq2.istate*`              | `datadumpdir/g2/` | eV<sup>2</sup>                               | Squared e-ph (1-phonon) vertices for every IBZ electron state. Binary.                                                                         |
| `Vm2.istate*`              | `datadumpdir/V2/` | eV<sup>2</sup>Å<sup>-6</sup>amu<sup>-3</sup> | Squared ph-ph (3-phonon) vertices for every IBZ phonon state. Binary.                                                                          |
| `el(ph).dos`               | `./`              | eV<sup>-1</sup>                              | Band resolved electronic (phononic) density of states. `numbands` (`numbranches`) columns of reals.                                            |
| `el(ph).ens_ibz`           | `./`              | eV                                           | IBZ electronic (phononic) band energies. `numbands` (`numbranches`) columns of reals.                                                          |
| `el.inwindow_states_ibz`   | `./`              | none                                         | IBZ electronic states (wave vector index, band index) within the transport active window. 2 columns of integers.                               |
| `el(ph).vels_ibz`          | `./`              | Kms<sup>-1</sup>                             | IBZ electronic (phononic) band (branch) velocities. In each row, there are 3 (Cartesian direction) sets of `numbands` (`numbranches`) numbers. |
| `el(ph).wavevecs_ibz[fbz]` | `./`              | crystal                                      | IBZ \[FBZ\] electronic (phononic) wave vectors. For the electrons, these are only within the transport window.                                 |
| `ph.W_rta_phiso[subs]`     | `./`              | THz                                          | IBZ RTA ph-iso \[subs\] scattering rates. `numbranches` columns of reals.                                                                      |
| `el.ens_kpath`             | `./`              | eV                                           | Electron energies along the given k-path.                                                                                                      |
| `ph.ens_qpath`             | `./`              | eV                                           | Phonon energies along the given q-path.                                                                                                        |
| `gk_qpath`                 | `./`              | eV                                           | Absolute value of the e-ph matrix elements (averaged over the degenerate bands and branches) for the given k-vector and q-path.                |

### Finite temperature data

| File name                          | Directory            | Units                          | Description                                                                                                                                            |
| ---------------------------------- | -------------------- | ------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `Xchimp.istate*`                   | `datadumpdir/mu*/X/` | THz                            | Transition probability for e-ch. imp. processes for every IBZ electron state. Binary.                                                                  |
| `Xminus[plus].istate*`             | `datadumpdir/mu*/X/` | THz                            | Transition probability for e-ph (1-phonon) minus \[plus\] processes for every IBZ electron state. Binary.                                              |
| `Y.istate*`                        | `datadumpdir/mu*/Y/` | THz                            | Transition probability for ph-e (1-phonon) processes for every IBZ phonon state. Binary.                                                               |
| `Wm[p].istate*`                    | `datadumpdir/T*/W/`  | THz                            | Transition probability for ph-ph (3-phonon) minus \[plus\] processes for every IBZ phonon state. Binary.                                               |
| `el.W_rta_eph[chimp]`              | `./T*/`              | THz                            | IBZ RTA el-ph \[ch. imp.\] scattering rates. `numbands` columns of reals. Identically zero for bands outside the transport window.                     |
| `ph.W_rta_3ph[phe]`                | `./T*/`              | THz                            | IBZ RTA ph-ph \[e\] scattering rates. `numbranches` columns of reals.                                                                                  |
| `drag[nodrag]_el_sigma_*`          | `./T*/`              | *Ω*<sup>-1</sup>m<sup>-1</sup> | Band resolved (`_<integer>`) and total (`_tot`) charge conductivity tensor at every iteration step.                                                    |
| `drag[nodrag]_el_alphabyT_*`       | `./T*/`              | Am<sup>-1</sup>K<sup>-1</sup>  | Band resolved (`_<integer>`) and total (`_tot`) electronic Peltier(-ish) coefficient tensor at every iteration step.                                   |
| `drag[nodrag]_el_kappa0_*`         | `./T*/`              | Wm<sup>-1</sup>K<sup>-1</sup>  | Band resolved (`_<integer>`) and total (`_tot`) electronic thermal conductivity (zero E-field) tensor at every iteration step.                         |
| `drag[nodrag]_el_sigmaS_*`         | `./T*/`              | Am<sup>-1</sup>K<sup>-1</sup>  | Band resolved (`_<integer>`) and total (`_tot`) electronic thermopower times conductivity tensor at every iteration step.                              |
| `drag_ph_alphabyT_*`               | `./T*/`              | Am<sup>-1</sup>K<sup>-1</sup>  | Branch resolved (`_<integer>`) and total (`_tot`) phonon Peltier(-ish) coefficient tensor at every iteration step.                                     |
| `drag[nodrag]_ph_kappa_*`          | `./T*/`              | Wm<sup>-1</sup>K<sup>-1</sup>  | Branch resolved (`_<integer>`) and total (`_tot`) phonon thermal conductivity tensor at every iteration step.                                          |
| `RTA{nodrag}(partdcpl)[drag]_I0_*` | `./T*/`              | nmeVK<sup>-1</sup>             | Band resolved (`_<integer>`) and total (`_tot`) electronic response function to ∇ T-field in the RTA {dragless} (partially decoupled) \[drag\] theory. |
| `RTA{nodrag}(partdcpl)[drag]_J0_*` | `./T*/`              | nmC                            | Band resolved (`_<integer>`) and total (`_tot`) electronic response function to E-field in the RTA {dragless} (partially decoupled) \[drag\] theory.   |
| `RTA{nodrag}[drag]_F0_*`           | `./T*/`              | nmeVK<sup>-1</sup>             | Branch resolved (`_<integer>`) and total (`_tot`) phononic response function to ∇ T-field in the RTA {dragless} \[fully coupled\] theory.              |
| `drag_G0_*`                        | `./T*/`              | nmC                            | Branch resolved (`_<integer>`) and total (`_tot`) phononic response function to E-field in fully coupled theory.                                       |

### Postprocessing (runlevel 2)

| File name                                                           | Directory | Units                                         | Description                                                                                                                                                                                  |
| ------------------------------------------------------------------- | --------- | --------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `RTA{nodrag}(partdcpl)[drag]_{([iterated_el])}_sigma_spectral_*`    | `./T*/`   | *Ω*<sup>-1</sup>m<sup>-1</sup>eV<sup>-1</sup> | Band resolved (`_<integer>`) and total (`_tot`) spectral charge conductivity tensor in the RTA {(\[iterated\])} {dragless} (partially decoupled) \[drag\] theory.                            |
| `RTA{nodrag}(partdcpl)[drag]_{([iterated_el])}_alphabyT_spectral_*` | `./T*/`   | Am<sup>-1</sup>K<sup>-1</sup>eV<sup>-1</sup>  | Band resolved (`_<integer>`) and total (`_tot`) spectral electronic Peltier(-ish) coefficient tensor in the RTA {(\[iterated\])} {dragless} (partially decoupled) \[drag\] theory.           |
| `RTA{nodrag}(partdcpl)[drag]_{([iterated_el])}_kappa0_spectral_*`   | `./T*/`   | Wm<sup>-1</sup>K<sup>-1</sup>eV<sup>-1</sup>  | Band resolved (`_<integer>`) and total (`_tot`) spectral electronic thermal conductivity (zero E-field) tensor in the RTA {(\[iterated\])} {dragless} (partially decoupled) \[drag\] theory. |
| `RTA{nodrag}(partdcpl)[drag]_{([iterated_el])}_sigmaS_spectral_*`   | `./T*/`   | Am<sup>-1</sup>K<sup>-1</sup>eV<sup>-1</sup>  | Band resolved (`_<integer>`) and total (`_tot`) spectral electronic thermopower times conductivity tensor in the RTA {(\[iterated\])} {dragless} (partially decoupled) \[drag\] theory.      |
| `drag_iterated_ph_alphabyT_spectral_*`                              | `./T*/`   | Am<sup>-1</sup>K<sup>-1</sup>eV<sup>-1</sup>  | Branch resolved (`_<integer>`) and total (`_tot`) spectral phonon Peltier(-ish) coefficient tensor in the iterated drag theory.                                                              |
| `RTA{nodrag}[drag]_{[iterated_ph]}_kappa_spectral_*`                | `./T*/`   | Wm<sup>-1</sup>K<sup>-1</sup>eV<sup>-1</sup>  | Branch resolved (`_<integer>`) and total (`_tot`) spectral phonon thermal conductivity tensor in the RTA {\[iterated\]} {dragless} \[drag\] theory.                                          |
| `el[ph].en_grid`                                                    | `./`      | eV                                            | Uniform electron \[phonon\] energy mesh for spectral coefficient calculation.                                                                                                                |
