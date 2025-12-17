# V3gpu

[![CUDA](https://img.shields.io/badge/CUDA-Fortran-76B900?logo=nvidia)](https://developer.nvidia.com/cuda-fortran)

V3gpu is a **GPU-accelerated implementation** of three-phonon interaction matrix elements (`V3` / `Vm2`) used in lattice dynamics and Boltzmann Transport Equation (BTE) calculations.

The goal is to offload heavy `Vm2` computations from CPU to NVIDIA GPUs using CUDA Fortran, achieving **order-of-magnitude speedups** for large q-meshes while maintaining compatibility with modern phonon transport workflows (e.g., Elphbolt-style solvers).
---

## Key Features

- **CUDA Fortran kernels** for fast 3-phonon interaction vertex evaluation
- **Two-stage computation pipeline**: 
  - CPU-side energy/momentum selection rule filtering
  - GPU-side matrix element calculation for valid processes only
- **Dynamic memory management**: Expandable/shrinkable triplet lists
- **Optimized memory layout**: Flat 1D arrays for maximum GPU bandwidth
- **Energy conservation methods**: Triangle and tetrahedron integration schemes
- **Scalable**: Handles large q-meshes (20³, 28³, 30³+)
- **Test**: Tested on different materials as Si, MgB2, SiC-{2H, 4H, 6H}, Graphene ...
- **Library interface**: Builds as `libgpu_v3.so` for integration into existing Fortran codes

---

## Project Structure

```
V3gpu/
├── src/
│   ├── gpu_3ph_kernels.f90   # CUDA kernels (device functions, Vm2 calculation)
│   ├── gpu_driver.f90        # Host driver (memory management, kernel launch)
│   └── gpu_wrapper.f90       # Public Fortran API (interface for external codes) 
├── tests/
│   └── V3_gpu_test.f90       # test program
├── Makefile                  # Build configuration
├── README.md                 # This file
└── LICENSE                   # Not selected 
```

### What each file does

**`gpu_3ph_kernels.f90`**
- Device functions: `Vm2_3ph_reference_dev` (matrix element calculation)
- Utility functions: `mux_state`, `demux_state`, `mux_vector`, `imod`
- Delta functions: `delta_fn_triang`, `delta_fn_tetra` (energy conservation)
- Main kernel: `calculate_Vm2_kernel` (parallel Vm2 computation)

**`gpu_driver.f90`**
- `compute_V2_on_gpu`: Main public interface
- `build_tripletlist`: CPU-side triplet filtering with MSR/ESR
- Memory management: `expand_int64_2d`, `shrink_int64_2d`
- Dynamic allocation for optimal memory usage

**`gpu_wrapper.f90`**
- Provides a **C-compatible interface** to the Fortran GPU driver (`bind(C)`)
- Returns results via **C pointers** (`c_ptr`)
- Handles conversion between C and Fortran types and logicals 

**`V3_gpu_test.f90`**
- Standalone test program demonstrating usage
- Validates momentum and energy selection rules
- Benchmarking utilities

### Integration with ElphBolt

V3gpu is designed to be called from phonon transport codes.

In the **ElphBolt** codebase, a separate module (`gpu_interface`) calls
`compute_V2_on_gpu_c` from `gpu_wrapper.f90`.  
This interface layer translates ElphBolt-specific data structures into the
C-compatible arrays expected by V3gpu and retrieves the resulting Vm2 values
and triplet lists via pointers.

---

## Requirements

### Software
- **NVIDIA HPC SDK** 25.7-0 version with CUDA Fortran compiler (`nvfortran`)
- **CUDA Toolkit** CUDA Version: 12.4
- **GPU**: NVIDIA GPU with compute capability: 8.9 (Ada)
- **OS**: Linux (tested on Ubuntu 24.04.2 LTS)

### Hardware
- Recommended: >= 20 GB GPU memory for production runs (can push into the ~40^3 range, depending on material complexity).
- Minimum: 4 GB GPU memory for small test cases like 16^3 or smaller

---

## Installation

### 1. Clone the repository
```bash
git clone https://github.com/sallyi95/V3gpu.git
cd V3gpu
```

### 2. Set up environment
```bash
# Load NVIDIA HPC SDK module (adjust for your system)
module load nvhpc

# Verify nvfortran is available
nvfortran --version
```

### 3. Build the library and test program
```bash
make all
```

This creates:
- `libgpu_v3.so` - Shared library
- `V3_gpu_test` - Test executable

### 4. Run tests
```bash
./V3_gpu_test
```

Expected output:
```
==========================================
GPU 3-Phonon Calculation Test
==========================================

Utility Tests:
 - mux_state/demux_state ............. PASSED
 - int_div ........................... PASSED
 - imod (positive remainder) ......... PASSED
 - mux_vector (0-based) .............. PASSED
 - mux_vector (1-based) .............. PASSED

All utility tests passed.

------------------------------------------
S-list Construction (Host)
------------------------------------------
Initial possible transitions:        5,511,240
Valid transitions kept (S_count):       66,672
Filtered out:                        5,444,568
Percentage kept:                        1.21 %
Reduction factor:                        82.66×

CPU time to build S_list:               0.225 s

------------------------------------------
GPU Kernel Launch
------------------------------------------
Blocks:                                 261
Threads per block:                      256
Total threads:                        66,816
Triplets processed:                   66,672

GPU kernel time (Vm2):                  0.612 s
Total compute_V2_on_gpu time:           ~0.839 s

------------------------------------------
Output Statistics
------------------------------------------
V2 array size:                        66,672
Non-zero entries:                     66,672
V2(1):                        2.1359767E+03
V2(100):                      4.7184829E+05
V2(last):                     1.4854808E+06
V2 Norm:                      3.93220947E+08

------------------------------------------
TEST PASSED — No kernel crashes!
All transitions processed successfully.
==========================================

```

---

## Usage

### Basic Example

```fortran
program my_phonon_calc
  use gpu_3ph_driver, only: compute_V2_on_gpu
  
  ! Declare arrays (see V3_gpu_test.f90 for full example)
  complex(real64), allocatable :: evecs(:, :, :)
  real(real64), allocatable :: V2(:)
  integer(int64), allocatable :: triplet_list(:, :)
  
  ! ... Initialize phonon data (eigenvectors, IFCs, etc.) ...
  
  ! Call GPU computation
  call compute_V2_on_gpu( &
    nb, nwv_irred, nwv, ntrip, &
    evecs, Index_i, Index_j, Index_k, ifc3, &
    indexlist_irred, wavevecs, wvmesh, reclatt, Rj, Rk, &
    simplex_map, simplex_count, simplex_evals, ens, use_tetra, &
    V2, triplet_list, triplet_count, istat)
  
  if(istat == 0) then
    print*, "Success! Computed", triplet_count, "matrix elements"
    print*, "V2 norm:", norm2(V2)
  else
    print*, "GPU computation failed with error:", istat
  end if
end program
```

### Input Requirements

| Parameter | Type | Description |
|-----------|------|-------------|
| `evecs` | `complex(real64)` | Phonon eigenvectors `(nwv, nb, nb)` |
| `ifc3` | `real(real64)` | 3rd-order force constants `(3,3,3,ntrip)` |
| `ens` | `real(real64)` | Phonon frequencies `(nwv, nb)` |
| `wavevecs` | `real(real64)` | q-vectors in fractional coords `(nwv, 3)` |
| `wvmesh` | `integer(int64)` | q-mesh dimensions `(3)` |
| `simplex_map/count/evals` | - | Triangle/tetrahedron integration data |

### Output

| Parameter | Type | Description |
|-----------|------|-------------|
| `V2` | `real(real64)` | Squared matrix elements `\|V\|²` for valid processes |
| `triplet_list` | `integer(int64)` | Compact list of valid (lambda1, lambda2 , lambda3) triplets |
| `triplet_count` | `integer(int64)` | Number of valid processes found |

---

## Performance

### Benchmark: Silicon 30³ q-mesh

| Configuration | Processes | Vm2 CPU Time | Vm2 GPU Time | Speedup |
|--------------|-----------|----------|----------|---------|
| Si 10³ (6 bands) | ~ 845,872 valid | ~ 3.405 min | ~ 0.106 min | **32x** |
| Si 15³ (6 bands) | ~ 5.2M valid | ~ 23.31 min | ~ 0.648 min | **36x** |

*Tested on Intel Xeon Gold 6248R + NVIDIA Ada 20GB*

### Memory Usage

For a 25³ mesh with 6 bands:
- **Triplet list**: 1.3864 GB (58M processes)
- **Vm2 array memory**: ~ 0.4621 GB 
- **Total GPU memory**: ~ 4-5 GB

---

## Algorithm Overview

### Two-Stage Pipeline

**Stage 1: CPU-side Filtering (MSR + ESR)**
```
For each possible transition (q1, q2, s1, s2, s3):
  1. Momentum Selection Rule (MSR): q3 = ±(q1 ± q2)
  2. Energy Selection Rule (ESR): delta(omega1 ± omega2 ± omega3) > 0
  3. If satisfied: add to triplet_list
```

**Stage 2: GPU-side Computation**
```
Launch CUDA kernel with N threads (N = valid processes):
  Each thread computes Vm2 for one triplet:
    Vm2 = |∑_t ∑_{αβγ} Φ(a,b,c) e₁ e₂* e₃* exp(iq·R)|²
```

### Energy Conservation

Two methods implemented:
- **Triangle method**: Kurganskii et al. Phys. Stat. Sol.(b) 129, 293 (1985)
- **Tetrahedron method**: Standard linear interpolation scheme

Both evaluate δ(E) for energy conservation checking.

---

## Troubleshooting

### Common Issues

**1. Error 700 (Illegal Memory Access)**
```
Error: Kernel execution failed! Error code: 700
```
**Solution**: Check that `indexlist_irred` correctly maps IBZ → full BZ indices.
Verify with: `print*, "Max IBZ index:", maxval(indexlist_irred), "Expected:", nwv`

**2. Out of GPU Memory**
```
cudaMalloc failed or kernel launch fails
```
**Solution**: Reduce q-mesh size or increase GPU heap:
```fortran
heapsize = 512_cuda_count_kind*1024_cuda_count_kind*1024_cuda_count_kind
ierr = cudaDeviceSetLimit(cudaLimitMallocHeapSize, heapsize)
```
---

## Validation

Compare GPU results against CPU reference implementation:

```fortran
! Compute relative error
max_error = maxval(abs(V2_gpu - V2_cpu) / (abs(V2_cpu) + 1e-12))
print*, "Max relative error:", max_error

! Typical acceptable error: < 1e-10 (numerical precision)
```
---

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Commit your changes (`git commit -m 'describe your change'`)
4. Push to the branch (`git push origin`)
5. Open a Pull Request

### Development Guidelines

- Follow Fortran 2008+ standards
- Add tests for new features
- Document GPU memory requirements
- Profile performance with `nvprof` or `nsys`

---

## Citation

If you use V3gpu in your research, please cite:

```bibtex
@software{v3gpu,
  title = {V3gpu: GPU-Accelerated Three-Phonon Interactions},
  author = {Sally Issa},
  year = {2025},
  url = {https://github.com/sallyi95/V3gpu}
}
```

---

## License

This project is licensed under the ...License - see the [not selected yet](LICENSE) file for details.

---

## Related Projects

- [elphbolt](https://github.com/nakib/elphbolt) - electron-phonon Boltzmann transport solver
- [CUDA Fortran](https://developer.nvidia.com/cuda-fortran) - NVIDIA's Fortran compiler

---

## Contact

For questions, bug reports, or feature requests:
- **Issues**: [GitHub Issues](https://github.com/sallyi95/V3gpu/issues)
- **Email**: sally.issa@physik.hu-berlin.de

---

## Acknowledgments

- Elphbolt developers for algorithmic inspiration

---
