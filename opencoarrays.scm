(use-modules
 (guix packages)
 (guix git-download)
 (guix build-system cmake)
 (guix licenses)
 (gnu packages mpi)
 (gnu packages gcc)
 (gnu packages perl))

(define-public opencoarrays
  (package
   (name "opencoarrays")

   (version "2.9.2")

   (home-page "http://www.opencoarrays.org/")

   (source
    (origin
     (method git-fetch)
     (uri (git-reference
	   (url "https://github.com/sourceryinstitute/OpenCoarrays")
	   (commit "15dc8c3b4707fd4c25742c6091978fd2e30d4f52")))
     (sha256
      (base32
       "1zp8jdw1dks289gjfbm7f7l77fsb0acshq9l3c37s10ppka4p5w9"))
     (file-name (git-file-name name version))))

   (build-system cmake-build-system)

   (arguments
    `(#:phases (modify-phases %standard-phases (add-after
						'build 'mpi-setup ,%openmpi-setup))))

   (propagated-inputs
    `(("openmpi" ,openmpi)
      ("gfortran" ,gfortran)
      ("perl" ,perl)))

   (synopsis "OpenCoarrays is an open-source software project that produces an application binary interface (ABI) used by the GNU Compiler Collection (GCC) Fortran front-end to build executable programs that leverage the parallel programming features of Fortran 2018.")

   (description
    "OpenCoarrays supports the GNU Compiler Collection (GCC) Fortran compiler (gfortran) by providing a parallel application binary interface (ABI) that abstracts away the underlying communication library. OpenCoarrays thus enables gfortran to support Fortran's parallel programming features, often called 'Coarray Fortran,' without making direct reference to the back-end communication library: the Message Passing Interface (MPI). This ensures that Fortran programs and Fortran compilers may take advantage of other communication libraries without costly refactoring. Work is underway on the Caffeine project to support alternative communication libraries and alternative compilers by defining a compiler-independent parallel ABI atop the GASNet-EX exascale networking middleware.

OpenCoarrays provides a compiler wrapper (caf), a parallel runtime library (libcaf_mpi), and a program launcher (cafrun). The wrapper and launcher provide a uniform abstraction for compiling and executing Coarray Fortran without direct reference to the underlying MPI layer.")

   (license bsd-3)))

opencoarrays
