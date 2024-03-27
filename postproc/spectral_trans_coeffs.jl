include("parameters.jl")
include("statistics.jl")

using ArgParse
using DelimitedFiles

function calculate_mfp_cumulative_kappa(rundir, outdir, T)
    # TODO ...

    #TODO read info from input.nml and calculate volume
    volume = 0.39409804E-01 #nm^3, primitive cell volume
    
    #TODO generalize this to allow "el" also
    species = "ph"
    
    #TODO Generate T-dependent directory from T
    Tdir = "T0.300E+03/"
    
    #Read full-Brillouin zone energies and velocities. Reshape the latter appropriately.
    println("ϟ Reading full-Brillouin zone energies...")
    εs = readdlm(rundir*species*".ens_fbz") #eV
    εs_shape = size(εs)
    nq_fbz, nb = εs_shape[1], εs_shape[2] #number of FBZ points, bands
    
    #Read full-Brillouin zone energies and velocities. Reshape the latter appropriately.
    println("ϟ Reading full-Brillouin zone mean-free-paths (mfps)...")
    λs = readdlm(rundir*Tdir*"nodrag_iterated_ph_mfps_ibz") #nm

    #Create mfp sampling grid (log scale)
    low = -6 #1e-6 nm, a really tiny number I'd say.
    high = log10(1.5*maximum(λs)) #50% higher than largest value in λs

    N = 25 # TODO should read this from user input
    #λ' = exp10.(range(low, high, length = N))
    λp = exp10.(range(low, high, length = N))
    
    #Common multiplicative factor
    fac = 1.0e21/kB/T/volume/nq_fbz
end

function parse_commandline()
    args = ArgParseSettings()
    @add_arg_table args begin
        "--rundir"
        help = "elphbolt run directory"
        arg_type = String
        default = "./"

        "--particle"
        help = "(ph)onon or (el)ectron spectral coefficients"
        arg_type = String
        default = "ph"

        "--chempot"
        help = "Chemical potential"
        arg_type = Float64
        default = 0.0

        "--T"
        help = "Temperature"
        arg_type = Float64
        default = 0.0
    end
    return parse_args(args)
end

function main()
    println("ϟ Welcome to the elphbolt post-processor tool spectral_trans_coeffs.jl")
    
    println("ϟ I'll be using $(Threads.nthreads()) threads.")

    #Parse command line arguments
    parsed_args = parse_commandline()

    #Set elphbolt run directory
    rundir = parsed_args["rundir"]

    #Set postproc output directory
    outdir = rundir*"postproc_results/"
    mkpath(outdir)

    #Set chemical potential
    chempot = parsed_args["chempot"]

    #Set chemical potential
    chempot = parsed_args["T"]

    #TODO...
    T = 300
    calculate_mfp_cumulative_kappa(rundir, outdir, T)

    println("ϟ All done!")
end

main()
