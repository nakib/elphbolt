include("parameters.jl")
include("statistics.jl")

using ThreadsX
using ArgParse
using DelimitedFiles

function Heaviside(x::Float64)
    (1.0 + sign(x))/2.0
end

function calculate_mfp_cumulative_kappa(rundir, outdir, T)
    # TODO ...

    #TODO read info from input.nml and calculate volume
    volume = 0.39409804E-01 #nm^3, primitive cell volume
    
    #TODO generalize this to allow "el" also
    species = "ph"
    
    #TODO Generate T-dependent directory from T
    Tdir = "T0.300E+03/"

    #Read map from full Brillouin zone to the irreducible irreducible one.
    println("ϟ Reading FBZ->IBZ map...")
    fbz2ibz = readdlm(rundir*species*".fbz2ibz_map", Int32)
    
    #Read full Brillouin zone energies and velocities. Reshape the latter appropriately.
    println("ϟ Reading FBZ energies and velocities...")
    εs = readdlm(rundir*species*".ens_fbz") #eV
    vs = readdlm(rundir*species*".vels_fbz") #Km/s
    vs_shape = size(vs)
    vs = reshape(vs, (vs_shape[1], vs_shape[2]÷3, 3))

    εs_shape = size(εs)
    nq_fbz, nb = εs_shape[1], εs_shape[2] #number of FBZ points, bands
    
    #Read full Brillouin zone energies.
    println("ϟ Reading full-Brillouin zone mean-free-paths (mfps)...")
    λs = readdlm(rundir*Tdir*"nodrag_iterated_ph_mfps_ibz") #nm

    #Read full Brillouin zone reponse functions.

    response = zeros(nb, nq_fbz, 3)
    for ib in 1:nb
        fname = "nodrag_F0_"*string(ib)
        response[ib, :, :] = readdlm(rundir*Tdir*fname)
    end
    
    #Create mfp sampling grid (log scale)
    #low = -6 
    #high = log10(1.5*maximum(λs)) 
    nsamp = 25 # TODO should read this from user input
    λsamp = exp10.(range(-6, #1e-6 nm, a really tiny number I'd say.
                         log10(1.5*maximum(λs)) , #50% higher than largest value in λs
                         nsamp))

    #Calculate band resolved DOS
    println("ϟ Calculating phonon thermal conductivity accumulation wrt mean-free-path...")
    
    #kappa_accum = reshape(ThreadsX.map(l -> mfp_cumulative_kappa(l, εs), λsamp), [nsamp, nb])
    
    kappa_accum = zeros(nsamp, nb, 3, 3)
    for isamp in 1:nsamp
        for ib in 1:nb
            for iq_fbz in 1:nq_fbz
                iq_ibz = fbz2ibz[iq_fbz]
                if(λs[iq_ibz, ib] <= λsamp[isamp])    
                    ε = εs[iq_fbz, ib]

                    population_fac = Bose(ε, T)
                    population_fac *= 1.0 + population_fac

                    for icart in 1:3
                        kappa_accum[isamp, ib, icart, :] +=
                            ε*population_fac*vs[iq_fbz, ib, icart]*response[ib, iq_fbz, :]
                    end
                end
            end
        end
    end
    #Common multiplicative factor
    fac = qe*1.0e21/kB/T/volume/nq_fbz
    kappa_accum *= fac
    
    #Write to file
    #println("ϟ Writing results out to file...")
    #open(outdir*species*".kappa_accum_wrt_mfp", "w") do w
    #    writedlm(w, kappa_accum)
    #end
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
