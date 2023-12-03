include("statistics.jl")

using ArgParse

function calculate_spectral_coeffs(species, rundir, outdir, T, chempot)
    # TODO...

    #TODO Set temperature-dependent directory
    #Tdir = rundir*...
    
    #Read full-Brillouin zone energies and velocities. Reshape the latter appropriately.
    println("ϟ Reading full-Brillouin zone energies and velocities...")
    εs = readdlm(rundir*species*".ens_fbz") #eV
    vs = readdlm(rundir*species*".vels_fbz") #Km/s
    vs_shape = size(vs)
    vs = reshape(vs, (vs_shape[1], vs_shape[2]÷3, 3))

    #Read reciprocal lattice vectors and their discretization grid
    println("ϟ Reading reciprocal lattice vector data...")
    reclattvecs_data = readdlm(rundir*species*".reclattvecs") #nm^-1
    Bs = reclattvecs_data[1:3, :]
    wvmesh = convert.(Int32, reclattvecs_data[4, :])

    #TODO Read response function
    # RTA
    println("ϟ Reading reciprocal lattice vector data...")
    response_fn = readdlm(Tdir*"RTA_F0_tot") #nm.eV/K
    # TODO Full
    #...

    #TODO (Can postpone this for later)
    # Talk to spglib to generate symmetry operations
    # to symmetrize the transport coefficients.
    # This might be an overkill since elphbolt already
    # symmetrizes the response functions and the velocities.
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

    println("ϟ All done!")
end

main()
