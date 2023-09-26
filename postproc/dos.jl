include("parameters.jl")

using DelimitedFiles
using LinearAlgebra
using Plots
using LaTeXStrings
using ThreadsX
using ArgParse

function δ_gaussian(ω, ε, σ)
    #Returns the gaussian representation of δ(ω - ε)
    #for a given smearing σ.

    return exp(-((ω - ε)/σ)^2/2)/σ/√(2π)
end

function dos(ω, εs, vs, Bs, wvmesh)
    #Returns the density of states (DOS) calculated
    #on a uniform mesh at a given energy.
    #
    #The adaptive gaussian smearing method is used here.

    #Number of wave vectors and bands
    nq, nb = size(εs)

    #Hard-code a minimum smearing value.
    σ_min = 1.0e-6 # 1 μeV
    
    Qs = zeros(Float64, 3, 3)
    for dim ∈ 1:3
        Qs[dim, :] = Bs[dim, :]/wvmesh[dim]
    end
    
    dos = 0.0
    for iq′ ∈ 1:nq
        for ib′ ∈ 1:nb
            #Calculate the adaptive smearing parameter
            σ = 0.0
            for dim ∈ 1:3
                σ += (vs[iq′, ib′, :]⋅Qs[dim, :])^2
            end
            σ = max(hbar_eVps*√(σ/12), σ_min)
            
            dos += δ_gaussian(ω, εs[iq′, ib′], σ)
        end
    end
    return dos/nq
end

function parse_commandline()
    args = ArgParseSettings()
    @add_arg_table args begin
        "--rundir"
            help = "elphbolt run directory"
            arg_type = String
            default = "./"
    end
    return parse_args(args)
end

function main()
    println("ϟ Welcome to the elphbolt post-processor!")
    
    println("ϟ I'll be using $(Threads.nthreads()) threads.")

    #Parse command line arguments
    parsed_args = parse_commandline()

    #Set elphbolt run directory
    rundir = parsed_args["rundir"]

    #Set postproc output directory
    outdir = rundir*"postproc_results/"
    mkpath(outdir)

    #Read full-Brillouin zone energies and velocities. Reshape the latter appropriately.
    println("ϟ Reading full-Brillouin zone energies and velocities...")
    εs = readdlm(rundir*"ph.ens_fbz") #eV
    vs = readdlm(rundir*"ph.vels_fbz") #Km/s
    vs_shape = size(vs)
    vs = reshape(vs, (vs_shape[1], vs_shape[2]÷3, 3))

    #Read irreducible Brillouin zone energies and flatten.
    println("ϟ Reading irreducible Brillouin zone energies...")
    εs_irred = readdlm(rundir*"ph.ens_ibz") #eV
    εs_irred_shape = size(εs_irred)
    ωs = reshape(εs_irred, (prod(size(εs_irred))))

    #Read reciprocal lattice vectors and their discretization grid
    println("ϟ Reading reciprocal lattice vector data...")
    reclattvecs_data = readdlm(rundir*"ph.reclattvecs") #nm^-1
    Bs = reclattvecs_data[1:3, :]
    wvmesh = reclattvecs_data[4, :]
        
    #Calculate band resolved phonon DOS
    println("ϟ Calculating phonon density of states...")
    ph_dos = reshape(ThreadsX.map(ω -> dos(ω, εs, vs, Bs, wvmesh), ωs), εs_irred_shape)
    ph_dos[1, 1:3] .= 0.0 #Take care of the 3 Γ-point acoustic phonons

    #Write DOS to file
    println("ϟ Writing results out to file...")
    open(outdir*"ph.dos_branches", "w") do w
        writedlm(w, ph_dos)
    end
    
    #Plot in meV^-1-mev units
    println("ϟ Plotting results...")
    dos_plot = plot(1.0e3*εs_irred, 1.0e-3*ph_dos,
                    legend = false, seriestype = :scatter, mc = :green, ma = 0.6,
                    xlabel = "phonon energy [meV]",
                    ylabel = "DOS [meV]"*L"$^{-1}$")
    savefig(dos_plot, outdir*"phdos_bands.pdf")

    #display(dos_plot)

    println("ϟ All done!")
end

main()
