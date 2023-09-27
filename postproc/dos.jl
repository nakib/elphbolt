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
    return dos/prod(wvmesh)
end

function calculate_dos(species, rundir, outdir, chempot)
    #Calculate the density of states of a given 
    #species ∈ ["ph", "el"].

    #Read full-Brillouin zone energies and velocities. Reshape the latter appropriately.
    println("ϟ Reading full-Brillouin zone energies and velocities...")
    εs = readdlm(rundir*species*".ens_fbz") #eV
    vs = readdlm(rundir*species*".vels_fbz") #Km/s
    vs_shape = size(vs)
    vs = reshape(vs, (vs_shape[1], vs_shape[2]÷3, 3))

    #Read irreducible Brillouin zone energies and flatten.
    println("ϟ Reading irreducible Brillouin zone energies...")
    εs_irred = readdlm(rundir*species*".ens_ibz") #eV
    εs_irred_shape = size(εs_irred)
    ωs = reshape(εs_irred, (prod(εs_irred_shape)))

    #Read reciprocal lattice vectors and their discretization grid
    println("ϟ Reading reciprocal lattice vector data...")
    reclattvecs_data = readdlm(rundir*species*".reclattvecs") #nm^-1
    Bs = reclattvecs_data[1:3, :]
    wvmesh = convert.(Int32, reclattvecs_data[4, :])
    
    #Calculate band resolved DOS
    println("ϟ Calculating density of states...")
    species_dos = reshape(ThreadsX.map(ω -> dos(ω, εs, vs, Bs, wvmesh), ωs), εs_irred_shape)
    if(species == "ph")
        species_dos[1, 1:3] .= 0.0 #Take care of the 3 Γ-point acoustic phonons
    end

    #Write DOS to file
    println("ϟ Writing results out to file...")
    open(outdir*species*".dos_branches", "w") do w
        writedlm(w, species_dos)
    end
    
    #Plot in meV^-1-mev units
    println("ϟ Plotting results...")
    if(species == "ph")
        dos_plot = plot(1.0e3*εs_irred, 1.0e-3*species_dos,
                        legend = false, seriestype = :scatter, mc = :green, ma = 0.6,
                        xlabel = "phonon energy [meV]",
                        ylabel = "DOS [meV]"*L"$^{-1}$")
    else
        #Set hald of Fermi shell width
        half_Fermi_shell = 0.2 #eV

        #Find (closest) Cartesian index to the chemical potential
        chempot_index = argmin(abs.(εs_irred .- chempot))

        #Set vertical range
        padding = 1.25 #extra vertical space
        ymin = (2.0 - padding)*species_dos[chempot_index]
        ymax = padding*species_dos[chempot_index]
        
        dos_plot = plot(εs_irred .- chempot, species_dos,
                        xlims = [-half_Fermi_shell, half_Fermi_shell],
                        ylims = [ymin, ymax],
                        legend = false, seriestype = :scatter, mc = :green, ma = 0.6,
                        xlabel = "electron energy - "*L"\mu"*" [eV]",
                        ylabel = "DOS [eV.spin]"*L"$^{-1}$")
    end
    savefig(dos_plot, outdir*species*"dos_bands.pdf")

    #display(dos_plot)
end

function parse_commandline()
    args = ArgParseSettings()
    @add_arg_table args begin
        "--rundir"
        help = "elphbolt run directory"
        arg_type = String
        default = "./"

        "--dosof"
        help = "(ph)onon or (el)ectron density of states"
        arg_type = String
        default = "ph"

        "--chempot"
        help = "Chemical potential"
        arg_type = Float64
        default = 0.0
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

    #Set chemical potential
    chempot = parsed_args["chempot"]

    calculate_dos(parsed_args["dosof"], rundir, outdir, chempot)

    println("ϟ All done!")
end

main()
