using DelimitedFiles
using LinearAlgebra
using Plots
using LaTeXStrings
using ThreadsX
#using Distributed

qe = 1.602176634e-19 #C
hbar = 1.05457172647e-22 #J.ps
hbar_eVps = hbar/qe #ev.ps

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
            v_ib′ = vs[iq′, ib′, :]
            for dim ∈ 1:3
                σ += (v_ib′[:]⋅Qs[dim, :])^2
                #σ += dot(v_ib′[:], Qs[dim, :])^2
            end
            σ = max(hbar_eVps*√(σ/12), σ_min)
            
            dos += δ_gaussian(ω, εs[iq′, ib′], σ)
        end
    end
    return dos/nq
end

println(Threads.nthreads())

#Test on MgB2
#
#Path to elphbolt/superconda calculation
path = "/users/sol/nakib/elphbolt/examples/MgB2/superconda/gaussian_25c_0.35/"

#Read full-Brillouin zone energies and velocities. Reshape the latter appropriately.
εs = readdlm(path*"ph.ens_fbz") #eV
vs = readdlm(path*"ph.vels_fbz") #Km/s
vs_shape = size(vs)
vs = reshape(vs, (vs_shape[1], vs_shape[2]÷3, 3))

#Reciprocal lattice vectors 1/nm
Bs = [0.20380103E+02   0.11766458E+02  -0.00000000E+00;
      0.00000000E+00   0.23532916E+02   0.00000000E+00;
      0.00000000E+00  -0.00000000E+00   0.17844890E+02]

#Wave vector grid
wvmesh = [25 25 25]

#Number of equdistant energy sampling points
N = 400

#Sampling energies
ωs = collect(LinRange(0, 1.05*maximum(εs), N)) #eV

#Phonon DOS
ph_dos = ThreadsX.map(ω -> dos(ω, εs, vs, Bs, wvmesh), ωs)
ph_dos[1] = 0.0 #Take care of the Γ-point

#Plot in meV-meV^-1 units
dos_plot = plot(1.0e3*ωs, 1.0e-3*ph_dos, legend = false)
plot!(xlabel = "phonon energy [meV]")
plot!(ylabel = "DOS [meV]"*L"$^{-1}$")
savefig(dos_plot, path*"phdos.pdf")

#display(dos_plot)
