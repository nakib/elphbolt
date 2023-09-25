using DelimitedFiles
using LinearAlgebra
using Plots
using LaTeXStrings

qe = 1.602176634e-19 #C
hbar = 1.05457172647e-22 #J.ps
hbar_eVps = hbar/qe #ev.ps

function δ_gaussian(ω, ε, σ)
    #Returns the gaussian representation of δ(ω - ε)
    #for a given smearing σ.

    return exp(-((ω - ε)/σ)^2/2)/σ/√(2π)
end

function dos(ωs, εs, vs, Bs, wvmesh)
    #Returns the density of states (DOS) calculated
    #on a uniform mesh.
    #
    #The adaptive gaussian smearing method is used here.

    N = length(ωs)

    #Number of wave vectors and bands
    nq, nb = size(εs)

    #Hard-code a minimum smearing value.
    σ_min = 1.0e-6 # 1 μeV
    
    Qs = zeros(3, 3)
    for dim ∈ 1:3
        Qs[dim, :] = Bs[dim, :]/wvmesh[dim]
    end
    
    dos = zeros(N)
    for i ∈ 1:N
        for iqp ∈ 2:nq #skip the Γ-point
            for ibp ∈ 1:nb
                σ = 0.0
                v_ibp = vs[iqp, ibp, :]
                for dim ∈ 1:3
                    σ += (v_ibp[:]⋅Qs[dim, :])^2
                end                    
                σ = max(hbar_eVps*√(σ/12), 1.0e-4)
                
                dos[i] += δ_gaussian(ωs[i], εs[iqp, ibp], σ)
            end
        end
    end
    return dos/nq
end

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
N = 200

#Sampling energies
ωs = collect(LinRange(0, 1.05*maximum(εs), N)) #eV

#Phonon DOS
ph_dos = dos(ωs, εs, vs, Bs, wvmesh)

#Plot in meV-meV^-1 units
dos_plot = plot(1.0e3*ωs, 1.0e-3*ph_dos, legend = false)
plot!(xlabel = "phonon energy [meV]")
plot!(ylabel = "DOS [meV]"*L"$^{-1}$")
savefig(dos_plot, path*"phdos.pdf")

#display(dos_plot)
