include("parameters.jl")

function Bose(ħω, T)
    #Returns the Bose-Einstein distribution for energy ħω
    #in eV and temperature T in K.

    return 1.0/expm1(ħω/kB/T)
end

function Fermi(ε, μ, T)
    #Returns the Fermi-Dirac distribution for energy ε
    #in eV, chemical potential μ in eV, and temperature
    #T in K.

    return 1.0/(exp((ε - μ)/kB/T) + 1.0)
end
