### TODO ###
## compare to fine res integral from alan's paper
## change to neumann and check concentration of mass
## changing point on contour need diff hankel1?  show that hankel2 solves the ode.  show that hankelh1 does not.  show that this holds regardless of pos/neg root
## integral equation yukawa solver

using Plots
using SpecialFunctions
using FFTW
using LinearAlgebra
using IterativeSolvers
using MATLAB

includet("Talbot.jl")
includet("Fourier.jl")
includet("Yukawa.jl")

using .Talbot
using .Fourier
using .Yukawa

### BOUNDARY DISCRETIZATION ###
N = 128; θ = θ = reverse!([2π/(N-1) * k for k in 1:N])
dΩ = exp.(im*θ) .* (3 .+ sin.(5θ))
# dΩ = cos.(θ) + im * 3sin.(θ)

### EVALUATION GRID ###
Nᵣ = 32; Nᵩ = 128; Rmin = 1.1; Rmax = 10;
r = range(Rmin; stop=Rmax, length=Nᵣ)
φ = [2π/(Nᵩ-1) * k for k in 1:Nᵩ]
Ω = exp.(im*φ) .* (3 .+ sin.(5φ)) * r'
# Ω = r * (cos.(θ) + im * 3sin.(θ))'
# plot(dΩ, aspect_ratio=1)
# scatter!(Ω, color=:black, label=:none, markersize=1)

### BOUNDARY CONDITION ###
x₀ = [10+5im, -10+5im, -8im]
f(x, s) = sum([1/2π * besselk(0, sqrt(s) * norm(x - x₀)) for x₀ in x₀])

### TALBOT CONTOUR ###
Nᵧ = 32; γ = TalbotContour(
    # z(θ; t)
    (θ; t=nothing) -> length(θ) / 2t * (-0.6122 .+ 0.5017 * θ .* cot.(0.6407 * θ) + 0.6245 * im .* θ),
    # z′(θ; t)
    (θ; t=nothing) -> length(θ) / 2t * (-0.32143919 * θ .* csc.(0.6407*θ).^2 + 0.5017 * cot.(0.6407*θ) .+ 0.6245*im)
)

function BIESolve(t)
    
    @info "Time $t"

    ### GET CONTOUR POINTS ###
    s = Talbot.γ(γ, Nᵧ, t); s′ = Talbot.γ′(γ, Nᵧ, t);
    
    ### SOLVE SKIE AT TALBOT CONTOUR POINTS ###
    @time C = [YukawaSolve(dΩ, f, s) for s in s]
    
    ### BROMWICH INTEGRAL ###
    c(x,t) = -im / Nᵧ * sum(exp.(s*t) .* [C(x,s) for (C,s) in zip(C,s)] .* s′)

    ### EVALUATE ON TARGET GRID Ω ###
    return @time [c(x, t) for x in Ω]
    
end

### LOOP OVER TIME ###
# t = [0.1, 0.5, 2, 10, 20, 40, 60, 80, 100]
t = [0.5, 2, 20, 100]
c = [BIESolve(t) for t in t]

### PLOTTING ###
for i in 1:length(t)
    z = real(c[i]);
    mat"""
    i = $i;
    figure()
    surf(real($Ω), imag($Ω), $z)
    view(2)
    shading interp
    axis equal
    axis([-15 15 -15 15])
    colorbar
    exportgraphics(gcf,"star" + i + ".png",'Resolution',300)
    """
end
