# NORMAL DERIVATIVE ON BOUNDARY

using LinearAlgebra
using IterativeSolvers: gmres
using SpecialFunctions
using Plots

includet("Talbot.jl")
includet("Yukawa.jl")

using .Yukawa
using .Talbot

# N = 64; θ = [2π/N * k for k in 0:N-1];
# dΩ = exp.(im*θ);
# s = complex(1.0)
# x₀ = [0 - 8im, 0 + 8im]
# x₀ = 5.0 + 0im

f(x, s) = sum([1/2π * besselk(0, sqrt(s) * norm(x - x₀)) for x₀ in x₀])

function normal_derivative(dΩ, x₀, f; s)
    n̂, J, κ = Yukawa.GeomDerivs(dΩ);
    σ, gmres_log = gmres(-I/2 + Yukawa.DLP(dΩ, s), -f.(dΩ, s); reltol = 1e-12, maxiter=128, log=true);
    
    x = [real(dΩ) imag(dΩ)]
    K1(i, j) = besselk(1, sqrt(s) * norm(x[i,:] - x[j,:]))
    K2(i, j) = norm(x[i,:] - x[j,:])
    K3(i, j) = dot(x[i,:] - x[j,:], [real(n̂[j]) imag(n̂[j])])
    dK1(i, j) = sqrt(s) * (-besselk(0, sqrt(s) * norm(x[i,:] - x[j,:])) - besselk(1, sqrt(s) * norm(x[i,:] - x[j,:])) / norm(x[i,:] - x[j,:]) / sqrt(s)) * (x[i,:] - x[j,:])' / norm(x[i,:] - x[j,:])
    dK2(i, j) = (x[i,:] - x[j,:])' / norm(x[i,:] - x[j,:])
    dK3(i, j) = [real(n̂[j]) imag(n̂[j])]
    K(i, j) = sqrt(s) / 2π * dot([real(n̂[i]) imag(n̂[i])], (dK1(i, j) * K2(i, j) - K1(i, j) * dK2(i, j)) * K3(i, j) / K2(i, j) / K2(i, j) + K1(i, j) * dK3(i, j) / K2(i, j))

    # (K - λ) term
    dxdylog(i, j) = -real(dot(n̂[i], n̂[j])) / abs(dΩ[i] - dΩ[j])^2 + 2(real(dot(dΩ[i] - dΩ[j], n̂[i])) * real(dot(dΩ[i] - dΩ[j], n̂[j]))) / abs(dΩ[i] - dΩ[j])^4
    λ = -1 / 2π
    oddeven1(i) = sum([(K(i,j) - λ * dxdylog(i, j)) * σ[j] * J[j] * 2 * 2π / length(dΩ) for j in (iseven(i) ? (1:2:length(dΩ)) : (2:2:length(dΩ)))])

    # + λ term
    oddeven2(i) = sum([λ * dxdylog(i, j) * (σ[j] - σ[i]) * J[j] * 2 * 2π / length(dΩ) for j in (iseven(i) ? (1:2:length(dΩ)) : (2:2:length(dΩ)))])

    dĉdnₓ = [oddeven1(i) + oddeven2(i) for i in 1:length(dΩ)]

    # add normal derivative of particular term
    return dĉdnₓ .- sum([-sqrt(s)/2π * besselk.(1, sqrt(s) .* abs.(dΩ .- x₀)) .* real(dot.(dΩ .- x₀, n̂)) ./ abs.(dΩ .- x₀) for x₀ in x₀])

end

# invert laplace transform
Nᵧ = 32; γ = TalbotContour(
    # z(θ; t)
    (θ; t=nothing) -> length(θ) / 2t * (-0.6122 .+ 0.5017 * θ .* cot.(0.6407 * θ) + 0.6245 * im .* θ),
    # z′(θ; t)
    (θ; t=nothing) -> length(θ) / 2t * (-0.32143919 * θ .* csc.(0.6407*θ).^2 + 0.5017 * cot.(0.6407*θ) .+ 0.6245*im)
)
dcdnₓ = talbot_midpoint(s -> normal_derivative(dΩ, x₀, f; s), γ, Nᵧ)






#################
### UNIT BALL ###
#################

x₀ = 5.0 + 0im; N = 64; θ = [2π/N * k for k in 0:N-1]
dΩ = exp.(im*θ)

### PLOT TOTAL FLUX 
T = range(0.1; step=0.5, stop=7.0)
dcdnₓ = talbot_midpoint(s -> normal_derivative(dΩ, x₀, f; s), γ, Nᵧ)
plt = plot(T, [abs.(sum(dcdnₓ(t))) * 2π / N for t in T], label=:none)
xaxis!("Time"); yaxis!("Flux")
savefig(plt, "flux.png")

### PLOT POINTWISE FLUX
T = [0.1, 0.5, 1.0, 5.0, 10.0]
plt = plot()
for t in T
    plot!(θ, abs.(dcdnₓ(t)), label="t=$t")
end
display(plt)
savefig(plt, "pointwise.png")

### CONVERGENCE
# exact total flux from moth paper
dw=0.0001; flux_exact(t) = 2/π * sum([(besselj0(w) * bessely0(w * norm(x₀)) - besselj0(w * norm(x₀)) * bessely0(w)) / (bessely0(w)^2 + besselj0(w)^2) * w * exp(-w^2t) * dw for w in range(dw; stop=111, step=dw)])
# plot convergence
t_test = 1.0; err = []
for N in [64, 128, 256, 512, 1024]
    θ = [2π/N * k for k in 0:N-1]; dΩ = exp.(im*θ)
    dcdnₓ = talbot_midpoint(s -> normal_derivative(dΩ, x₀, f; s), γ, Nᵧ)
    fl_ex = flux_exact(t_test); fl_ap = abs(sum(dcdnₓ(t_test))) * 2π / N
    push!(err, abs(fl_ex - fl_ap) / abs(fl_ex))
end
plot([64, 128, 256, 512, 1024], err, xaxis=:log, yaxis=:log, marker=:square, label=:none, title="Error in Flux, Unit Ball")
xlabel!("# Boundary Points")
ylabel!("Error")
savefig("convergence.png")

### MAX TIME WRT X0
runs = [];
for x in range(2.5, stop=10.0, step=0.5)
    x₀ = complex(x)
    T = 10 .^ range( log10(0.1), log10(10), length = 30 )
    dcdnₓ = talbot_midpoint(s -> normal_derivative(dΩ, x₀, f; s), γ, Nᵧ)
    @time push!(runs, [abs.(sum(dcdnₓ(t))) * 2π / N for t in T])
end

plot(); xaxis!("Time"); yaxis!("Flux")
for run in runs
    plot!(T, run, label=:none, yaxis=:log, xaxis=:log)
end
savefig("multiple_flux.png")

scatter([T[argmax(x)] for x in runs], [x[argmax(x)] for x in runs], label=:none, xaxis=:log, yaxis=:log)
xaxis!("Time"); yaxis!("Max Flux")
linreg(x, y) = hcat(fill!(similar(x), 1), x) \ y
linreg(log10.([T[argmax(x)] for x in runs][1:12]), log10.([x[argmax(x)] for x in runs][1:12]))
savefig("maxtimes.png")

times = []
for run in runs
    
    # maxtime = T[argmax(run)]
    # push!(times, maxtime)
end
plot(range(2.5, stop=10.0, step=0.5), times, xaxis=:log, yaxis=:log, label=:none)
xaxis!("Distance Source to Origin"); yaxis!("Time at Max Flux")

savefig("maxtimes.png")


################
### STARFISH ###
################

x₀ = 5.0 + 0im; N = 128; θ = [2π/N * k for k in 0:N-1]
dΩ = exp.(im*θ) .* (3 .+ sin.(5θ))

### PLOT TOTAL FLUX 
T = range(0.1; step=0.05, stop=1.0)
dcdnₓ = talbot_midpoint(s -> normal_derivative(dΩ, x₀, f; s), γ, Nᵧ)
plt = plot(T, [abs.(sum(dcdnₓ(t))) * 2π / N for t in T], label=:none)
xaxis!("Time"); yaxis!("Flux")
savefig(plt, "flux.png")

### PLOT POINTWISE FLUX
T = [0.1, 0.5, 1.0, 5.0, 10.0]
plt = plot()
for t in T
    plot!(θ, abs.(dcdnₓ(t)), label="t=$t")
end
display(plt)
savefig(plt, "pointwise.png")

### PLOT CONVERGENCE
t_test = 1.0; N = 2048
θ = [2π/N * k for k in 0:N-1]; dΩ = exp.(im*θ) .* (3 .+ sin.(5θ))
dcdnₓ = talbot_midpoint(s -> normal_derivative(dΩ, x₀, f; s), γ, Nᵧ)
fl_ex = abs(sum(dcdnₓ(t_test))) * 2π / N
err = []
for N in [64, 128, 256, 512, 1024]
    θ = [2π/N * k for k in 0:N-1];
    dΩ = exp.(im*θ) .* (3 .+ sin.(5θ))
    dcdnₓ = talbot_midpoint(s -> normal_derivative(dΩ, x₀, f; s), γ, Nᵧ)
    fl_ap = abs(sum(dcdnₓ(t_test))) * 2π / N
    push!(err, abs(fl_ex - fl_ap) / abs(fl_ex))
end
plt = plot([64, 128, 256, 512, 1024], err, xaxis=:log, yaxis=:log, marker=:square, label=:none, title="Error in Flux, Starfish")
xlabel!("# Boundary Points")
ylabel!("Error")
savefig(plt, "convergence.png")


################
### CRESCENT ###
################

x₀ = complex(0.0); N = 128; θ = [2π/N * k for k in 0:N-1]
dΩ = [π/2 <= θ < 3π/2 ? exp(im*θ) : -0.5*cos(θ) + im*sin(θ) for θ in θ]

### PLOT TOTAL FLUX 
T = range(0.1; step=0.05, stop=1.0)
dcdnₓ = talbot_midpoint(s -> normal_derivative(dΩ, x₀, f; s), γ, Nᵧ)
plt = plot(T, [abs.(sum(dcdnₓ(t))) * 2π / N for t in T], label=:none)
xaxis!("Time"); yaxis!("Flux")
savefig(plt, "flux.png")

### PLOT POINTWISE FLUX
T = [0.1, 0.5, 1.0, 5.0, 10.0]
plt = plot()
for t in T
    plot!(θ, abs.(dcdnₓ(t)), label="t=$t")
end
display(plt)
savefig(plt, "pointwise.png")

### PLOT CONVERGENCE
t_test = 1.0; N = 2048
θ = [2π/N * k for k in 0:N-1]; dΩ = exp.(im*θ) .* (3 .+ sin.(5θ))
dcdnₓ = talbot_midpoint(s -> normal_derivative(dΩ, x₀, f; s), γ, Nᵧ)
fl_ex = abs(sum(dcdnₓ(t_test))) * 2π / N
err = []
for N in [64, 128, 256, 512, 1024]
    θ = [2π/N * k for k in 0:N-1];
    dΩ = [π/2 <= θ < 3π/2 ? exp(im*θ) : -0.5*cos(θ) + im*sin(θ) for θ in θ]
    dcdnₓ = talbot_midpoint(s -> normal_derivative(dΩ, x₀, f; s), γ, Nᵧ)
    fl_ap = abs(sum(dcdnₓ(t_test))) * 2π / N
    push!(err, abs(fl_ex - fl_ap) / abs(fl_ex))
end
plt = plot([64, 128, 256, 512, 1024], err, xaxis=:log, yaxis=:log, marker=:square, label=:none, title="Error in Flux, Crescent")
xlabel!("# Boundary Points")
ylabel!("Error")
savefig(plt, "convergence.png")




############
### BLOB ###
############

x₀ = 1.0 + im * 1.0; N = 128; θ = [2π/N * k for k in 0:N-1]
dΩ = (sin.(0.5θ).^3 .+ cos.(0.5θ).^3) .* exp.(im * 0.5θ)

### PLOT TOTAL FLUX 
T = 10 .^ range( log10(0.01), log10(0.2), length = 10 )
dcdnₓ = talbot_midpoint(s -> normal_derivative(dΩ, x₀, f; s), γ, Nᵧ)
plt = plot(T, [abs.(sum(dcdnₓ(t))) * 2π / N for t in T], label=:none)
xaxis!("Time"); yaxis!("Flux")
savefig(plt, "flux.png")

### PLOT POINTWISE FLUX
T = [0.1, 0.5, 1.0, 5.0, 10.0]
plt = plot()
for t in T
    plot!(θ, abs.(dcdnₓ(t)), label="t=$t")
end
display(plt)
savefig(plt, "pointwise.png")

### PLOT CONVERGENCE
t_test = 1.0; N = 2048
θ = [2π/N * k for k in 0:N-1]; dΩ = exp.(im*θ) .* (3 .+ sin.(5θ))
dcdnₓ = talbot_midpoint(s -> normal_derivative(dΩ, x₀, f; s), γ, Nᵧ)
fl_ex = abs(sum(dcdnₓ(t_test))) * 2π / N
err = []
for N in [64, 128, 256, 512, 1024]
    @show N
    θ = [2π/N * k for k in 0:N-1];
    dΩ = exp.(im*θ) .* (3 .+ sin.(5θ))
    dcdnₓ = talbot_midpoint(s -> normal_derivative(dΩ, x₀, f; s), γ, Nᵧ)
    fl_ap = abs(sum(dcdnₓ(t_test))) * 2π / N
    push!(err, abs(fl_ex - fl_ap) / abs(fl_ex))
end
plt = plot([64, 128, 256, 512, 1024], err, xaxis=:log, yaxis=:log, marker=:square, label=:none, title="Error in Flux, Bean")
xlabel!("# Boundary Points")
ylabel!("Error")
savefig(plt, "convergence.png")