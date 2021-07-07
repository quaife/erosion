ℱ = x -> fftshift(fft(x))
ℱ⁻¹ = x -> ifft(ifftshift(x))
modes = N -> N % 2 == 0 ? [0; 1-N÷2:N÷2-1] : collect(-(N-1)÷2:(N-1)÷2)
∂ = fx -> im * modes(length(fx)) .* fx

N = 64; θ = reverse!([2π/N * k for k in 1:N])
dΩ = exp.(im*θ) .* (3 .+ sin.(5θ))
n̂, J, κ = Yukawa.GeomDerivs(dΩ)
plot(dΩ, aspect_ratio=1, label="")
quiver!(real(dΩ[1:N÷8:end]), imag(dΩ[1:N÷8:end]), quiver=(real(n̂[1:N÷8:end]), imag(n̂[1:N÷8:end])))
quiver!(real(dΩ[1:N÷8:end]), imag(dΩ[1:N÷8:end]), quiver=(real(im*n̂[1:N÷8:end]), imag(im*n̂[1:N÷8:end])))


# BARYCENTRIC QUADRATURE
# rewrite the bessel term as
# K₁ = (K₁ + log(z)) + (-log(z))
# first term smooth, second term not
# apply trapezoid rule to smooth component, barycentric cauchy to log component
# DLP eval is the sum of the two intergrals

# SMOOTH COMPONENT
smooth = besselk.(1, sqrt(s) * abs.(x .- dΩ)) + log.(x .- dΩ)
trap = -sqrt(s)/length(dΩ) * sum(smooth .* real(dot.(x .- dΩ, n̂)) ./ abs.(x .- dΩ) .* σ .* J)

# LOGARITHMIC COMPONENT
σ′ = (ℱ⁻¹ ∘ ∂ ∘ ℱ)(σ)
# boundary data on v(x), where v(x) defined as D[σ](x) = Re(v(x))
v₀ = σ .- 1 / (2π * im) * (2π / length(dΩ)) * [sum(replace!(x -> isnan(x) ? σ₀′ : x, (σ₀ .- σ) ./ (x₀ .- dΩ))) for (x₀, σ₀, σ₀′) in zip(dΩ, σ, σ′)]
# cauchy integral
v = 1 / (2π*im) * sum(v₀ ./ (dΩ .- x) .* (im .* n̂ .* J))
a = sum(dΩ) / length(dΩ) # some point in Ωᶜ
one = 1 / (x - a) / (2π*im) * sum(1 ./ (dΩ .- x) ./ (dΩ .- a) .* (im .* n̂ .* J))
logint = v / one

# approximate f(a) by trapezoid rule
function cauchy_unit_circle(f, a, N)
    θ = 2π/N:2π/N:2π
    z = exp.(im*θ)
    return sum(z .* f.(z) ./ (z .- a)) / N
end

function bary_unit_circle(f, a, N)
    θ = 2π/N:2π/N:2π
    z = exp.(im*θ)
    cauchy = 1 / (2π*im) * sum(z .* f.(z) ./ (z .- a)) * 2π / N
    one = 1 / (2π*im) * sum(1 ./ (z .- a) * im .* z) * 2π / N
    return cauchy / one
end

# test function
f(z::Complex) = exp(z)sin(z^2) + 1
# range of discretization points and singularities to test
N = [4, 8, 16, 32, 64, 128, 256, 512]
a = [0, 0.1im, 0.5im, 0.9im, 0.99im]
# matrix to store results for plotting
vals = zeros(Complex, (length(a), length(N)))
baryvals = zeros(Complex, (length(a), length(N)))

# crunch integrals
for (i,a) in enumerate(a)
    for (j,N) in enumerate(N)
        vals[i,j] = cauchy_unit_circle(f, a, N)
        baryvals[i,j] = bary_unit_circle(f, a, N)
    end
end

# plot error
plt = plot(title="Convergence", xlabel="log(N)", ylabel="log(error)")
for (i, a) in enumerate(a)
    plot!(plt, N, abs.(vals[i,:] .- f(a)), label=string("a = ",imag(a), "i"), xaxis=:log, yaxis=:log, line=:dash)
    plot!(plt, N, abs.(baryvals[i,:] .- f(a)), label=string("a = ",imag(a), "i"), xaxis=:log, yaxis=:log)
end
display(plt)