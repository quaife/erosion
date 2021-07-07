module Fourier

    export FSolve, modes

    using LinearAlgebra: norm
    using SpecialFunctions: besselk, hankelh2
    using FFTW

    # functions for fancy FFTs
    ℱ = x -> fftshift(fft(x))
    ℱ⁻¹ = x -> ifft(ifftshift(x))
    modes = N -> N % 2 == 0 ? [0; 1-N÷2:N÷2-1] : collect(-(N-1)÷2:(N-1)÷2)
    ∂ = fx -> im * modes(length(fx)) .* fx

    function FSolve(N, x₀, r, s)

        # boundary discretization
        θ = 2π/N * (N % 2 == 0 ? [-N÷2; modes(N-1)] : modes(N))
        dΩ = exp.(im * θ);
        # boundary condition in laplace space
        f(x, s) = sum([1/2π * besselk(0, sqrt(s) * norm(x - x₀)) for x₀ in x₀])
        # solve ode on boundary
        u(r, n, s) = hankelh2(n, -im * r * sqrt(s))
        c₁(s) = ℱ(f.(dΩ, s)) ./ u.(abs.(dΩ), modes(N), s)
        # back into laplace space
        Cʰ(r, s) = ℱ⁻¹(c₁(s) .* u.(r, modes(N), s))

        return C(r, s) = Cʰ(r, s) .- f.(r * exp.(im * θ), s)

    end

end