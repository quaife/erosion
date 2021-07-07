module Yukawa

    export YukawaSolve, YukawaSolveBarycentric

    using FFTW
    using LinearAlgebra: norm, dot, I
    using SpecialFunctions: besselk
    using IterativeSolvers: gmres

    function YukawaSolve(dΩ, f, s)
            
        # geometry
        n̂, J, κ = Yukawa.GeomDerivs(dΩ)
        # solve SKIE
        σ, gmres_log = gmres(-I/2 + Yukawa.DLP(dΩ, s), -f.(dΩ, s); reltol = 1e-12, maxiter=128, log=true)
        # evaluate DLP at (x; s)
        Cʰ(x, s) = -sqrt(s)/length(dΩ) * sum(besselk.(1, sqrt(s) * abs.(x .- dΩ)) .* real(dot.(x .- dΩ, n̂)) ./ abs.(x .- dΩ) .* σ .* J)
        # construct C from Cp, Ch
        C(x, s) = Cʰ(x, s) + f(x, s)

        return C

    end

    function YukawaSolveBarycentric(dΩ, f, s)
        
        ℱ = x -> fftshift(fft(x))
        ℱ⁻¹ = x -> ifft(ifftshift(x))
        modes = N -> N % 2 == 0 ? [0; 1-N÷2:N÷2-1] : collect(-(N-1)÷2:(N-1)÷2)
        ∂ = fx -> im * modes(length(fx)) .* fx

        # geometry
        n̂, J, κ = GeomDerivs(dΩ)
        # solve SKIE
        σ, gmres_log = gmres(-I/2 + DLP(dΩ, s), -f.(dΩ, s); reltol = 1e-12, maxiter=128, log=true)
        # SMOOTH COMPONENT
        smooth(x, s) = besselk.(1, sqrt(s) * abs.(x .- dΩ)) + log.(x .- dΩ)
        trap(x, s) = -sqrt(s)/length(dΩ) * sum(smooth(x, s) .* real(dot.(x .- dΩ, n̂)) ./ abs.(x .- dΩ) .* σ .* J)

        # LOGARITHMIC COMPONENT
        σ′ = (ℱ⁻¹ ∘ ∂ ∘ ℱ)(σ)
        # boundary data on v(x), where v(x) defined as D[σ](x) = Re(v(x))
        v₀ = σ .- 1 / (2π * im) * (2π / length(dΩ)) * [sum(replace!(x -> isnan(x) ? σ₀′ : x, (σ₀ .- σ) ./ (x₀ .- dΩ))) for (x₀, σ₀, σ₀′) in zip(dΩ, σ, σ′)]
        # cauchy integral
        v⁺(x) = 1 / (2π*im) * 2π / length(dΩ) * sum(v₀ ./ (dΩ .- x) .* (im .* n̂ .* J))
        # a = sum(dΩ) / length(dΩ) # some point in Ωᶜ
        cauchy_one(x, a) = 1 / (x - a) / (2π*im) * 2π / length(dΩ) * sum(1 ./ (dΩ .- x) ./ (dΩ .- a) .* (im .* n̂ .* J))
        v(x) = v⁺(x) / cauchy_one(x, sum(dΩ) / length(dΩ))
        Cʰ(x, s) = trap(x, s) + real(v(x))

        # construct C from Cp, Ch
        C(x, s) = Cʰ(x, s) + f(x, s)

        return C

    end

    # compute normal vector, jacobian, curvature on dΩ
    function GeomDerivs(dΩ::Array{Complex{Float64}, 1})

        # shifted k vector for fourier diff
        k = (N = length(dΩ)) % 2 == 0 ? [0; 1-N/2:N/2-1] : collect(-(N-1)/2:(N-1)/2)
        # functions for fancy FFTs ;)
        ℱ = x -> fftshift(fft(x))
        ∂ = fx -> im * k .* fx
        ℱ⁻¹ = x -> ifft(ifftshift(x))

        # parametrizing dΩ as r(Θ), crunch derivatives in fourier space
        r′ = (ℱ⁻¹ ∘ ∂ ∘ ℱ)(dΩ)
        r′′ = (ℱ⁻¹ ∘ ∂ ∘ ℱ)(r′)
        # jacobian
        J = norm.(r′)
        # unit outward normal
        n̂ = -im * r′ ./ J
        # curvature
        κ = imag(conj.(r′) .* r′′) ./ J.^3

        return n̂, J, κ
    end

    # generate double layer potential matrix
    function DLP(dΩ::Array{Complex{Float64}, 1}, s::Complex{Float64})

        N = length(dΩ)
        
        # derivatives at boundary points
        n̂, J, κ = GeomDerivs(dΩ)
        
        # DLP off-diagonal terms
        K(i, j) = -sqrt(s) * besselk(1, sqrt(s) * abs(dΩ[i] - dΩ[j])) * real((dot(dΩ[i] - dΩ[j], n̂[j]))) / abs(dΩ[i] - dΩ[j])
        K_diag(i) = 1/2 * κ[i]
        
        return 2π/N * [i == j ? 1/2π * complex(K_diag(i)) * J[i] : 1/2π * complex(K(i, j)) * J[j] for i in 1:N, j in 1:N]
    end

    function solveSKIE(dΩ, f, s)
        return gmres(-I/2 + DLP(dΩ, s), f.(dΩ, s); reltol = 1e-12, maxiter=20, log=true)
    end

end