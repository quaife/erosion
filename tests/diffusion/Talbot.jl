module Talbot

    export TalbotContour, talbot_midpoint, get_talbot_points, get_talbot_derivative

    """
    A curve of the form z(θ) = σ + μ(cotθ + iνθ), for θ in (-π, π).
    """
    struct TalbotContour
        z::Function
        z′::Function
        TalbotContour(z::Function, z′::Function) = new(
            (θ; t=nothing) -> z(θ; t=t),
            (θ; t=nothing) -> z′(θ; t=t)
        )
        TalbotContour(σ::Float64, μ::Float64, ν::Float64) = new(
            (θ; t=nothing) -> σ .+ μ * (θ .* cot.(θ) + im * ν .* θ),
            (θ; t=nothing) -> μ * (cot.(θ) - θ .* csc.(θ).^2 .+ im * ν)
        )
        TalbotContour(σ::Function, μ::Function, ν::Function) = new(
            (θ; t=nothing) -> σ(θ; t=t) .+ μ(θ; t=t) .* (θ .* cot.(θ) + im * ν(θ; t=t) .* θ),
            (θ; t=nothing) -> μ(θ; t=t) * (cot.(θ) - θ .* csc.(θ).^2 .+ im * ν(θ; t=t))
        )
    end

    """
        talbot_midpoint(Y, s, N)

    Construct y(t), the inverse laplace transform of Y(s) using the N-point
    midpoint rule on the Talbot contour s(θ; t) parametrized by θ ∈ [-π,π].
    """
    function talbot_midpoint(Y::Function, s::TalbotContour, N)
        dθ = 2π/N; θ = -π+dθ/2:dθ:π-dθ/2
        return t -> -im / N * sum(exp.(s.z(θ; t=t)*t) .* Y.(s.z(θ; t=t)) .* s.z′(θ; t=t))
    end

    function γ(s::TalbotContour, N, t)
        dθ = 2π/N; θ = -π+dθ/2:dθ:π-dθ/2
        return s.z(θ; t=t)
    end

    function γ′(s::TalbotContour, N, t)
        dθ = 2π/N; θ = -π+dθ/2:dθ:π-dθ/2
        return s.z′(θ; t=t)
    end

end