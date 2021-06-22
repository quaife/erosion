# basic.jl: Basic routines such as datatypes.
using DelimitedFiles

# readvec: Read a vector from a text file.
function readvec(filename::AbstractString)
	iostream = open(filename, "r")
	invec = readdlm(iostream, comments=true)[:,1]
	close(iostream)
	return invec
end

# ThetaLenType: Includes the geometry data and stress of a single body.
mutable struct ThetaLenType
	theta::Vector{Float64}; len::Float64; xsm::Float64; ysm::Float64; matau::Float64
end

# ThLenDenType: Includes the vector of all thlens and the density function.
mutable struct ThLenDenType
	thlenvec::Vector{ThetaLenType}; tt::Float64;
	density::Vector{Float64}; denrot::Vector{Float64}; 
end

# Create a new ThLenDenType variable.
new_thlenden(thlenvec::Vector{ThetaLenType}) = ThLenDenType(thlenvec, 0., [],[])