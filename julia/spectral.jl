# spectral.jl: Collection of spectral routines 

# fftnice: shifted fft.
function fftnice(fx::Vector)
	fh = fft(fx)
	return fftshift(fh)
end
# ifftnice: shifted ifft.
function ifftnice(fh::Vector)
	fhsh = ifftshift(fh)
	return ifft(fhsh)
end
# kvec: Construct the vector of k values.
# In the even case, sets the highest mode to zero or not depending on hmode.
function kvec(nn::Integer, hmode=0)
	if iseven(nn)
		nsm = div(nn,2)-1
		return [hmode*(nsm+1); -nsm : nsm]
	else
		nsm = div(nn-1,2)
		return -nsm : nsm
	end
end
#= specdiff: Compute the derivative of fx spectrally.
Default: assumes that the length of the interval is 1.0 =#
function specdiff(fx::Vector, intvlen::Float64=2*pi)	
	fh = fftnice(fx)
	kv = kvec(length(fx), 0)
	dfh = 2*pi*im * kv.*fh
	df = ifftnice(dfh)
	imagtest(df)
	return real(df)/intvlen
end
#= specint: Compute the antiderivative of fx spectrally.
Assumes that the input fx has mean-zero and 
forces the output Fx to have mean-zero. =#
function specint(fx::Vector, intvlen::Float64=2*pi)
	fh = fftnice(fx)
	kv = kvec(length(fx), 0)
	Fh = fh./(2*pi*im*kv)
	# This kills the zeroth mode and, if nn is even, also the highest mode.
	Fh[kv.==0] = 0.
	Fx = ifftnice(Fh)
	imagtest(Fx)
	return real(Fx)*intvlen	
end
# imagtest: Test that the imaginary part is negligible.
function imagtest(fx::Vector, relthold::Float64=1e-8)
    maximrel = maxabs(imag(fx))/maxabs(fx)
    if maximrel > relthold
        warn("imag part too big: ", maximrel)
    end
    return
end


# expsmooth: Equivalent to a Gaussian filter.
function expsmooth(fx::Vector, factor::Float64)
    fh = fftnice(fx)
    kv = kvec(length(fx), 1)
    fh .*= exp(-factor*abs(kv).^2)
    fs = ifftnice(fh)
    imagtest(fs)
    return real(fs) 
end
# expsmooth(fx, -0.5*sigma^2) = gaussfilter(fx, sigma)



#------------------OBSELETE BELOW------------------#
#=
# gaussfilter: Apply a Gaussian filter of width sigma.
function gaussfilter(fx::Vector, sigma::Float64)
   fh = fftnice(fx)
   kv = kvec(length(fx), 1)
   fh .*= exp(-0.5*sigma^2 * abs(kv).^2)
   fs = ifftnice(fh)
   return real(fs)
end
# krasnyfilter: Apply a Krasny filter to the spectrum.
# The Krasny filter does not delay the shape instability, so I voided it.
function krasnyfilter(fx::Vector, relthold::Float64=1e-8)
	fh = fftnice(fx)
	absthold = relthold*maxabs(fh)
	fh[abs(fh) .< absthold] = 0.0
	fx = ifftnice(fh)
	imagtest(fx)
	return real(fx)
end 
=#
