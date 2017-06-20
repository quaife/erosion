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
# imagtest: Test that the imaginary part is negligible.
function imagtest(fx::Vector, thold::Float64=1e-10)
	maxim = maxabs(imag(fx))
	if maxim > thold
		warn("imag part too big: ", maxim)
	end
	return
end
#= specdiff: Compute the derivative of fx spectrally.
Default: assumes that the length of the interval is 1.0 =#
function specdiff(fx::Vector, intvlen::Float64=1.0)	
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
function specint(fx::Vector, intvlen::Float64=1.0)
	fh = fftnice(fx)
	kv = kvec(length(fx), 0)
	Fh = fh./(2*pi*im*kv)
	# This kills the zeroth mode and, if nn is even, also the highest mode.
	Fh[kv.==0] = 0.
	Fx = ifftnice(Fh)
	imagtest(Fx)
	return real(Fx)*intvlen	
end
# gaussfilter: Apply a Gaussian filter of width sigma.
function gaussfilter(fx::Vector, sigma::Float64)
	fh = fftnice(fx)
	kv = kvec(length(fx), 1)
	fh .*= exp(-0.5*sigma^2 * abs(kv).^2)
	fs = ifftnice(fh)
	return real(fs)
end

#function krasnyfilter!(uu, thold)
#	uu[abs(uu)<thold] = zero(eltype(uu))
#end