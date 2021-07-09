# OBJECTIVE: Collection of spectral routines 

module SpectralMethods

	export specdiff, specint, expsmooth, gaussfilter

	using FFTW: fft, ifft, fftshift, ifftshift

	# Shifted fft and ifft.
	fftnice(fx::Vector) = fftshift( fft(fx) )
	ifftnice(fh::Vector) = ifft( ifftshift(fh) )

	# Construct the vector of k values.
	# In the even case, sets the highest mode to zero or not depending on hmode.
	function kvec(nn::Integer, hmode = 0)
		if iseven(nn)
			nsm = div(nn,2)-1
			return [hmode*(nsm+1); -nsm : nsm]
		else
			nsm = div(nn-1,2)
			return -nsm : nsm
		end
	end

	# Test that the imaginary part is negligible.
	function imagtest(fx::Vector, relthold::Float64 = 1e-8)
		maximrel = maximum(abs.(imag(fx)))/maximum(abs.(fx))
		if maximrel > relthold
			@warn string("imag part too big: ", maximrel)
		end
	end

	# Compute the derivative of fx spectrally.
	# Default: assumes that the length of the interval is 1.0.
	function specdiff(fx::Vector, intvlen::Float64 = 2*pi)	
		fh = fftnice(fx)
		kv = kvec(length(fx), 0)
		df = ifftnice(2*pi*im * kv.*fh)
		imagtest(df)
		return real(df)/intvlen
	end

	# Compute the antiderivative of fx spectrally.
	# Assumes that the input fx has mean-zero and forces the output Fx to have mean-zero.
	function specint(fx::Vector, intvlen::Float64 = 2*pi)
		fh = fftnice(fx)
		kv = kvec(length(fx), 0)
		Fh = fh./(2*pi*im*kv)
		# This kills the zeroth mode and, if nn is even, also the highest mode.
		Fh[kv.==0] .= 0.
		Fx = ifftnice(Fh)
		imagtest(Fx)
		return real(Fx)*intvlen	
	end

	# Smooth the function fx in Fourier space by a given factor.
	function expsmooth(fx::Vector, factor::Float64)
		fh = fftnice(fx)
		kv = kvec(length(fx), 1)
		fh .*= exp.(-factor*abs.(kv).^2)
		fs = ifftnice(fh)
		imagtest(fs)
		return real(fs) 
	end

	# Apply a Gaussian filter of width sigma.
	gaussfilter(fx::Vector, sigma::Float64) = expsmooth(fx, 0.5*sigma^2)

end