# Test calling a Fortran code.
workspace()

aa = 1.8
bb = ccall((:__mymod_MOD_foo, "libmymod.so"), 
		Float64, (Ptr{Float64},), &aa)

xx = [1.2, 0.4]
yy = [0.0, 0.0]
ccall((:__mymod_MOD_add, "libmymod.so"), 
	Void, (Ptr{Float64}, Ptr{Float64}), xx, yy)
