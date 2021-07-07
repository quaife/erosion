using ScatteredInterpolation

# get grid
X = [collect(x) for x in collect(product(unique(flow.xtar),unique(flow.ytar)))]
drop = isBody.([xi[1] for xi in X[:]], [xi[2] for xi in X[:]], Ref(thlenden.thlenvec)) # drop points in bods
points = deleteat!(X[:], drop)
points = hcat(points...)
# get vals on grid
samples = hcat(deleteat!(u[:], drop), deleteat!(v[:], drop))
# interpolate
itp = ScatteredInterpolation.interpolate(ThinPlate(), points, samples)
# extrapolate
xtp = [evaluate(itp, p) for p in X]
u_itp = [xi[1] for xi in xtp]
v_itp = [xi[2] for xi in xtp]

# plot interp error on grid
u_err = abs.(u - u_itp)
u_err[drop] .= 0
v_err = abs.(v - v_itp)
v_err[drop] .= 0
plt = quiver(x, y, (u_err, v_err), aspect_ratio=1, transpose=true)
plot_bodies(plt, thlenden.thlenvec)
display(plt)

# plot extension
α = 0.1; τ = 10
x_mod = flow.xtar[drop]
y_mod = flow.ytar[drop]
u_itp_norm = u_itp # u_itp_norm = u_itp ./ sqrt.(u_itp.^2 + v_itp.^2)
u_itp_norm = α * u_itp_norm[drop]
v_itp_norm = v_itp # v_itp_norm = v_itp ./ sqrt.(u_itp.^2 + v_itp.^2)
v_itp_norm = α * v_itp_norm[drop]
quiver(x_mod[1:τ:length(x_mod)], y_mod[1:τ:length(x_mod)], quiver=(u_itp_norm[1:τ:length(x_mod)], v_itp_norm[1:τ:length(x_mod)]), transpose=true, aspect_ratio=1, xlims=(-1,1), ylims=(-1,1))
plt = contourf(x, y, u_itp, levels=15, transpose=true, aspect_ratio=1, xlims=(-1,1), ylims=(-1,1))
plot_bodies(plt, thlenden.thlenvec)
display(plt)







### MODIFIED ADVECT() FROM RUN1.JL TO USE TEST RBF INTERP FROM SCATTEREDINTERPOLATION.JL
