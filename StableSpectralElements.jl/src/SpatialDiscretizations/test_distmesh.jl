using StableSpectralElements
using DistMesh
fd(pp) = drectangle(pp, -1, 1, -1, 1)
mesh = distmesh2d(fd, huniform, 0.2, ((-1,-1), (1,1)))
# mesh is a struct
#- `p::Vector{Point2d}`: The node positions.
#- `t::Vector{Index3}`: The triangle connectivity indices.

println(length(mesh.p))
println(length(mesh.t))

VX = zeros(length(mesh.p))
VY = zeros(length(mesh.p))
EToV = zeros(Int, length(mesh.p), 3)
for i in 1:length(mesh.p)
    VX[i] = mesh.p[i][1]
    VY[i] = mesh.p[i][2]
end

for i in 1:length(mesh.t)
    EToV[i,:] = mesh.t[i]
end