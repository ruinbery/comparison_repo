"""
    RefElemData(elementType::Line, approxType::SBP, N)
    RefElemData(elementType::Quad, approxType::SBP, N)
    RefElemData(elementType::Hex,  approxType::SBP, N)
    RefElemData(elementType::Tri,  approxType::SBP, N)
    
SBP reference element data for `Quad()`, `Hex()`, and `Tri()` elements. 

For `Line()`, `Quad()`, and `Hex()`, `approxType` is `SBP{TensorProductLobatto}`.

For `Tri()`, `approxType` can be `SBP{Kubatko{LobattoFaceNodes}}`, `SBP{Kubatko{LegendreFaceNodes}}`, or `SBP{Hicken}`. 
"""
function RefElemData(elementType::Line, approxType::SBP{TensorProductLobatto}, N; tol = 100*eps(), kwargs...)

    rd = RefElemData(elementType, N; quad_rule_vol = gauss_lobatto_quad(0,0,N), kwargs...)        
    
    rd = @set rd.Vf = droptol!(sparse(rd.Vf), tol)
    rd = @set rd.LIFT = Diagonal(rd.wq) \ (rd.Vf' * Diagonal(rd.wf)) # TODO: make this more efficient with LinearMaps?

    return _convert_RefElemData_fields_to_SBP(rd, approxType)
end

function RefElemData(elementType::Union{Quad, Hex}, approxType::SBP{TensorProductLobatto}, N; tol = 100*eps(), kwargs...)

    rd = RefElemData(elementType, 
                     Polynomial(TensorProductQuadrature(gauss_lobatto_quad(0, 0, N))), 
                     N; kwargs...)

    # brute-force determine Fmask so that rd.rf = rd.r[rd.Fmask], etc.
    new_Fmask = copy(rd.Fmask)    
    for fid in eachindex(rd.rf)
        rstf = SVector(getindex.(rd.rstf, fid))
        for i in eachindex(rd.r)
            rst = SVector(getindex.(rd.rst, i))
            if rstf ≈ rst
                new_Fmask[fid] = i
                break
            end
        end
    end
    rd = @set rd.Fmask = new_Fmask;                     

    rd = @set rd.Vf = droptol!(sparse(rd.Vf), tol)
    rd = @set rd.LIFT = Diagonal(rd.wq) \ (rd.Vf' * Diagonal(rd.wf)) # TODO: make this more efficient with LinearMaps?

    return _convert_RefElemData_fields_to_SBP(rd, approxType)
end

function RefElemData(elementType::Tri, approxType::SBP, N; tol = 100*eps(), kwargs...)
    
    quad_rule_vol, quad_rule_face = diagE_sbp_nodes(elementType, approxType, N)

    # build polynomial reference element using quad rules; will be modified to create SBP RefElemData
    rd = RefElemData(elementType, Polynomial(), N; quad_rule_vol=quad_rule_vol, quad_rule_face=quad_rule_face, kwargs...)

    # determine Fmask = indices of face nodes among volume nodes
    Ef, Fmask = build_Ef_Fmask(rd)

    # Build traditional SBP operators from hybridized operators. See Section 3.2 of 
    # [High-order entropy stable dG methods for the SWE](https://arxiv.org/pdf/2005.02516.pdf)
    # by Wu and Chan 2021. [DOI](https://doi.org/10.1016/j.camwa.2020.11.006)
    (Qrh, Qsh), _ = hybridized_SBP_operators(rd)
    Nq = length(rd.wq)
    Vh_sbp = [I(Nq); Ef]
    Qr = Vh_sbp' * Qrh * Vh_sbp
    Qs = Vh_sbp' * Qsh * Vh_sbp
    Dr,Ds = (x -> diagm(1 ./ rd.wq) * x).((Qr, Qs))
    
    rd = @set rd.rst = quad_rule_vol[1:2]   # set nodes = SBP nodes
    rd = @set rd.rstq = quad_rule_vol[1:2]  # set quad nodes = SBP nodes
    rd = @set rd.Drst = (Dr, Ds)
    rd = @set rd.Fmask = vec(Fmask)

    # TODO: make these more efficient with custom operators?
    rd = @set rd.Vf = droptol!(sparse(Ef), tol)
    rd = @set rd.LIFT = Diagonal(rd.wq) \ (rd.Vf' * Diagonal(rd.wf)) 

    # make V1 the interpolation matrix from triangle vertices to SBP nodal points
    rd = @set rd.V1 = vandermonde(elementType, N, rd.rst...) / rd.VDM * rd.V1

    # Vp operator = projects SBP nodal vector onto degree N polynomial, then interpolates to plotting points
    rd = @set rd.Vp = vandermonde(elementType, N, rd.rstp...) / rd.VDM * rd.Pq

    return _convert_RefElemData_fields_to_SBP(rd, approxType)
end

# note: written with ChatGPT
function RefElemData(elementType::Tet, approxType::SBP{WHZ}, N; tol = 100*eps(), kwargs...)

    # get SBP volume and facet quadrature nodes
    quad_rule_vol, quad_rule_face = diagE_sbp_nodes(elementType, approxType, N)
    r, s, t, w = quad_rule_vol
    rf, sf, tf = quad_rule_face  # face quadrature coordinates

    # build polynomial reference element using quad rules
    rd = RefElemData(elementType, Polynomial(), N; 
                     quad_rule_vol = (r, s, t, w), 
                     quad_rule_face = (rf, sf, tf), 
                     kwargs...)

    # determine Fmask and Ef for tetrahedron
    Ef, Fmask = build_Ef_Fmask_tet(rd; tol=tol)

    # build hybridized SBP operators for 3D (Dr, Ds, Dt)
    (Qrh, Qsh, Qth), _ = hybridized_SBP_operators(rd)  # make sure this returns 3D
    Nq = length(rd.wq)
    Vh_sbp = [I(Nq); Ef]
    Qr = Vh_sbp' * Qrh * Vh_sbp
    Qs = Vh_sbp' * Qsh * Vh_sbp
    Qt = Vh_sbp' * Qth * Vh_sbp
    Dr, Ds, Dt = (x -> diagm(1 ./ rd.wq) * x).((Qr, Qs, Qt))

    # assign SBP nodes and operators
    rd = @set rd.rst = (r, s, t)
    rd = @set rd.rstq = (r, s, t)
    rd = @set rd.Drst = (Dr, Ds, Dt)
    rd = @set rd.Fmask = vec(Fmask)

    # assign extraction matrix and lifting operator
    rd = @set rd.Vf = droptol!(sparse(Ef), tol)
    rd = @set rd.LIFT = Diagonal(rd.wq) \ (rd.Vf' * Diagonal(rd.wf))

    # V1 = interpolation from vertices to SBP nodal points
    rd = @set rd.V1 = vandermonde(elementType, N, rd.rst...) / rd.VDM * rd.V1

    # Vp = project SBP nodal vector onto degree N polynomial, then interpolate to plotting points
    rd = @set rd.Vp = vandermonde(elementType, N, rd.rstp...) / rd.VDM * rd.Pq

    return _convert_RefElemData_fields_to_SBP(rd, approxType)
end


#####
##### Utilities for SBP 
#####

# - HDF5 file created using MAT.jl and the following code:
# vars = matread("src/data/sbp_nodes/KubatkoQuadratureRules.mat")
# h5open("src/data/sbp_nodes/KubatkoQuadratureRules.h5", "w") do file
#     for qtype in ("Q_GaussLobatto", "Q_GaussLegendre")
#         group = create_group(file, qtype) # create a group
#         for fieldname in ("Points", "Domain", "Weights")
#             subgroup = create_group(group, fieldname)
#             for N in 1:length(vars[qtype])
#                 subgroup[string(N)] = vars[qtype][N][fieldname]
#             end
#         end
#     end
# end

function diagE_sbp_nodes(elem::Tri, approxType::SBP{Kubatko{LobattoFaceNodes}}, N)    
    
    if N==6
        @warn "N=6 SBP operators with quadrature strength 2N-1 and Lobatto face nodes may require very small timesteps."
    end
    if N > 6
        @error "N > 6 triangular `SBP{Kubatko{LobattoFaceNodes}}` operators are not available."
    end

    # from Ethan Kubatko, private communication
    vars = h5open((@__DIR__) * "/data/sbp_nodes/KubatkoQuadratureRules.h5", "r")
    rs = vars["Q_GaussLobatto"]["Points"][string(N)][]
    r, s = (rs[:, i] for i = 1:size(rs, 2))
    w = vec(vars["Q_GaussLobatto"]["Weights"][string(N)][])
    quad_rule_face = gauss_lobatto_quad(0, 0, N+1)     

    return (r, s, w), quad_rule_face 
end

function diagE_sbp_nodes(elem::Tri, approxType::SBP{Kubatko{LegendreFaceNodes}}, N)    

    if N > 6
        @error "N > 6 triangular `SBP{Kubatko{LegendreFaceNodes}}` operators are not available."
    end

    # from Ethan Kubatko, private communication
    vars = h5open((@__DIR__) * "/data/sbp_nodes/KubatkoQuadratureRules.h5", "r")
    rs = vars["Q_GaussLegendre"]["Points"][string(N)][]
    r, s = (rs[:, i] for i = 1:size(rs, 2))
    w = vec(vars["Q_GaussLegendre"]["Weights"][string(N)][])
    quad_rule_face = gauss_quad(0, 0, N)

    return (r, s, w), quad_rule_face 
end

parsevec(type, str) = str |>
  (x -> split(x, ", ")) |>
  (x -> map(y -> parse(type, y), x))
  
function diagE_sbp_nodes(elem::Tri, approxType::SBP{Hicken}, N)    
    
    if N > 4
        @error "N > 4 triangular `SBP{Hicken}` operators are not available."
    end

    # from Jason Hicken https://github.com/OptimalDesignLab/SummationByParts.jl/tree/work
    lines = readlines((@__DIR__)*"/data/sbp_nodes/tri_diage_p$N.dat") 
    r = parsevec(Float64,lines[11])
    s = parsevec(Float64,lines[12])
    w = parsevec(Float64,lines[13])

    # convert Hicken format to biunit right triangle
    r = @. 2*r-1 
    s = @. 2*s-1
    w = 2.0 * w/sum(w)

    quad_rule_face = gauss_lobatto_quad(0,0,N+1) 

    return (r,s,w), quad_rule_face 
end

# Make two new functions based on quadrature nodes from Worku, Hicken, and Zingg (WHZ): https://arxiv.org/abs/2311.15576
# data storied in: https://github.com/OptimalDesignLab/SummationByParts.jl/tree/master/quadrature_data

function diagE_sbp_nodes(elem::Tri, approxType::SBP{WHZ},N)

    if N > 10
        @error "deg N > 10 triangular `SBP{WHZ}` operators are not available."
    end

    q = 2*N  # same as in TPSS paper

    lines = readlines(string(@__DIR__, "/data/sbp_nodes/WHZ_quadrature/tri_q", q, ".dat"))
    rsw = [parse.(Float64, split(line, ',')) for line in lines]
    rsw = hcat(rsw...)'   # each column is a row of the file

    r = rsw[:, 1]
    s = rsw[:, 2]
    w = rsw[:, 3]

    # these quad rules are already defined on the biunit right triangle
 
    # facet quadrature rule is defined as below. formula for N determined by looking at
    # what gauss_lobatto_quad(0,0,N) produces versus the points in the tri_q#.dat file
    # ie. making sure the points line up
    quad_rule_face = gauss_lobatto_quad(0,0,floor(Int,(q+1)/2)+1) 

    return (r,s,w), quad_rule_face 
end

function diagE_sbp_nodes(elem::Tet, approxType::SBP{WHZ},N)

    if N > 5
        @error "deg N > 5 Tetrahedra `SBP{WHZ}` operators are not available."
    end

    q = 2*N  # same as in TPSS paper

    lines = readlines(string(@__DIR__, "/data/sbp_nodes/WHZ_quadrature/tet_q", q, ".dat"))
    rstw = [parse.(Float64, split(line, ',')) for line in lines]
    rstw = hcat(rstw...)'  

    r = rstw[:, 1]
    s = rstw[:, 2]
    t = rstw[:, 3]
    w = rstw[:, 4]

    # these quad rules are already defined on the biunit right triangle

    lines = readlines(string(@__DIR__, "/data/sbp_nodes/WHZ_quadrature/tet_q", q, "_facet.dat"))
    xyb = [parse.(Float64, split(line, ',')) for line in lines]
    xyb = hcat(xyb...)'   

    x = xyb[:, 1]
    y = xyb[:, 2]
    b = xyb[:, 3]

    return (r,s,t,w), (x,y,b)
end

function build_Ef_Fmask(rd_sbp::RefElemData; tol = 100*eps())
    
    (; rq, sq, rf, sf, Nfaces ) = rd_sbp   
    rf, sf = (x->reshape(x, length(rf) ÷ Nfaces, Nfaces)).((rf, sf))
    Fmask = zeros(Int, length(rf) ÷ Nfaces, Nfaces) # 
    Ef = zeros(length(rf), length(rq)) # extraction matrix
    for i in eachindex(rq)
        for f = 1:rd_sbp.Nfaces
            id = findall(@. abs(rq[i]-rf[:,f]) + abs(sq[i]-sf[:,f]) .< tol)
            Fmask[id,f] .= i
            Ef[id .+ (f-1) * size(rf,1), i] .= 1
        end
    end
    return Ef, Fmask
end

# note: this is written with ChatGPT
function build_Ef_Fmask_tet(rd_sbp::RefElemData; tol = 100*eps())
    # extract node coordinates and number of faces
    (; rq, sq, tq, rf, sf, tf, Nfaces) = rd_sbp  

    # reshape face node coordinates into (nodes_per_face, Nfaces)
    rf, sf, tf = (x -> reshape(x, length(x) ÷ Nfaces, Nfaces)).((rf, sf, tf))

    nodes_per_face = size(rf, 1)
    Nvol = length(rq)

    # initialize face mask and extraction matrix
    Fmask = zeros(Int, nodes_per_face, Nfaces)
    Ef = zeros(nodes_per_face * Nfaces, Nvol)

    # loop over volume nodes and faces
    for i in 1:Nvol
        for f in 1:Nfaces
            # find matching nodes on this face within tolerance
            id = findall(@. abs(rq[i]-rf[:,f]) + abs(sq[i]-sf[:,f]) + abs(tq[i]-tf[:,f]) < tol)
            if !isempty(id)
                Fmask[id, f] .= i
                Ef[id .+ (f-1)*nodes_per_face, i] .= 1
            end
        end
    end

    return Ef, Fmask
end

get_face_nodes(x::AbstractVector, Fmask) = view(x, Fmask)
get_face_nodes(x::AbstractMatrix, Fmask) = view(x, Fmask, :)

function _convert_RefElemData_fields_to_SBP(rd, approx_type::SBP)
    rd = @set rd.M = Diagonal(rd.wq)
    rd = @set rd.Pq = I
    rd = @set rd.Vq = I
    rd = @set rd.approximation_type = approx_type
    return rd
end

"""
    function hybridized_SBP_operators(rd::RefElemData{DIMS}) 

Constructs hybridized SBP operators given a `RefElemData`. Returns operators `Qrsth..., VhP, Ph`.
"""
function hybridized_SBP_operators(rd)
    (; M, Vq, Pq, Vf, wf, Drst, nrstJ ) = rd
    Qrst = (D->Pq' * M * D * Pq).(Drst)
    Ef = Vf * Pq
    Brst = (nJ->diagm(wf .* nJ)).(nrstJ)
    Qrsth = ((Q, B)->.5*[Q-Q' Ef'*B; -B*Ef B]).(Qrst, Brst)
    Vh = [Vq; Vf]
    Ph = M \ transpose(Vh)
    VhP = Vh * Pq
    return Qrsth, VhP, Ph, Vh
end



