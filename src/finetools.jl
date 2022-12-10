using Random
using FinEtools
using Metis

"""
    make_topo_faces(t2v, e2t)

Make topological faces.

Returns
- `t2v`: the attribute `tfid` lists the numbers of the topological faces.
"""
function make_topo_faces(fens::FENodeSet, fes::E, crease_ang = 30/180*pi) where {E<:AbstractFESet} 
    t2v = to_core(fens, fes)
    orientedt2v, orientable = orient_surface_mesh(t2v)
    orientedt2v = make_topo_faces(orientedt2v)
    setlabel!(fes, [f for f in attribute(orientedt2v.left, "tfid")])
    return fens, fes
end

function to_core(fens::FENodeSet, fes::E) where {E<:AbstractFESet} 
    isa(fes, FinEtools.FESetModule.FESetT3) || error("Only 3-node triangles accepted")
    v = ShapeColl(P1, count(fens), "fens")
    N, T = size(fens.xyz, 2), eltype(fens.xyz)
    v.attributes["geom"] = VecAttrib([SVector{N, T}(fens.xyz[i, :]) for i in 1:size(fens.xyz, 1)])
    t = ShapeColl(T3, count(fes), "fes")
    IncRel(t, v, connasarray(fes))
end

function create_partitions(fens, fes, elem_per_partition = 50)
    fens, fes = make_topo_faces(fens, fes)
    surfids = deepcopy(fes.label)
    uniquesurfids = unique(surfids)
    nuniquesurfids = length(uniquesurfids)
    # This is the global partitioning of the surface elements
    partitionids = fill(0, count(fes))
    cpartoffset = 0
    for s in uniquesurfids
        el = selectelem(fens, fes, label = s)
        sfes = subset(fes, el)
        np = max(2, Int(round(length(el) / elem_per_partition)))
        femm1 = FEMMBase(IntegDomain(sfes, SimplexRule(2, 1)))
        C = dualconnectionmatrix(femm1, fens, 2)
        g = Metis.graph(C; check_hermitian=true)
        # @info "Surface $(s), $(np) partitions"
        spartitioning = Metis.partition(g, np; alg=:KWAY)
        for k in eachindex(el)
            partitionids[el[k]] = spartitioning[k] + cpartoffset
        end
        cpartoffset += length(unique(spartitioning))
    end 
    # Randomize the surface and partition ids
    pnumbers = unique(surfids)
    p = randperm(length(pnumbers))
    for k in eachindex(surfids)
        surfids[k] = p[surfids[k]]
    end
    pnumbers = unique(partitionids)
    p = randperm(length(pnumbers))
    for k in eachindex(partitionids)
        partitionids[k] = p[partitionids[k]]
    end
    # npanelgroups = length(unique(partitionids))
    return surfids, partitionids
end