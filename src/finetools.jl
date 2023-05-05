using Random
using FinEtools
using Metis
using LinearAlgebra
using FinEtools.MeshExportModule: VTK

"""
    make_topo_faces(t2v, e2t)

Make topological faces.

Returns
- `t2v`: the attribute `tfid` lists the numbers of the topological faces.
"""
function make_topo_faces(fens::FENodeSet, fes::E, crease_ang = 30/180*pi) where {E<:AbstractFESet} 
    t2v = _to_core(fens, fes)
    orientedt2v, orientable = orient_surface_mesh(t2v)
    orientedt2v = make_topo_faces(orientedt2v, crease_ang)
    setlabel!(fes, [f for f in attribute(orientedt2v.left, "tfid")])
    return fens, fes
end

function _to_core(fens::FENodeSet, fes::E) where {E<:AbstractFESet} 
    isa(fes, FinEtools.FESetModule.FESetT3) || error("Only 3-node triangles accepted")
    v = ShapeColl(P1, count(fens), "fens")
    N, T = size(fens.xyz, 2), eltype(fens.xyz)
    v.attributes["geom"] = VecAttrib([SVector{N, T}(fens.xyz[i, :]) for i in 1:size(fens.xyz, 1)])
    t = ShapeColl(T3, count(fes), "fes")
    IncRel(t, v, connasarray(fes))
end

"""
    create_partitions(fens, fes, elem_per_partition = 50; crease_ang = 30/180*pi, cluster_max_normal_deviation = 2 * crease_ang)

Create partitions of the triangulated boundary into clusters.

# Input

- `fens`, `fes` = finite element mesh,
- `elem_per_partition` = desired number of elements per partition,
- `crease_ang` = crease angle to determine boundaries between topological faces,
- `cluster_max_normal_deviation` = maximum deviation of the normal within the cluster.

# Output

- `surfids` = array of surface identifiers, one for each boundary facet
- `partitionids` = array of partition identifiers (i.e. cluster identifiers),  one for each boundary facet
- `surface_elem_per_partition` = dictionary of cluster sizes, indexed by the surface id
"""
function create_partitions(fens, fes, elem_per_partition = 50; crease_ang = 30/180*pi, cluster_max_normal_deviation = 2 * crease_ang)
    fens, fes = make_topo_faces(fens, fes, crease_ang)
    surfids = deepcopy(fes.label)
    uniquesurfids = unique(surfids)
    normals = _compute_normals(fens, fes)
    # This is the global partitioning of the surface elements
    partitionids = fill(0, count(fes))
    surface_elem_per_partition = Dict()
    cpartoffset = 0
    for surf in uniquesurfids
        el = selectelem(fens, fes, label = surf)
        spartitioning, adjusted_elem_per_partition = _partition_surface(fens, fes, surf, el, elem_per_partition, cluster_max_normal_deviation, normals)
        for k in eachindex(el)
            partitionids[el[k]] = spartitioning[k] + cpartoffset
        end
        cpartoffset += length(unique(spartitioning))
        surface_elem_per_partition[surf] = adjusted_elem_per_partition
    end 
    # Randomize the surface ids
    prmtd = randperm(length(uniquesurfids))
    for k in eachindex(surfids)
        surfids[k] = prmtd[surfids[k]]
    end
    permuted_surface_elem_per_partition = Dict()
    for k in uniquesurfids
        permuted_surface_elem_per_partition[prmtd[k]] = surface_elem_per_partition[k]
    end
    surface_elem_per_partition = permuted_surface_elem_per_partition
    permuted_surface_elem_per_partition = nothing
    # Randomize the partition ids
    pnumbers = unique(partitionids)
    prmtd = randperm(length(pnumbers))
    for k in eachindex(partitionids)
        partitionids[k] = prmtd[partitionids[k]]
    end
    # npanelgroups = length(unique(partitionids))
    return surfids, partitionids, surface_elem_per_partition
end


function _compute_normals(fens, fes)
    sn = SurfaceNormal(size(fens.xyz, 2))
    parametric_centroid = centroidparametric(fes)
    N = bfun(fes, parametric_centroid);
    gradNpar = bfundpar(fes, parametric_centroid);
    conns = connasarray(fes)
    normals = fill(0.0, count(fes), 3)
    for j in 1:count(fes)
        c = N' * view(fens.xyz, conns[j, :], :);
        J = view(fens.xyz, conns[j, :], :)' * gradNpar;
        n = updatenormal!(sn, c, J, 0)
        n ./= norm(n);
        normals[j, :] = n
    end 
    return normals
end

function _clusters_sufficiently_flat(fens, fes, el, surf, normals, spartitioning, cluster_max_normal_deviation)
    upids = unique(spartitioning)
    cmnd = cos(cluster_max_normal_deviation)
    # Now go through the partitions
    for p in upids
        pel = [el[k] for k in eachindex(spartitioning) if spartitioning[k] == p]
        for i in 1:length(pel)
            for j in i+1:length(pel)
                if dot(view(normals, pel[i], :), view(normals, pel[j], :)) < cmnd
                    return false
                end
            end
        end
    end
    return true
end

function _estimate_number_of_partitions(nel, np)
    fraction = 0.6
    nc = Int(floor(nel / np))
    nc0 = nc
    while (nel - nc * np < fraction * nc) || !(nc < nc0)
        np += 1
        nc = Int(floor(nel / np))
        if np == nel 
            break
        end
    end
    return min(nel, np)
end

function _partition_surface(fens, fes, surf, el, elem_per_partition, cluster_max_normal_deviation, normals)
    sfes = subset(fes, el)
    nel = count(sfes)
    femm1 = FEMMBase(IntegDomain(sfes, SimplexRule(2, 1)))
    C = dualconnectionmatrix(femm1, fens, 2)
    g = Metis.graph(C; check_hermitian=true)
    np = min(nel, Int(ceil(nel / elem_per_partition)))
    np = _estimate_number_of_partitions(nel, np)
    while true
        elem_per_partition = Int(round(nel / np))
        if np > 1
            spartitioning = Metis.partition(g, np; alg=:KWAY)
        else
            spartitioning = ones(Int32, nel)
        end
        upids = sort(unique(spartitioning))
        @assert upids[1] > 0
        if _clusters_sufficiently_flat(fens, fes, el, surf, normals, spartitioning, cluster_max_normal_deviation)
            # Make the partitioning numbers contiguous
            newpids = fill(0, maximum(upids))
            p = 1
            for k in eachindex(upids)
                newpids[upids[k]] = p
                p += 1
            end
            for j in eachindex(spartitioning)
                spartitioning[j] = newpids[spartitioning[j]]
            end
            return spartitioning, elem_per_partition
        else
            if elem_per_partition  <= 1
                return collect(1:length(el)), 1
            end
        end
        np = _estimate_number_of_partitions(nel, np+1)
    end
end
