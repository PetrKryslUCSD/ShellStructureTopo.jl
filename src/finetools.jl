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
    orientedt2v = make_topo_faces(orientedt2v)
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

function create_partitions(fens, fes, elem_per_partition = 50; min_num_partitions = 4, max_normal_deviation = pi/2)
    fens, fes = make_topo_faces(fens, fes)
    surfids = deepcopy(fes.label)
    uniquesurfids = unique(surfids)
    nuniquesurfids = length(uniquesurfids)
    normals = _compute_normals(fens, fes)
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
        spartitioning = Metis.partition(g, np; alg=:KWAY)
        # Now we have the surface partitioning based purely on topology. The
        # next thing is to consider the geometry so that the partitions are
        # reasonably flat.
        spartitioning = _adjust_for_curvature(fens, fes, el, s, normals, spartitioning, elem_per_partition, min_num_partitions, max_normal_deviation)
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


function _adjust_for_curvature(fens, fes, el, surf, normals, spartitioning, elem_per_partition, min_num_partitions, max_normal_deviation)
    upids = unique(spartitioning)
    cmnd = cos(max_normal_deviation)
    # Now go through the partitions
    for p in upids
        pel = [el[k] for k in eachindex(spartitioning) if spartitioning[k] == p]
        pass = true
        for i in 1:length(pel)
            for j in i+1:length(pel)
                if dot(view(normals, pel[i], :), view(normals, pel[j], :)) < cmnd
                    pass = false; break
                end
            end
        end
        if !pass
            VTK.vtkexportmesh("mt017_s=$(surf)-p=$(p)-fails.vtk", connasarray(subset(fes, pel)), fens.xyz, VTK.T3)
        end
    end
    return spartitioning
end

