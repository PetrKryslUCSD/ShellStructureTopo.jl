using FinEtools


"""
    orient_surface_mesh(t2v)

Orient surface mesh.

Return
- `orientedt2v`: oriented incidence relation
- `orientable`: Can the surface mesh be oriented? Bool flag.
    This would be false for a Mobius strip, for instance.

This is a purely topological operation. Creases in the surface are not
recognized as edges. The operation proceeds by flooding the surface across
manifold edges, stopping at sheet or non-manifold edges.

The incidence relation `orientedt2v` includes an attribute for each triangle
that provides the surface id to which the triangle belongs
(`orientedt2v.left.attributes["surfid"]`).
"""
function orient_surface_mesh(fens, fes)
    orientable = true
    all_oriented = false
    surfid = fill(0, nshapes(t2v.left))
    e2v = ir_skeleton(t2v)
    t2e = ir_bbyfacets(t2v, e2v)
    ot2ev = [SVector(t2e[i]) for i in 1:nshapes(t2e.left)]
    e2t = ir_transpose(t2e)
    while !all_oriented
        all_oriented = true
        for i in 1:nshapes(t2v.left)
            if surfid[i] == 0
                orientable = orientable && _orient_surface!(ot2ev, surfid, t2v, i, t2e, e2t);
                all_oriented = false;
                break;
            end
        end
    end
    elements = ShapeColl(T3, nshapes(t2v.left))
    vrts = ShapeColl(P1, nshapes(t2v.right))
    oc = [_vertfromedges(ot2ev[i], e2v) for i in 1:length(ot2ev)]
    orientedt2v = IncRel(elements, vrts, oc)
    orientedt2v.left.attributes["surfid"] = VecAttrib(surfid)
    orientedt2v.right.attributes["geom"] = attribute(t2v.right, "geom")
    return orientedt2v, orientable
end

"""
    make_topo_faces(t2v, e2t)

Make topological faces.

Returns
- `t2v`: the attribute `tfid` lists the numbers of the topological faces.
"""
function make_topo_faces(fens, fes, crease_ang = 30/180*pi)
    e2v = ir_skeleton(t2v)
    e2v.right.attributes["geom"] = t2v.right.attributes["geom"]
    t2e = ir_bbyfacets(t2v, e2v)
    e2t = ir_transpose(t2e)
    e2t = _classify_mesh_edges(t2v, e2t, crease_ang)
    all_done = false
    tfid = fill(0, nshapes(t2v.left))
    while !all_done
        all_done = true
        for i in 1:nshapes(t2v.left)
            if tfid[i] == 0
                _collect_tf!(tfid, t2v, i, t2e, e2t)
                all_done = false;
                break;
            end
        end
    end
    t2v.left.attributes["tfid"] = VecAttrib(tfid)
    return t2v
end

function to_core(fens, fes)
    isa(fes, FinEtools.FESetModule.FESetT3) || error("Only 3-node triangles accepted")
    v = ShapeColl(P1, count(fens), "fens")
    N, T = size(fens.xyz, 2), eltype(fens.xyz)
    v.attributes["geom"] = VecAttrib([SVector{N, T}(fens.xyz[i, :]) for i in 1:size(fens.xyz, 1)])
    t = ShapeColl(T3, count(fes), "fes")
    IncRel(t, v, connasarray(fes))
end
