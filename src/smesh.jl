using StaticArrays
using MeshCore: nshapes, VecAttrib
using MeshCore: ir_skeleton, ir_transpose, ir_bbyfacets
using MeshCore: attribute
using MeshCore: P1, L2, T3, ShapeColl, manifdim, nfacets, facetdesc, nshapes
using MeshCore: IncRel

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
function orient_surface_mesh(t2v)
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
    return orientedt2v, orientable
end

function _vertfromedges(e, e2v)
    e1id = abs(e[1])
    e2id = abs(e[2])
    e3id = abs(e[3])
    _commonv(ea, eb) = (ea[1] == eb[1] ? ea[1] : (ea[2] == eb[1] ? ea[2] : (ea[2] == eb[2] ? ea[2] : ea[1])))
    if e[1] > 0
        return SVector(e2v[e1id][1], e2v[e1id][2], _commonv(e2v[e2id], e2v[e3id]))
    else
        return SVector(e2v[e1id][2], e2v[e1id][1], _commonv(e2v[e2id], e2v[e3id]))
    end
end

function _iscompatible(je, kes)
    for e in 1:3
        if je == -kes[e]
            return true
        end
    end
    return false
end

function _deduce_edge_order(je, kes)
    if _iscompatible(je, kes)
        return kes
    else 
        return SVector(-kes[3], -kes[2], -kes[1])
    end
end

function _orient_surface!(ot2ev, surfid, t2v, i, t2e, e2t)
    orientable = true
    currsurfid = maximum(surfid) + 1
    ot2ev[i] = t2e[i] # First Triangle: Accept its default orientation
    surfid[i] = currsurfid
    stack = Int[]
    push!(stack, i)
    while (!isempty(stack))
        j = pop!(stack)
        for e in 1:3
            eid = abs(ot2ev[j][e])
            if length(e2t[eid]) == 2 # is this a manifold edge?
                k = e2t[eid][1] == j ? e2t[eid][2] : e2t[eid][1]
                if surfid[k] == 0
                    ot2ev[k] = _deduce_edge_order(ot2ev[j][e], ot2ev[k]) 
                    surfid[k] = currsurfid
                    push!(stack, k)
                else
                    orientable = orientable && _iscompatible(ot2ev[j][e], ot2ev[k])
                end
            end
        end
    end
    return orientable
end
