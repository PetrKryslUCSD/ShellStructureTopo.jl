using StaticArrays
using MeshCore: nshapes, VecAttrib
using MeshCore: ir_skeleton, ir_transpose, ir_bbyfacets
using MeshCore: attribute
using MeshCore: P1, L2, T3, ShapeColl, manifdim, nfacets, facetdesc, nshapes
using MeshCore: IncRel

function _normal(c, geom)
    n = StaticArrays.cross(geom[c[2]] - geom[c[1]], geom[c[3]] - geom[c[1]])
    n /= StaticArrays.norm(n)
    return n
end

function identify_tes(t2v, e2t, crease_ang)
    cada = cos(crease_ang);
    geom = attribute(t2v.right, "geom")
    surfid = attribute(t2v.left, "surfid")
    # e2v = ir_skeleton(t2v)
    # t2e = ir_bbyfacets(t2v, e2v)
    # ot2ev = [SVector(t2e[i]) for i in 1:nshapes(t2e.left)]
    # e2t = ir_transpose(t2e)
    isonte = fill(false, nshapes(e2t.left))
    for e in 1:nshapes(e2t.left)
        numt_at_e = length(e2t[e])
        @assert numt_at_e >= 1
        if numt_at_e == 2 # manifold edge
            # if the triangles belong to different surfaces, this is a TE
            if surfid[e2t[e][1]] != surfid[e2t[e][2]] 
                println("different surfaces, $(surfid[e2t[e][1]]) vs $(surfid[e2t[e][2]])")
                isonte[e] = true
            else
                # check if it is a crease
                n1 = _normal(t2v[e2t[e][1]], geom)
                n2 = _normal(t2v[e2t[e][2]], geom)
                if StaticArrays.dot(n1, n2) < cada
                    println("crease, $(StaticArrays.dot(n1, n2))")
                    isonte[e] = true
                end
            end
        elseif numt_at_e == 1 # sheet edge
            isonte[e] = true
        else # non manifold edge
            isonte[e] = true
        end
    end
    e2t.left.attributes["isonte"] = VecAttrib(isonte)
    return e2t
end
