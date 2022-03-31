
module mt007
using Random
using StaticArrays
using MeshCore: VecAttrib
using MeshCore: P1, L2, T3, ShapeColl, manifdim, nfacets, facetdesc, nshapes
using MeshCore: Q4ShapeDesc, shapedesc, n1storderv, nridges, nshifts, nvertices
using MeshCore: IncRel, attribute
using MeshCore: ir_skeleton, ir_transpose, ir_bbyfacets
using MeshSteward: T3block
using MeshSteward: vtkwrite
using MeshSteward: import_ABAQUS, vtkwrite, export_MESH, import_MESH
using ShellStructureTopo: orient_surface_mesh, identify_tes
using Test
function normal(c, geom)
    StaticArrays.cross(geom[c[2]] - geom[c[1]], geom[c[3]] - geom[c[1]])
end
function test()
    connectivities = import_ABAQUS(joinpath("../models", "extrusion.inp"))
    @test length(connectivities) == 1
    t2v = connectivities[1]
    for k in 1:nshapes(t2v.left)
        t2v._v[k] = shuffle!(Vector(t2v._v[k]))
    end
    e2v = ir_skeleton(t2v)
    # println(summary(e2v))
    t2e = ir_bbyfacets(t2v, e2v)
    orientedt2v, orientable = orient_surface_mesh(t2v)
    @test length(unique(attribute(orientedt2v.left, "surfid"))) == 3
    orientedt2v.right.attributes["geom"] = t2v.right.attributes["geom"]
    geom = orientedt2v.right.attributes["geom"]
    orientedt2v.right.attributes["normals"] = VecAttrib([normal(orientedt2v[i], geom) for i in 1:nshapes(orientedt2v.left)])
    vtkwrite("mt007-consistent", orientedt2v, [(name = "normals", )])
    # try rm("mt007-consistent.vtu"); catch end
    t2v = orientedt2v
    e2v = ir_skeleton(t2v)
    e2v.right.attributes["geom"] = t2v.right.attributes["geom"]
    # println(summary(e2v))
    t2e = ir_bbyfacets(t2v, e2v)
    e2t = ir_transpose(t2e)
    e2t = identify_tes(orientedt2v, e2t, 30*pi/180)
    isonte = attribute(e2t.left, "isonte")
    c = [e2v[idx] for idx in 1:nshapes(e2t.left) if isonte[idx]]
    lft = ShapeColl(shapedesc(e2v.left), length(c))
    tee2t = IncRel(lft, e2v.right, c)
    vtkwrite("mt007-te", tee2t)
    true
end
end
using .mt007
mt007.test()