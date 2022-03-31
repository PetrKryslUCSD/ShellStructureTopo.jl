
module mt008
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
using ShellStructureTopo: orient_surface_mesh, _normal, make_topo_faces
using Test
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
    t2v, orientable = orient_surface_mesh(t2v)
    @test orientable == true
    @test length(unique(attribute(t2v.left, "surfid"))) == 3
    geom = t2v.right.attributes["geom"]
    t2v.right.attributes["normals"] = VecAttrib([_normal(t2v[i], geom) for i in 1:nshapes(t2v.left)])
    vtkwrite("mt008-normals", t2v, [(name = "normals", )])
    # try rm("mt008-consistent.vtu"); catch end
    t2v = make_topo_faces(t2v)
    vtkwrite("mt008-classified", t2v, [(name = "tfid", )])
    true
end
end
using .mt008
mt008.test()