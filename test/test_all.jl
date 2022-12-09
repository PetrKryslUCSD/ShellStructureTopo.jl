module mt001
using StaticArrays
using MeshCore: VecAttrib
using MeshCore: P1, L2, T3, ShapeColl, manifdim, nfacets, facetdesc, nshapes
using MeshCore: Q4ShapeDesc, shapedesc, n1storderv, nridges, nshifts, nvertices
using MeshCore: IncRel
using MeshCore: ir_skeleton, ir_transpose, ir_bbyfacets
using Test
function test()
    a = VecAttrib([1, 2, 4]);
    # @test show(a) != ""
    c = [(1, 2, 3), (4, 1, 3)]
    cc = [SVector{nvertices(T3)}(c[idx]) for idx in 1:length(c)]
    elements = ShapeColl(T3, length(c))
    vrts = ShapeColl(P1, 4)
    t2v = IncRel(elements, vrts, cc)
    e2v = ir_skeleton(t2v)
    # println(summary(e2v))
    t2e = ir_bbyfacets(t2v, e2v)
    # println(summary(t2e))
    # @show e2v
    @test t2e[1] == [4, -2, 1] 
    # t2v = ShapeColl(T3, 6, "t3s")
    # vrts = ShapeColl(P1, 12, "vrts")
    # ir = IncRel(q4s, vrts, cc)
    # @test summary(ir.left) == "q4s = 6 x Q4"
    # @test summary(ir.right) == "vrts = 12 x P1"
    # @test summary(ir) == "(q4s, vrts): q4s = 6 x Q4, vrts = 12 x P1" 
    # @show Q4
    # @show q4s
    true
end
end
using .mt001
mt001.test()

module mt002
using StaticArrays
using MeshCore: VecAttrib
using MeshCore: P1, L2, T3, ShapeColl, manifdim, nfacets, facetdesc, nshapes
using MeshCore: Q4ShapeDesc, shapedesc, n1storderv, nridges, nshifts, nvertices
using MeshCore: IncRel, attribute
using MeshCore: ir_skeleton, ir_transpose, ir_bbyfacets
using ShellStructureTopo: orient_surface_mesh
using Test
function test()
    a = VecAttrib([1, 2, 4]);
    # @test show(a) != ""
    c = [(1, 2, 3), (4, 1, 3)]
    cc = [SVector{nvertices(T3)}(c[idx]) for idx in 1:length(c)]
    elements = ShapeColl(T3, length(c))
    vrts = ShapeColl(P1, 4)
    t2v = IncRel(elements, vrts, cc)
    t2v.right.attributes["geom"] = VecAttrib([SVector(0.0, 0.0) for _ in 1:4])
    e2v = ir_skeleton(t2v)
    # println(summary(e2v))
    t2e = ir_bbyfacets(t2v, e2v)
    orientedt2v, orientable = orient_surface_mesh(t2v)
    @test length(unique(attribute(orientedt2v.left, "surfid"))) == 1
    # @show orientedt2v, orientable
    true
end
end
using .mt002
mt002.test()

module mt003
using StaticArrays
using MeshCore: VecAttrib
using MeshCore: P1, L2, T3, ShapeColl, manifdim, nfacets, facetdesc, nshapes
using MeshCore: Q4ShapeDesc, shapedesc, n1storderv, nridges, nshifts, nvertices
using MeshCore: IncRel, attribute
using MeshCore: ir_skeleton, ir_transpose, ir_bbyfacets
using ShellStructureTopo: orient_surface_mesh
using Test
function test()
    c = [(1, 2, 3), (4, 3, 1)]
    cc = [SVector{nvertices(T3)}(c[idx]) for idx in 1:length(c)]
    elements = ShapeColl(T3, length(c))
    vrts = ShapeColl(P1, 4)
    t2v = IncRel(elements, vrts, cc)
    t2v.right.attributes["geom"] = VecAttrib([SVector(0.0, 0.0) for _ in 1:4])
    e2v = ir_skeleton(t2v)
    # println(summary(e2v))
    t2e = ir_bbyfacets(t2v, e2v)
    orientedt2v, orientable = orient_surface_mesh(t2v)
    @test length(unique(attribute(orientedt2v.left, "surfid"))) == 1
    # @show t2v._v, orientedt2v._v, orientable
    true
end
end
using .mt003
mt003.test()

module mt004
using StaticArrays
using MeshCore: VecAttrib
using MeshCore: P1, L2, T3, ShapeColl, manifdim, nfacets, facetdesc, nshapes
using MeshCore: Q4ShapeDesc, shapedesc, n1storderv, nridges, nshifts, nvertices
using MeshCore: IncRel, attribute
using MeshCore: ir_skeleton, ir_transpose, ir_bbyfacets
using MeshSteward: T3block
using MeshSteward: vtkwrite
using ShellStructureTopo: orient_surface_mesh
using Test
function test()
    t2v = T3block(1.0, 2.0, 3, 6)
    e2v = ir_skeleton(t2v)
    # println(summary(e2v))
    t2e = ir_bbyfacets(t2v, e2v)
    orientedt2v, orientable = orient_surface_mesh(t2v)
    @test length(unique(attribute(orientedt2v.left, "surfid"))) == 1
    # @show t2v._v, orientedt2v._v, orientable
    orientedt2v.right.attributes["geom"] = t2v.right.attributes["geom"]
    vtkwrite("mt004", orientedt2v)
    try rm("mt004.vtu"); catch end
    true
end
end
using .mt004
mt004.test()


module mt005
using StaticArrays
using MeshCore: VecAttrib
using MeshCore: P1, L2, T3, ShapeColl, manifdim, nfacets, facetdesc, nshapes
using MeshCore: Q4ShapeDesc, shapedesc, n1storderv, nridges, nshifts, nvertices
using MeshCore: IncRel, attribute
using MeshCore: ir_skeleton, ir_transpose, ir_bbyfacets
using MeshSteward: T3block
using MeshSteward: vtkwrite
using MeshSteward: import_ABAQUS, vtkwrite, export_MESH, import_MESH
using ShellStructureTopo: orient_surface_mesh
using Test
function test()
    connectivities = import_ABAQUS(joinpath("../models", "extrusion.inp"))
    @test length(connectivities) == 1
    t2v = connectivities[1]
    e2v = ir_skeleton(t2v)
    # println(summary(e2v))
    t2e = ir_bbyfacets(t2v, e2v)
    orientedt2v, orientable = orient_surface_mesh(t2v)
    @test length(unique(attribute(orientedt2v.left, "surfid"))) == 3
    orientedt2v.right.attributes["geom"] = t2v.right.attributes["geom"]
    vtkwrite("mt005", orientedt2v)
    try rm("mt005.vtu"); catch end
    true
end
end
using .mt005
mt005.test()


module mt006
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
using ShellStructureTopo: orient_surface_mesh
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
    t2v.right.attributes["normals"] = VecAttrib([normal(t2v[i], geom) for i in 1:nshapes(t2v.left)])
    vtkwrite("mt006-inconsistent", t2v, [(name = "normals", )])
    # try rm("mt006-inconsistent.vtu"); catch end
    orientedt2v.right.attributes["normals"] = VecAttrib([normal(orientedt2v[i], geom) for i in 1:nshapes(orientedt2v.left)])
    vtkwrite("mt006-consistent", orientedt2v, [(name = "normals", )])
    # try rm("mt006-consistent.vtu"); catch end
    true
end
end
using .mt006
mt006.test()

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
using ShellStructureTopo: orient_surface_mesh, _classify_mesh_edges
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
    e2t = _classify_mesh_edges(orientedt2v, e2t, 30*pi/180)
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
    @test length(unique(attribute(t2v.left, "tfid"))) == 4
    vtkwrite("mt008-classified", t2v, [(name = "tfid", )])
    true
end
end
using .mt008
mt008.test()


module mt009
using Random
using FinEtools
using FinEtools.MeshExportModule: VTK
using FinEtools.MeshImportModule
using ShellStructureTopo: to_core
using MeshSteward: vtkwrite
using Test
function test()
    output = MeshImportModule.import_ABAQUS(joinpath("../models", "extrusion.inp"))
    fens, fes = output["fens"], output["fesets"][1]
    t2v = to_core(fens, fes)
    vtkwrite("mt009-faces", t2v)
    true
end
test()
end

module mt010
using Random
using FinEtools
using FinEtools.MeshExportModule: VTK
using FinEtools.MeshImportModule
using ShellStructureTopo: to_core
using MeshSteward: vtkwrite
using Test
function test()
    output = MeshImportModule.import_ABAQUS(joinpath("../models", "quarter-shell.inp"))
    fens, fes = output["fens"], output["fesets"][1]
    t2v = to_core(fens, fes)
    vtkwrite("mt010-faces", t2v)
    true
end
test()
end

