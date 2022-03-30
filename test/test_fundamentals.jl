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
    println(summary(e2v))
    t2e = ir_bbyfacets(t2v, e2v)
    println(summary(t2e))
    @show e2v
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
using MeshCore: IncRel
using MeshCore: ir_skeleton, ir_transpose, ir_bbyfacets
using ShellStructureTopo: smesh_orient
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
    println(summary(e2v))
    t2e = ir_bbyfacets(t2v, e2v)
    orientedt2v, orientable = smesh_orient(t2v)
    @show orientedt2v, orientable
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
using MeshCore: IncRel
using MeshCore: ir_skeleton, ir_transpose, ir_bbyfacets
using ShellStructureTopo: smesh_orient
using Test
function test()
    c = [(1, 2, 3), (4, 3, 1)]
    cc = [SVector{nvertices(T3)}(c[idx]) for idx in 1:length(c)]
    elements = ShapeColl(T3, length(c))
    vrts = ShapeColl(P1, 4)
    t2v = IncRel(elements, vrts, cc)
    e2v = ir_skeleton(t2v)
    println(summary(e2v))
    t2e = ir_bbyfacets(t2v, e2v)
    orientedt2v, orientable = smesh_orient(t2v)
    @show t2v._v, orientedt2v._v, orientable
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
using ShellStructureTopo: smesh_orient
using Test
function test()
    t2v = T3block(1.0, 2.0, 3, 6)
    e2v = ir_skeleton(t2v)
    println(summary(e2v))
    t2e = ir_bbyfacets(t2v, e2v)
    orientedt2v, orientable = smesh_orient(t2v)
    @show t2v._v, orientedt2v._v, orientable
    orientedt2v.right.attributes["geom"] = t2v.right.attributes["geom"]
    vtkwrite("mt004", orientedt2v)
    # try rm("mt004.vtu"); catch end
    true
end
end
using .mt004
mt004.test()