
module mt019_part
using Random
using FinEtools
using FinEtools.MeshExportModule: VTK, MESH
using FinEtools.MeshImportModule
using ShellStructureTopo: make_topo_faces, create_partitions
using MeshSteward: vtkwrite
using Metis
using Test
function test()
    output = MeshImportModule.import_ABAQUS(joinpath("../models", "cylinders-93k.inp"))
    fens, fes = output["fens"], output["fesets"][1]
    surfids, partitionids, surface_elem_per_partition  = create_partitions(fens, fes, 2000; max_normal_deviation = pi/6)
    @show surface_elem_per_partition
    VTK.vtkexportmesh("mt019_part.vtk", connasarray(fes), fens.xyz, VTK.T3; scalars=[("topological_face", surfids), ("partitioning", partitionids)]);
    true
end
test()
end

