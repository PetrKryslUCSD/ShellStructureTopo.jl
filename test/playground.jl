
module mt021_part
using Random
using FinEtools
using FinEtools.MeshExportModule: VTK, MESH
using FinEtools.MeshImportModule
using ShellStructureTopo: make_topo_faces, create_partitions
using MeshSteward: vtkwrite
using Metis
using Test
function test()
    output = MeshImportModule.import_ABAQUS(joinpath("./models", "cylinders-93k.inp"))
    fens, fes = output["fens"], output["fesets"][1]
    # bfes = meshboundary(fes)
    surfids, partitionids, surface_elem_per_partition  = create_partitions(fens, fes, 250; max_normal_deviation = pi/4)
    @show surface_elem_per_partition
    VTK.vtkexportmesh("mt021_part.vtk", connasarray(fes), fens.xyz, VTK.T3; scalars=[("topological_face", surfids), ("partitioning", partitionids)]);
    true
end
test()
end

