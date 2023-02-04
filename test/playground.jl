
module mt022_part
using Random
using FinEtools
using FinEtools.MeshExportModule: VTK, MESH
using FinEtools.MeshImportModule
using ShellStructureTopo: make_topo_faces, create_partitions
using MeshSteward: vtkwrite
using Metis
using Test
function test()
    output = MeshImportModule.import_ABAQUS(joinpath("./models", "heat_105_5mm.inp"))
    fens, fes = output["fens"], output["fesets"][1]
    bfes = meshboundary(fes)
    surfids, partitionids, surface_elem_per_partition  = create_partitions(fens, bfes, 250; max_normal_deviation = pi/4)
    @show surface_elem_per_partition
    VTK.vtkexportmesh("mt022_part.vtk", connasarray(bfes), fens.xyz, VTK.T3; scalars=[("topological_face", surfids), ("partitioning", partitionids)]);
    true
end
@time test()
end

