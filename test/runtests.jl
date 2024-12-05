using Test
using ZipFile
using DataDrop

using DataDrop

function maybeunzip(impfun, file)
    if !isfile(file)
        zf = DataDrop.with_extension(file, ".zip")
        if !isfile(zf)
            error("File $(file) not found, and neither was its archive")
        end
        exdir=""
        fileFullPath = isabspath(zf) ?  zf : joinpath(pwd(),zf)
        basePath = dirname(fileFullPath)
        outPath = (exdir == "" ? basePath : (isabspath(exdir) ? exdir : joinpath(pwd(),exdir)))
        isdir(outPath) ? "" : mkdir(outPath)
        zarchive = ZipFile.Reader(fileFullPath)
        for f in zarchive.files
            fullFilePath = joinpath(outPath,f.name)
            if (endswith(f.name,"/") || endswith(f.name,"\\"))
                mkdir(fullFilePath)
            else
                write(fullFilePath, read(f))
            end
        end
        close(zarchive)
    end
    return impfun(file)
end


# module mt021_part
# using Random
# using FinEtools
# using FinEtools.MeshExportModule: VTK, MESH
# using FinEtools.MeshImportModule
# using Main: maybeunzip
# using ShellStructureTopo: make_topo_faces, create_partitions
# using MeshSteward: vtkwrite
# using Metis
# using Test
# function test()
#     output = maybeunzip(MeshImportModule.import_ABAQUS, joinpath("../models", "cylinders-93k.inp"))
#     fens, fes = output["fens"], output["fesets"][1]
#     surfids, partitionids, surface_elem_per_partition  = create_partitions(fens, fes, 2000; cluster_max_normal_deviation = pi/4)
#     VTK.vtkexportmesh("mt021_part.vtk", connasarray(fes), fens.xyz, VTK.T3; scalars=[("topological_face", surfids), ("partitioning", partitionids)]);
#     true
# end
# test()
# end

@testset "Randomization Tests" begin
    include("test_random.jl")
end

@testset "All Tests" begin
    include("test_all.jl")
end
