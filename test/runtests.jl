using Test
using ZipFile

function unzip(file, exdir="")
    fileFullPath = isabspath(file) ?  file : joinpath(pwd(),file)
    basePath = dirname(fileFullPath)
    outPath = (exdir == "" ? basePath : (isabspath(exdir) ? exdir : joinpath(pwd(), exdir)))
    isdir(outPath) ? "" : mkdir(outPath)
    zarchive = ZipFile.Reader(fileFullPath)
    for f in zarchive.files
        fullFilePath = joinpath(outPath,f.name)
        if (endswith(f.name, "/") || endswith(f.name, "\\"))
            mkdir(fullFilePath)
        else
            write(fullFilePath, read(f))
        end
    end
    close(zarchive)
end

zips = [
"cylinders-24k.zip",  
"extrusion.zip",      
"stepped-cylinder.zip",                         
"cylinders-93k.zip",  
"quarter-shell.zip",  
]

for za in zips
    unzip(joinpath("../models", za))
end

@time @testset "All Tests" begin
    include("test_all.jl")
end
