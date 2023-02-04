using Test
using ZipFile
using DataDrop

@time @testset "All Tests" begin
    include("test_all.jl")
end
