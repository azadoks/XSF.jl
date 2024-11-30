using Test
using XCrySDenStructureFormat
XSF = XCrySDenStructureFormat

@testset "XCrySDenStructureFormat.jl" begin
    @testset "Just load and save some files" begin
        basedir   = joinpath(@__DIR__, "data")
        testfiles = filter!(endswith(".xsf"), readdir(basedir))
        for file in testfiles
            system = last(XSF.load_xsf(joinpath(basedir, file)))
            mktempdir() do d
                outfile = joinpath(d, file)
                XSF.save_xsf(outfile, system)
            end
        end
    end
end
