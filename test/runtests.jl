using ReactionSensitivity
using Test

@testset "ReactionSensitivity.jl" begin
    
    if Sys.isapple() || Sys.islinux()
        lib_dir = "lib/"
    elseif Sys.iswindows()
        lib_dir = "lib\\"
    end
    input_file = "sensitivity.xml"
    retcode = rxn_gsa(input_file, lib_dir)
    @test retcode == Symbol("Success")

end
