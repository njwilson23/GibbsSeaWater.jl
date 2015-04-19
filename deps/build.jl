
function build()
    run(`make clean`)
    run(`make library`)
end
cd(build, "src/gsw_c_v3.03/")

if ~isdir("build")
    mkdir("build")
end

ext = @linux? (".so" : @osx? ( ".dylib" : @windows? ( ".dll" : "" )))
if ext == ""
    error("Platform not linux, OS X, or Windows")
end

run(`cp $("src/gsw_c_v3.03/libgswteos-10"*ext) $("build/libgswteos-10"*ext)`)

# using BinDeps
# 
# @BinDeps.setup
# 
# gsw = library_dependency("gsw")
# provides(Sources,
#          URI("http://www.teos-10.org/software/gsw_c_v3.03.zip"),
#          gsw,
#          unpacked_dir="gsw")
# 
# @build_steps begin
#     GetSources(gsw)
#     CreateDirectory("build")
#     FileUnpacker("gsw_c_v3.03.zip", "build")
#     @build_steps begin
#         ChangeDirectory("build")
#         MakeTargets("gsw_c_v3.03")
#         FileRule(joinpath("build", "gsw-3.03.so"), @build_steps begin
#             `cp gsw_c_v3.03/gsw-3.03.so .`
#         end)
#     end
# end
# 
# @BinDeps.install

