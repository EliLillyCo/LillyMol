using Test
here = @__DIR__ 
push!(LOAD_PATH, here)

using LillyMol
using LillyMolSubstructure

function is_failure(reason)
  println(reason)
  false
end

function foobar()
  m = Molecule()
  q = LillyMolSubstructure.SubstructureQuery()
end

function ok_smarts()::Bool
  q = LillyMolSubstructure.SubstructureQuery()
  q.build_from_smarts("CC") || return is_failure("Bad smarts")
  true
end

function main()
foobar()

@test ok_smarts()
end

main()
