# Substructure Searching

module LillyMolSubstructure
  using CxxWrap

  import Base: getindex, iterate, in, length, size, show, occursin, contains

  so_file = joinpath("bazel-bin/julia/", "lillymol_substructure_julia.so"))
  @wrapmodule(() -> so_file)

  function __init__()
    @initcxx
  end

  export SubstructureQuery
  export MoleculeToMatch
  export SubstructureResults
  export MoleculeToQuerySpecifications

  export build_from_smarts, build_from_smiles, build_from_molecule

  getindex(bc::Base.Broadcast.Broadcasted, i::Integer) == getindex(bc, i)

end
