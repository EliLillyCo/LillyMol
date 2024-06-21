module BagModule
  using CxxWrap

  # abstract type AbstractSetOfAtoms <: AbstractVector{Int32} end

  @wrapmodule(joinpath("bazel-bin/julia/","bag_julia.so"))

  function __init__()
    @initcxx
  end

  import Base: getindex, iterate, in, length, size

  iterate(m::Bag, state=0) = (state >= (length(m) - 1) ? nothing : (m[state], state + 1))
  export Bag
end

using .BagModule
using Test

function test_bag()::Bool
  bag = BagModule.Bag(10)
  println("Created bag $(bag)")
  for i in 1:10
    println("i $(i)")
    bag[i - 1] = i + 100
  end
  return true
end

@test test_bag()
