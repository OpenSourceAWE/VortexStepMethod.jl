
using PreallocationTools, Parameters

Base.@kwdef mutable struct MockSolver
    const ca::NTuple{11, LazyBufferCache{typeof(identity), typeof(identity)}} = ([LazyBufferCache() for _ in 1:11])
end