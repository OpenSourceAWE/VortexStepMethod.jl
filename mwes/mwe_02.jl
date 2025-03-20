
# I tried to put a lazy cache (to avoid allocations for temporary vectors)
# into a struct. This does not work at all with the macro provided by the package 
# Parameters. It works (no syntax error) with the macro Base.@kwdef, but if you 
# make use of this cache it is extremely slow.
#
# This mwe can serve as starting point for further investigations. Currently the only
# cache that works is a cache in a global const variable.

using PreallocationTools, Parameters

Base.@kwdef mutable struct MockSolver
    const ca::NTuple{11, LazyBufferCache{typeof(identity), typeof(identity)}} = ([LazyBufferCache() for _ in 1:11])
end