module MegaQuiverTools

import Base.show

export Quiver
# are these useful?
export nvertices, narrows, indegree, outdegree, isacyclic, isconnected, issink, issource
export Eulerform, canonical_stability, iscoprime, slope
export isSchurroot
export allHNtypes, hassemistables, hasstables, isamplystable
export allTelemanbounds, all_weights_endomorphisms_universal_bundle, all_weights_universal_bundle, all_weights_irreducible_component_canonical
export mKroneckerquiver, loopquiver, subspacequiver, threevertexquiver

using Memoize, LinearAlgebra




"""
A quiver is represented by its adjacency
``n \\times n`` matrix ``adjacency = (a_{ij}),
``
where ``n`` is the number of vertices
and ``a_{ij}`` is the number of arrows i → j.

Attributes:

- `adjacency` is the adjacency matrix of the quiver
- `name` is the name of the quiver, defaults to `""`.
"""
mutable struct Quiver
    adjacency::Matrix{Int}
    name::String

    function Quiver(adjacency::Matrix{Int}, name::String)
        if !(size(adjacency)[1] == size(adjacency)[2])
            throw(DomainError(adjacency, "adjacency matrix must be square"))
        else 
            new(adjacency, name)
        end
    end 
    function Quiver(adjacency::Matrix{Int})
        if !(size(adjacency)[1] == size(adjacency)[2])
            throw(DomainError(adjacency, "adjacency matrix must be square"))
        else
            new(adjacency, "")
        end
    end
end

function show(io::IO, Q::Quiver)
    print(io, Q.name*", with adjacency matrix ")
    println(io, Q.adjacency)
end

# """
# Returns the adjacency matrix of the quiver.

# OUTPUT: A square matrix M whose entry M[i,j] is the number of arrows from the vertex i to the vertex j.
# """
# function adjacency_matrix(Q::Quiver)
#     return Q.adjacency
# end

"""
Returns the (necessarily symmetric) adjacency matrix of the underlying graph of the quiver.
"""
function underlying_graph(Q::Quiver)
    return Matrix{Int}(Q.adjacency + transpose(Q.adjacency) - diagm(diag(Q.adjacency)))
end

"""
Returns the number of vertices of the quiver.
"""
nvertices(Q::Quiver) = size(Q.adjacency)[1]

"""
Returns the number of arrows of the quiver.
"""
narrows(Q::Quiver) = sum(Q.adjacency)

"""
Checks wether the quiver is acyclic, i.e. has no oriented cycles.
"""
isacyclic(Q::Quiver) = all(entry == 0 for entry in Q.adjacency^nvertices(Q))

# TODO probably not needed
"""
Checks wether the underlying graph of the quiver is connected.

Examples:
```julia-repl
julia> Q = Quiver([0 1 0; 0 0 1; 1 0 0])
julia> isconnected(Q)
true

julia> Q = Quiver([0 1 0; 1 0 0; 0 0 2])
false

# The 4-Kronecker quiver:
julia> Q = mKroneckerquiver(4)
julia> isconnected(Q)
true

# The 4-loop quiver:
julia> Q = LoopQuiver(4)
julia> isconnected(Q)
true

# The 4-subspace quiver:
julia> Q = SubspaceQuiver(4)
julia> isconnected(Q)
true

# The A10 quiver:
julia> A10 = Quiver(   [0 1 0 0 0 0 0 0 0 0;
                        0 0 1 0 0 0 0 0 0 0;
                        0 0 0 1 0 0 0 0 0 0;
                        0 0 0 0 1 0 0 0 0 0;
                        0 0 0 0 0 1 0 0 0 0;
                        0 0 0 0 0 0 1 0 0 0;
                        0 0 0 0 0 0 0 1 0 0;
                        0 0 0 0 0 0 0 0 1 0;
                        0 0 0 0 0 0 0 0 0 1;
                        0 0 0 0 0 0 0 0 0 0] )
julia> isconnected(A10)
true

# The A10 quiver without one arrow:
julia> A10 = Quiver(   [0 1 0 0 0 0 0 0 0 0;
                        0 0 1 0 0 0 0 0 0 0;
                        0 0 0 1 0 0 0 0 0 0;
                        0 0 0 0 1 0 0 0 0 0;
                        0 0 0 0 0 1 0 0 0 0;
                        0 0 0 0 0 0 0 0 0 0;
                        0 0 0 0 0 0 0 1 0 0;
                        0 0 0 0 0 0 0 0 1 0;
                        0 0 0 0 0 0 0 0 0 1;
                        0 0 0 0 0 0 0 0 0 0] )
julia> isconnected(A10)
false
```
"""
function isconnected(Q::Quiver)
    paths = underlying_graph(Q)
    for i in 2:nvertices(Q) - 1
        paths += paths*underlying_graph(Q)
    end
    for i in 1:nvertices(Q), j in 1:nvertices(Q)
            if i != j && paths[i, j] == 0 && paths[j, i] == 0
                return false
            end
    end
    return true
end

"""
Returns the number of incoming arrows to the vertex j.

Examples:
```julia-repl
julia> Q = mKroneckerquiver(4)
julia> indegree(Q, 1)
0
julia> indegree(Q, 2)
4
```
"""
indegree(Q::Quiver, j::Int) = sum(Q.adjacency[:, j])

"""
Returns the number of outgoing arrows from the vertex i.

Examples:
```julia-repl
julia> Q = mKroneckerquiver(4)
julia> outdegree(Q, 1)
4
julia> outdegree(Q, 2)
0
```
"""
outdegree(Q::Quiver, i::Int) = sum(Q.adjacency[i, :])

"""
Checks if the vertex i is a source, i.e., a vertex with no incoming arrows.

Examples:
```julia-repl
julia> Q = mKroneckerquiver(4)
julia> issource(Q, 1)
true
julia> issource(Q, 2)
false
```
"""
issource(Q::Quiver, i::Int) = indegree(Q, i) == 0

"""
Checks if the vertex j is a sink, i.e., a vertex with no outgoing arrows.

Examples:
```julia-repl
julia> Q = mKroneckerquiver(4)
julia> issink(Q, 1)
false
julia> issink(Q, 2)
true
```
"""
issink(Q::Quiver, j::Int) = outdegree(Q, j) == 0

function arrows(Q::Quiver)
    return reduce(vcat, [[i,j] for k in 1:Q.adjacency[i,j]] for i in 1:nvertices(Q) for j in 1:nvertices(Q) if Q.adjacency[i,j] > 0)
end

"""
Returns the Euler matrix of the quiver.

The Euler matrix of a quiver Q is defined as 
```math
E = I - A,
```
where ``A`` is the adjacency matrix of Q and ``I`` is the identity matrix of the same size as ``A``.
"""
@memoize Dict euler_matrix(Q::Quiver) = Matrix{Int}(I, nvertices(Q), nvertices(Q)) - Q.adjacency

"""
Computes the Euler form of the quiver for vectors x and y.

The Euler form is defined as the bilinear form
```math
\\langle x,y\\rangle = x^T * E * y,
```
where E is the Euler matrix of the quiver.
"""
Eulerform(Q::Quiver, x::Vector{Int}, y::Vector{Int}) = x'*euler_matrix(Q)*y

"""
The canonical stability parameter for the couple ``(Q,d)`` is given by ``<d,- > - < - ,d>``
"""
canonical_stability(Q::Quiver, d::Vector{Int})::Vector{Int} = -(-transpose(euler_matrix(Q)) + euler_matrix(Q))*d

"""Checks wether the given dimension vector ``d`` is `theta`-coprime for the stability parameter `theta`."""
function iscoprime(d::Vector{Int}, theta::Vector{Int})
    return all(e -> theta'*e != 0, all_proper_subdimension_vectors(d))
end

"""Checks if the gcd of all the entries of d is ``1``."""
function iscoprime(d::Vector{Int})
    return gcd(d) == 1
end

function slope(d::Vector{Int}, theta::Vector{Int}, slope_denominator::Function = sum)
    return (theta'*d)//slope_denominator(d)
end

@memoize Dict function all_forbidden_subdimension_vectors(d::Vector{Int}, theta::Vector{Int}, slope_denominator::Function = sum)
    return filter(e -> slope(e, theta, slope_denominator) > slope(d, theta, slope_denominator), all_proper_subdimension_vectors(d))
end


"""Checks if there is a theta-semistable representation of dimension vector d.

Examples:
```julia-repl
julia> A2 = mKroneckerquiver(1)
julia> theta = [1,-1]
julia> d = [1,1]
julia> hassemistables(A2, d, theta)
true

julia> d = [2,2]
julia> hassemistables(A2, d, theta)
true

julia> d = [1,2]
julia> hassemistables(A2, d, theta)
false

julia> d = [0,0]
julia> hassemistables(A2, d, theta)
true

# The 3-Kronecker quiver:
julia> K3 = mKroneckerquiver(3)
julia> theta = [3,-2]
julia> d = [2,3]
julia> hassemistables(K3, d, theta)
true

julia> d = [1,4]
julia> hassemistables(K3, d, theta)
false
```
"""
@memoize Dict function hassemistables(Q::Quiver, d::Vector{Int}, theta::Vector{Int}, slope_denominator::Function = sum; algorithm="schofield")
    if all(di == 0 for di in d)
        return true
    else
        # collect the list of all subdimension vectors e of bigger slope than d
        slope_d = slope(d, theta, slope_denominator)
        subdimensionsBiggerSlope = filter(e -> slope(e, theta, slope_denominator) > slope_d, all_proper_subdimension_vectors(d))
        # to have semistable representations, none of the vectors above must be generic subdimension vectors.
        return all(e -> !is_generic_subdimension_vector(Q, e, d, algorithm="schofield"), subdimensionsBiggerSlope)
    end
end

# TODO write examples with properly semistable representations and with sst nonempty and st empty.
# TODO put these in unit tests
"""Checks if Q has a theta-stable representation of dimension vector d."""
@memoize Dict function hasstables(Q::Quiver, d::Vector{Int}, theta::Vector{Int}; slope_denominator::Function = sum)
    if all(di == 0 for di in d)
        return true
    else
        # collect the list of all subdimension vectors e of bigger slope than d
        slope_d = slope(d, theta, slope_denominator)
        subdimensions_bigger_or_equal_slope = filter(e -> slope(e, theta, slope_denominator) >= slope_d, all_proper_subdimension_vectors(d))
        # to have semistable representations, none of the vectors above must be generic subdimension vectors.
        return all(e -> !is_generic_subdimension_vector(Q, e, d), subdimensions_bigger_or_equal_slope)
    end
end

"""
Checks if d is a Schur root for Q.
By a lemma of Schofield (See Lemma 4.2 of https://arxiv.org/pdf/0802.2147.pdf),
this is equivalent to the existence of a stable representation of dimension vector d."""
isSchurroot(Q::Quiver, d::Vector{Int}) = hasstables(Q, d, canonical_stability(Q, d))

"""Checks if e is a generic subdimension vector of d.

A dimension vector e is called a generic subdimension vector of d if a generic representation
of dimension vector d possesses a subrepresentation of dimension vector e.

By a result of Schofield (see Thm. 5.3 of https://arxiv.org/pdf/0802.2147.pdf)
e is a generic subdimension vector of d if and only if
``<e',d-e> \\geq 0``
for all generic subdimension vectors e' of e.
"""
@memoize Dict function is_generic_subdimension_vector(Q::Quiver, e::Vector{Int}, d::Vector{Int}; algorithm::String = "schofield")
    if e == d
        return true
    elseif all(ei==0 for ei in e)
        return false
    else
    # considering subdimension vectors that violate the numerical condition
    # TODO this filtering is inefficent. For fixed d-e, this is a LINEAR form, we KNOW which eprimes violate the condition. We should just check those.
        euler_matrix_temp = euler_matrix(Q) * (d-e) #to speed up computation of <eprime,d-e>
        subdimensions = filter(eprime -> eprime'*euler_matrix_temp < 0, all_nonzero_subdimension_vectors(e))
        # none of the subdimension vectors violating the condition should be generic
        return !any(eprime -> is_generic_subdimension_vector(Q, eprime, e, algorithm="schofield"), subdimensions)
    end
end

"""
Returns the list of all generic subdimension vectors of d.
"""
function all_generic_subdimension_vectors(Q::Quiver, d::Vector{Int}) 
    return filter(e -> is_generic_subdimension_vector(Q, e, d), all_subdimension_vectors(d))
end


"""
Returns a list of all the Harder Narasimhan types of representations of Q with dimension vector d, with respect to the slope function theta/slope_denominator.
"""
@memoize Dict function allHNtypes(Q::Quiver, d::Vector{Int}, theta::Vector{Int}, slope_denominator::Function=sum; ordered=true)

    if all(di == 0 for di in d)
        return [[d]]
    else
        # We consider just proper subdimension vectors which admit a semistable representation and for which μ(e) > μ(d)
        # Note that we also eliminate d by the following
        subdimensions = filter(e -> hassemistables(Q,  e, theta, slope_denominator), all_forbidden_subdimension_vectors(d, theta, slope_denominator)) 
        
        # We sort the subdimension vectors by slope because that will return the list of all HN types in ascending order with respect to the partial order from Def. 3.6 of https://mathscinet.ams.org/mathscinet-getitem?mr=1974891
        if ordered
            subdimensions = sort(subdimensions, by = e -> slope(e, theta, slope_denominator))
        end

        # The HN types which are not of the form (d) are (e,f^1,...,f^s) where e is a proper semistable subdimension vector with μ(e) > μ(d), (f^1,...,f^s) is a HN type of f = d-e and μ(e) > μ(f^1) holds.

        alltypes = [[e, efstar...] for e in subdimensions for efstar in filter(fstar -> slope(e, theta, slope_denominator) > slope(fstar[1], theta, slope_denominator), allHNtypes(Q, d-e, theta, slope_denominator, ordered=ordered))]

        # Possibly add d again, at the beginning, because it is smallest with respect to the partial order from Def. 3.6
        if hassemistables(Q, d, theta, slope_denominator)
            pushfirst!(alltypes, [d])
        end
        return alltypes
    end
end

"""
Returns the codimension of the given HN stratum.
"""
function codimensionHNstratum(Q::Quiver, stratum::Vector{Vector{Int}})
    if length(stratum) == 1
        return 0
    else
        return -sum(Eulerform(Q, stratum[i], stratum[j]) for i in 1:length(stratum)-1 for j in i+1:length(stratum))
    end
end

"""
Checks wether the dimension vector d is amply stable with respect to the slope function theta/denominator.
    
This means that the codimension of the unstable locus in the parameter space is at least 2.
"""
function isamplystable(Q::Quiver, d::Vector{Int}, theta::Vector{Int}, slope_denominator::Function = sum)
    # We say that representations of a given dimension vector d are amply stable (for any notion of stability) if the codimension of the semistable locus is at least 2.
    # We verify this by computing the codimension of each HN stratum.
    HN = filter(hntype -> hntype != [d], allHNtypes(Q, d, theta, slope_denominator))
    return all(stratum -> codimensionHNstratum(Q, stratum) >= 2, HN)
end

######################################################################
# Weights of various standard vector bundles for the HN stratification
######################################################################


""" Computes the weight on ``\\det(N_{S/R}|_Z)`` of the 1-PS ``\\lambda``
corresponding to the given HN type."""
function Telemanbound_onstratum(Q::Quiver, hntype::Vector{Vector{Int}}, theta::Vector{Int}, slope_denominator::Function = sum)::Int
    if length(hntype) == 1
        throw(ArgumentError("Weight not defined for HN type of length 1."))
    end
    slopes = map(h -> slope(h, theta, slope_denominator), hntype)
    slopes = lcm(denominator.(slopes)) .* slopes
    return sum((slopes[t] - slopes[s])*Eulerform(Q, hntype[s], hntype[t]) for s in 1:length(hntype)-1 for t in s+1:length(hntype))
end

""" Computes the weight on ``\\det(N_{S/R}|_Z)`` of the 1-PS corresponding to each
HN type for the given Q, d, theta and denominator."""
function allTelemanbounds(Q::Quiver, d::Vector{Int}, theta::Vector{Int}, slope_denominator::Function = sum)
    #This is only relevant on the unstable locus
    HN = filter(hntype -> hntype != [d], allHNtypes(Q, d, theta, slope_denominator))
    return Dict([hntype,Telemanbound_onstratum(Q, hntype, theta, slope_denominator)] for hntype in HN)
end

"""Returns the weights of a universal bundle ``U_i(a)`` for the linearization ``a``
for the 1-PS corresponding to the given HN type."""
function weights_universal_bundle_onstratum(theta::Vector{Int}, a::Vector{Int}, hntype, slope_denominator::Function = sum)::Vector{Int}
    slopes = map(h -> slope(h, theta, slope_denominator), hntype)
    slopes *= lcm(denominator.(slopes))

    constant_term = sum(slopes[i]* (a' * hntype[i]) for i in eachindex(hntype))

    return -constant_term .+ slopes
end
"""Computes the weights of the universal bundle ``U_i(a)`` for the linearization ``a``
on all the non-dense Harder-Narasimhan strata for each 1-PS corresponding to each HN type."""
function all_weights_universal_bundle(Q::Quiver, d::Vector{Int}, theta::Vector{Int}, a::Vector{Int}, slope_denominator::Function = sum)
    HN = filter(hntype -> hntype != [d], allHNtypes(Q, d, theta, slope_denominator))
    return Dict([hntype, weights_universal_bundle_onstratum(theta, a, hntype, slope_denominator)] for hntype in HN)
end


"""Computes the weight of the irreducible component of ``\\omega_R|_Z``
on a Harder-Narasimhan stratum for the 1-PS corresponding to each HN type.
More explicitly, if ``\\omega_X = O(rH)``, this returns the weight of
the pullback of O(H) on the given stratum."""
function weight_irreducible_component_canonical_on_stratum(Q::Quiver,d::Vector{Int},hntype::Vector{Vector{Int}},theta::Vector{Int},slope_denominator::Function = sum)::Int
    kweights = map(di -> slope(di,theta, slope_denominator), hntype)
    kweights = kweights * lcm(denominator.(kweights))

    dd = sum( kweights[m] .* hntype[m] for m in 1:length(hntype))
    # The Fano paper shows that under appropriate conditions,
    # the canonical bundle is given by linearizing with minus
    # the canonical stability parameter.
    can = canonical_stability(Q,d)
    can /= gcd(can)
    return can' * dd
end

"""Computes the weights of the irreducible component of ``\\omega_R|_Z``
on all the non-dense Harder-Narasimhan strata for each 1-PS relative to the HN type.
More explicitly, if ``\\omega_X = O(rH)``, this returns the weights of the pullback of O(H) on each stratum."""
function all_weights_irreducible_component_canonical(Q::Quiver,d::Vector{Int},theta::Vector{Int}, slope_denominator::Function = sum)
    HN = filter(hntype -> hntype != [d], allHNtypes(Q,d,theta))
    return Dict([hntype, weight_irreducible_component_canonical_on_stratum(Q, d, hntype, theta, slope_denominator)] for hntype in HN)
end

"""Computes the weights of the endomorphism of the universal bundle ``U_i \\otimes U_j``
on the given Harder-Narasimhan stratum for the 1-PS relative to the HN type."""
function weights_endomorphism_universal_bundle_on_stratum(hntype::Vector{Vector{Int}},theta::Vector{Int},slope_denominator::Function = sum)::Vector{Int}
    # the maximum weight of the tensors of the universal bundles U_i^\vee \otimes U_j is slope of first term in the HN type - slope of the last term in the HN type
    kweights = map(di -> slope(di,theta, slope_denominator), hntype)
    kweights = kweights * lcm(denominator.(kweights))
    # return kweights[1] - kweights[end] # this is the largest one
    return [kweights[i] - kweights[j] for i in 1:length(hntype) for j in 1:length(hntype)]
end

"""Computes the weights of the endomorphisms of the universal bundles ``U_i \\otimes U_j``
on all the non-dense Harder-Narasimhan strata for each 1-PS relative to the HN type."""
function all_weights_endomorphisms_universal_bundle(Q::Quiver,d::Vector{Int},theta::Vector{Int}, slope_denominator::Function = sum)
    HN = filter(hntype -> hntype != [d], allHNtypes(Q, d, theta, slope_denominator))
    return Dict([hntype, weights_endomorphism_universal_bundle_on_stratum(hntype, theta, slope_denominator)] for hntype in HN)
end


#################
# Constructors
#################

function mKroneckerquiver(m::Int)
    return Quiver([0 m; 0 0], string(m)*"-Kronecker quiver")
end

function threevertexquiver(m12::Int, m13::Int, m23::Int)
    return Quiver([0 m12 m13; 0 0 m23; 0 0 0], "An acyclic 3-vertex quiver")
end

function loopquiver(m::Int)
    return Quiver(Matrix{Int}(reshape([m], 1, 1)), string(m)*"-loop quiver")
end

function subspacequiver(m::Int)
    A = zeros(Int, m+1, m+1)
    for i in 1:m
        A[i, m+1] = 1
    end
    return Quiver(A, string(m)*"-subspace quiver")
end


 

#################
# Technical tools
#################

function thin_dimension_vectors(Q::Quiver)
    return ones(Int, nvertices(Q))
end

@memoize Dict function all_subdimension_vectors(d::Vector{Int})
    return collect.(Iterators.product(map(di -> 0:di, d)...))
end
@memoize Dict function all_nonzero_subdimension_vectors(d::Vector{Int})::Vector{Vector{Int}}
    return filter(e->!all(ei == 0 for ei in e), all_subdimension_vectors(d))
end

@memoize Dict function all_proper_subdimension_vectors(d::Vector{Int})::Vector{Vector{Int}}
    return filter(e -> any(ei != 0 for ei in e) && e != d, all_subdimension_vectors(d))
end

include("EN-weights.jl")

end