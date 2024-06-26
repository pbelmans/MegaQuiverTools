# Tutorial

## Installation

At the moment the package is not registered,
so you can install it by running the following command in the Julia REPL:

```julia-repl
pkg> add https://github.com/Catullo99/MegaQuiverTools.git
```

## Examples

To start using MegaQuiverTools in the REPL, one first must import it.

```julia-repl
julia> using MegaQuiverTools
```

Quivers can be built by passing the adjacency matrix to the Quiver() constructor:

```julia-repl
julia> Quiver([0 3; 0 0])
Quiver([0 3; 0 0], "")
```

The constructor accepts an optional string for naming the quiver:

```julia-repl
julia> MyQ = Quiver([0 3; 0 0], "My personal quiver")
Quiver([0 3; 0 0], "My personal quiver")
```

MegaQuiverTools has several constructors in place for many common examples:

```julia-repl
julia> mKroneckerquiver(4)
Quiver([0 4; 0 0], "4-Kronecker quiver")
julia> loopquiver(5)
Quiver([5;;], "5-loop quiver")
julia> subspacequiver(3)
Quiver([0 0 0 1; 0 0 0 1; 0 0 0 1; 0 0 0 0], "3-subspace quiver")
julia> threevertexquiver(1,6,7)
Quiver([0 1 6; 0 0 7; 0 0 0], "An acyclic 3-vertex quiver")
```


Dimension vectors and stability parameters are represented by `Vector{Int}` objects:

```julia-repl
julia> Q = mKroneckerquiver(3); d = [2,3];
julia> θ = canonical_stability(Q,d)
2-element Vector{Int64}:
  9
 -6
julia> iscoprime(d,θ)
true
```
Here, `iscoprime()` checks if ``d`` is θ-coprime, i.e., if any
proper subdimension vector ``0 \neq d' \nleq d`` satisfies ``\theta \cdot d' = 0``.

The bilinear Euler form relative to a quiver Q of any two vectors
in ``\mathbb{Z}Q`` can be computed:

```julia-repl
julia> Q = mKroneckerquiver(3); d = [2,2]; e = [3,4];
julia> Eulerform(Q,d,e)
-10
julia> Eulerform(Q,e,d)
-4
```

One can check if semistable, respectively stable representations
exist for a given dimension vector and stability parameter:

```julia-repl
julia> Q = mKroneckerquiver(3);
julia> d = [2,3];
julia> θ = [3,-2];
julia> hassemistables(Q,d,θ)
true
julia> hasstables(Q,d,θ)
true

julia> K2 = mKroneckerquiver(2)
Quiver([0 2; 0 0], "2-Kronecker quiver")
julia> hasstables(K2,[2,2],[1,-1])
false
julia> hassemistables(K2,[2,2],[1,-1])
true
```

One can also determine whether stable representations exist at all
for a given dimension vector by checking if it is a Schur root:

```julia-repl
julia> Q = mKroneckerquiver(3); d = [2,2];
julia> MegaQuiverTools.isSchurroot(Q,d)
true
julia> K2 = mKroneckerquiver(2);
julia> MegaQuiverTools.isSchurroot(K2,d)
false
```

To investigate the Harder-Narasimhan stratification of the parameter space
``\mathrm{R}(Q,\mathbf{d})``, the module provides a recursive closed formula.

```julia-repl
julia> Q = mKroneckerquiver(3); d = [2,3]; θ = [3,-2];
julia> allHNtypes(Q,d,θ)
8-element Vector{Vector{Vector{Int64}}}:
 [[2, 3]]
 [[1, 1], [1, 2]]
 [[2, 2], [0, 1]]
 [[2, 1], [0, 2]]
 [[1, 0], [1, 3]]
 [[1, 0], [1, 2], [0, 1]]
 [[1, 0], [1, 1], [0, 2]]
 [[2, 0], [0, 3]]
julia> isamplystable(Q,d,θ)
true
```

The method `isampystable()` determines whether the codimension of the θ-semistable locus,
``\mathrm{R}^{\theta-sst}(Q,\mathbf{d})\subset\mathrm{R}(Q,\mathbf{d})``, is at least 2.

The method `allHNtypes()` provides a list of all the Harder-Narasimhan types that appear in the problem.

The method `allTelemanbounds()` computes the bounds to apply Teleman quantization on the non-dense strata.
The output is a dictionary whose keys are the HN types and whose values are the weights themselves.

```julia-repl
julia> Q = mKroneckerquiver(3); d = [2,3]; theta = [3,-2];
julia> allTelemanbounds(Q,d,theta)
Dict{Vector{Vector{Int64}}, Int64} with 7 entries:
  [[2, 2], [0, 1]]         => 20
  [[2, 1], [0, 2]]         => 100
  [[1, 0], [1, 2], [0, 1]] => 100
  [[1, 0], [1, 3]]         => 120
  [[1, 0], [1, 1], [0, 2]] => 90
  [[1, 1], [1, 2]]         => 15
  [[2, 0], [0, 3]]         => 90
```

### Use cases

The following are examples of use cases for MegaQuiverTools.jl

**Verify Teleman inequalities**

In the following example, for each ``i,j`` and on each Harder-Narasimhan stratum,
we compute the weight of ``\mathcal{U}_i^\vee \otimes \mathcal{U}_j`` relative to the
1-PS corresponding to the HN stratum. These are then compared to the Teleman bounds.

```julia-repl
julia> Q = mKroneckerquiver(3); d = [2,3]; theta = [3,-2];
julia> hn = allTelemanbounds(Q,d,theta)
Dict{Vector{Vector{Int64}}, Int64} with 7 entries:
  [[2, 2], [0, 1]]         => 20
  [[2, 1], [0, 2]]         => 100
  [[1, 0], [1, 2], [0, 1]] => 100
  [[1, 0], [1, 3]]         => 120
  [[1, 0], [1, 1], [0, 2]] => 90
  [[1, 1], [1, 2]]         => 15
  [[2, 0], [0, 3]]         => 90

julia> endom = MegaQuiverTools.all_weights_endomorphisms_universal_bundle(Q,d,theta)
Dict{Vector{Vector{Int64}}, Int64} with 7 entries:
  [[2, 2], [0, 1]]         => 5
  [[2, 1], [0, 2]]         => 10
  [[1, 0], [1, 2], [0, 1]] => 15
  [[1, 0], [1, 3]]         => 15
  [[1, 0], [1, 1], [0, 2]] => 10
  [[1, 1], [1, 2]]         => 5
  [[2, 0], [0, 3]]         => 5

julia> check = all(endom[key] < hn[key] for key in keys(hn))
true
```

The fact that all of these inequalities are satisfied allows to conclude that the higher cohomology of
``\mathcal{U}_i^\vee \otimes \mathcal{U}_j`` vanishes.


## Bundle library

MegaQuiverTools contains a basic implementation of Bundle objects.

These are meant to be used as containers to perform computations with weights
or linearizations.

### Functionality

One can create a Bundle object by passing a list of weights and the rank to the constructor:

```julia-repl
julia> Bundle([1, 2, 3], 3)
Bundle of rank 3 with weights: [1, 2, 3]
```

`Note: the rank is not computed automatically from the weights. This is a design choice, to force a sanity check on the user.
Should there be the necessity to bypass this, the expected rank is the length of the list of weights.`

Direct sum, tensor product, wedge and symmetric product are implemented:

```julia-repl
julia> U = Bundle([1, 2, 3], 3); V = Bundle([4, 5], 2);
julia> U ⊕ V
Bundle of rank 5 with weights: [1, 2, 3, 4, 5]

julia> U ⊗ V
Bundle of rank 6 with weights: [5, 6, 6, 7, 7, 8]

julia> det(U)
Bundle of rank 1 with weights: [6]

julia> wedge(U, 0)
Bundle of rank 1 with weights: [0]
julia> wedge(U, 1)
Bundle of rank 3 with weights: [1, 2, 3]
julia> wedge(U, 2)
Bundle of rank 3 with weights: [3, 4, 5]
julia> wedge(U, 3)
Bundle of rank 1 with weights: [6]
julia> wedge(U, 4)
Bundle of rank 0 with weights: Int64[]

julia> symm(U,0)
Bundle of rank 1 with weights: [0]
julia> symm(U,1)
Bundle of rank 3 with weights: [1, 2, 3]
julia> symm(U,2)
Bundle of rank 6 with weights: [2, 3, 4, 4, 5, 6]
julia> symm(U,3)
Bundle of rank 10 with weights: [3, 4, 5, 5, 6, 7, 6, 7, 8, 9]
julia> symm(U,4)
Bundle of rank 15 with weights: [4, 5, 6, 6, 7, 8, 7, 8, 9, 10, 8, 9, 10, 11, 12]
```

A basic implementation of box products is also available:
  
```julia-repl

julia> U = Bundle([1, 2, 3], 3); V = Bundle([4, 5], 2);
julia> a = U ⊠ V
Bundle of rank 6 with weights: [[1, 4], [1, 5], [2, 4], [2, 5], [3, 4], [3, 5]]
julia> b = V ⊠ V
Bundle of rank 4 with weights: [[4, 4], [4, 5], [5, 4], [5, 5]]



julia> a ⊗ b
Bundle of rank 24 with weights: [[5, 8], [5, 9], [6, 8], [6, 9], [5, 9], [5, 10], [6, 9], [6, 10], [6, 8], [6, 9], [7, 8], [7, 9], [6, 9], [6, 10], [7, 9], [7, 10], [7, 8], [7, 9], [8, 8], [8, 9], [7, 9], [7, 10], [8, 9], [8, 10]]

julia> a ⊕ b
Bundle of rank 10 with weights: [[1, 4], [1, 5], [2, 4], [2, 5], [3, 4], [3, 5], [4, 4], [4, 5], [5, 4], [5, 5]]

julia> wedge(b,2)
Bundle of rank 6 with weights: [[8, 9], [9, 8], [9, 9], [9, 9], [9, 10], [10, 9]]

julia> wedge(b,3)
Bundle of rank 4 with weights: [[13, 13], [13, 14], [14, 13], [14, 14]]

julia> wedge(b,5)
Bundle of rank 0 with weights: Vector{Int64}[]

julia> symm(a,2)
Bundle of rank 21 with weights: [[2, 8], [2, 9], [3, 8], [3, 9], [4, 8], [4, 9], [2, 10], [3, 9], [3, 10], [4, 9], [4, 10], [4, 8], [4, 9], [5, 8], [5, 9], [4, 10], [5, 9], [5, 10], [6, 8], [6, 9], [6, 10]]
```

One use case is creating Bundle objects containing weights
of linearisations with respect to 1-PSs and then verifying Teleman inequality for
objects built with complex combinations of tensor, box and exterior products.

Another use case is defining line bundles on a projective space and then computing
the effects of direct summands, tensor products and exterior products on the
linearization. One such example is the computation of the Eagon-Northcott complex
for P^n, done below.

```julia-repl
julia> Q = mKroneckerquiver(2);
julia> U = [Bundle([-1], 1),Bundle([0], 1)]
2-element Vector{Bundle}:
  Bundle of rank 1, with weights: [-1]
  Bundle of rank 1, with weights: [0]

julia> EagonNorthcottcomplex(Q,U)
1-element Vector{Bundle}:
  Bundle of rank 1, with weights: [[-1, -1]]


julia> Q = mKroneckerquiver(5);
julia> EagonNorthcottcomplex(Q,U)
4-element Vector{Bundle}:
  Bundle of rank 10, with weights: [[-1, -1], [-1, -1], [-1, -1], [-1, -1], [-1, -1], [-1, -1], [-1, -1], [-1, -1], [-1, -1], [-1, -1]]
  Bundle of rank 20, with weights: [[-1, -2], [-2, -1], [-1, -2], [-2, -1], [-1, -2], [-2, -1], [-1, -2], [-2, -1], [-1, -2], [-2, -1], [-1, -2], [-2, -1], [-1, -2], [-2, -1], [-1, -2], [-2, -1], [-1, -2], [-2, -1], [-1, -2], [-2, -1]]
  Bundle of rank 15, with weights: [[-1, -3], [-2, -2], [-3, -1], [-1, -3], [-2, -2], [-3, -1], [-1, -3], [-2, -2], [-3, -1], [-1, -3], [-2, -2], [-3, -1], [-1, -3], [-2, -2], [-3, -1]]
  Bundle of rank 4, with weights: [[-1, -4], [-2, -3], [-3, -2], [-4, -1]]
```

Each Bundle object contains a rank and a list of weights.

