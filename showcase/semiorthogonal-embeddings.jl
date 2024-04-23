# This script checks how many copies of D^b(Q) via the "Fourier-Mukai-like" functor in D^b(M) are semiorthogonal.
# It shows which Teleman vanishings can be applied and which Euler characteristics are zero.

using MegaQuiverTools, IntersectionTheory

function eulerform_UidualUjtimesH(Q::Quiver,d::Vector{Int},twist::Int = 1)

    if gcd(d) != 1
        throw(ArgumentError("gcd(d,e) != 1; universal bundles don't exist"))
    end
    if !hassemistables(Q, d, canonical_stability(Q,d))
        throw(ArgumentError("Moduli space is empty"))
    end

    r = gcd(canonical_stability(Q,d))
    M = IntersectionTheory.matrix_moduli(narrows(Q),d[1],d[2])

    # coefficients for the canonical bundle
    c1 = Eulerform(Q,d,[1,0]) - Eulerform(Q,[1,0],d)
    c2 = Eulerform(Q,d,[0,1]) - Eulerform(Q,[0,1],d)

    #e1 is c_1(\mathcal{U}_2), f1 is c_1(\mathcal{U}_1)
    e1,f1 = gens(M.ring.R)[1], gens(M.ring.R)[4]
    h = - twist//r * (c2*M.ring(e1) + c1*M.ring(f1))
    H = OO(M,h)

    return all(chi(u*IntersectionTheory.dual(v)*H) == 0 for u in M.bundles for v in M.bundles)
end

function eulerform_UitimesH(Q::Quiver, d::Vector{Int}, twist::Int = 1)

    if gcd(d) != 1
        throw(ArgumentError("gcd(d,e) != 1; universal bundles don't exist"))
    end
    if !hassemistables(Q, d, canonical_stability(Q,d))
        throw(ArgumentError("Moduli space is empty"))
    end

    r = gcd(canonical_stability(Q,d))
    M = IntersectionTheory.matrix_moduli(narrows(Q),d[1],d[2])

    if twist == 0
        return all(chi(u) == 0 for u in M.bundles)
    else
    c1 = Eulerform(Q,d,[1,0]) - Eulerform(Q,[1,0],d)
        c2 = Eulerform(Q,d,[0,1]) - Eulerform(Q,[0,1],d)

        e1,f1 = gens(M.ring.R)[1], gens(M.ring.R)[4]
        h = - twist//r * (c2*M.ring(e1) + c1*M.ring(f1))
        H = OO(M,h)

        return all(chi(u*H) == 0 for u in M.bundles)
    end
end

function eulerform_UidualtimesH(Q::Quiver, d::Vector{Int}, twist::Int = 1)

    if gcd(d) != 1
        throw(ArgumentError("gcd(d,e) != 1; universal bundles don't exist"))
    end
    if !hassemistables(Q, d, canonical_stability(Q,d))
        throw(ArgumentError("Moduli space is empty"))
    end

    r = gcd(canonical_stability(Q,d))
    M = IntersectionTheory.matrix_moduli(narrows(Q),d[1],d[2])

    
    if twist == 0
        return all(chi(IntersectionTheory.dual(u)) == 0 for u in M.bundles)
    else
        c1 = Eulerform(Q,d,[1,0]) - Eulerform(Q,[1,0],d)
        c2 = Eulerform(Q,d,[0,1]) - Eulerform(Q,[0,1],d)
    
        e1,f1 = gens(M.ring.R)[1], gens(M.ring.R)[4]
    
        h = - twist//r * (c2*M.ring(e1) + c1*M.ring(f1))
        H = OO(M,h)
        
        return all(chi(IntersectionTheory.dual(u)*H) == 0 for u in M.bundles)
    end
end

function eulerform_H(Q::Quiver, d::Vector{Int}, twist::Int = 1)

    if gcd(d) != 1
        throw(ArgumentError("gcd(d,e) != 1; universal bundles don't exist"))
    end
    if !hassemistables(Q, d, canonical_stability(Q,d))
        throw(ArgumentError("Moduli space is empty"))
    end

    r = gcd(canonical_stability(Q,d))
    M = IntersectionTheory.matrix_moduli(narrows(Q),d[1],d[2])

    c1 = Eulerform(Q,d,[1,0]) - Eulerform(Q,[1,0],d)
    c2 = Eulerform(Q,d,[0,1]) - Eulerform(Q,[0,1],d)

    e1,f1 = gens(M.ring.R)[1], gens(M.ring.R)[4]
    h = - twist//r * (c2*M.ring(e1) + c1*M.ring(f1))
    H = OO(M,h)

    return chi(H) == 0
end


"""Checks wether the Teleman inequality holds for ``U_i \\otimes U_j \\otimes O(twist * H)."""
function teleman_twisted_endomorphisms_universal_bundle(Q::Quiver,d::Vector{Int},theta::Vector{Int},twist::Int = 1, slope_denominator::Function = sum)
    uidualuj = MegaQuiverTools.all_weights_endomorphisms_universal_bundle(Q, d, theta, slope_denominator)
    # this is actually -h because in the paper we twist by fractions of the anticanonical.
    h =  MegaQuiverTools.all_weights_irreducible_component_canonical(Q, d, theta, slope_denominator)
    bound = allTelemanbounds(Q, d, theta, slope_denominator)
    return all(maximum(uidualuj[key]) - twist * h[key] < bound[key] for key in keys(bound))
end

function teleman_twisted_universal_bundle(Q::Quiver,d::Vector{Int},theta::Vector{Int}, a::Vector{Int}, twist::Int = 1, slope_denominator::Function = sum)
    u = MegaQuiverTools.all_weights_universal_bundle(Q, d, theta, a, slope_denominator)
    # this is actually -h because in the paper we twist by fractions of the anticanonical.
    h =  MegaQuiverTools.all_weights_irreducible_component_canonical(Q, d, theta, slope_denominator)
    bound = allTelemanbounds(Q, d, theta, slope_denominator)
    return all(maximum(u[key]) - twist * h[key] < bound[key] for key in keys(bound))
end

function teleman_twisted_dual_universal_bundle(Q::Quiver,d::Vector{Int},theta::Vector{Int}, a::Vector{Int}, twist::Int = 1, slope_denominator::Function = sum)
    u = MegaQuiverTools.all_weights_universal_bundle(Q, d, theta, a, slope_denominator)
    # this is actually -h because in the paper we twist by fractions of the anticanonical.
    h =  MegaQuiverTools.all_weights_irreducible_component_canonical(Q, d, theta, slope_denominator)
    bound = allTelemanbounds(Q, d, theta, slope_denominator)
    return all(maximum(- u[key]) - twist * h[key] < bound[key] for key in keys(bound))
end

function teleman_H(Q::Quiver,d::Vector{Int},theta::Vector{Int},twist::Int = 1, slope_denominator::Function = sum)
    h =  MegaQuiverTools.all_weights_irreducible_component_canonical(Q, d, theta, slope_denominator)
    bound = allTelemanbounds(Q, d, theta, slope_denominator)
    return all(- twist * h[key] < bound[key] for key in keys(bound))
end

function main()

    # file_output = open("semiorthogonal-embeddings-output.txt", "w")
    # write(file_output, "This is the output of the semiorthogonal-embeddings.jl script.\n\n")
    # write(file_output, "It helps to check how many copies of D^b(Q) via the \"Fourier-Mukai-like\" functor in D^b(M) are semiorthogonal.\n")
    # write(file_output, "Concretely, it shows which Teleman vanishings can be applied and which Euler characteristics are zero.\n\n")
    # close(file_output)

    c = 50
    cases = [(mKroneckerquiver(m),[2,3]) for m in 3:c]
    a = [-1,1]

    for case in cases
        Q = case[1]
        d = case[2]
        theta = canonical_stability(Q, d)
        r = gcd(theta)


        Teleman_UidualUjH = [teleman_twisted_endomorphisms_universal_bundle(Q, d, theta, - twist) for twist in 1:r-2]
        Teleman_UidualtimesH = [teleman_twisted_dual_universal_bundle(Q, d, theta, a, - twist) for twist in 1:r-2]
        Teleman_H = [teleman_H(Q, d, theta, - twist) for twist in 1:r-2]
        Teleman_UitimesH = [teleman_twisted_universal_bundle(Q, d, theta, a, - twist) for twist in 0:r-2]


        Eulerchar_UidualUjH = [eulerform_UidualUjtimesH(Q, d, - twist) for twist in 1:r-1]
        Eulerchar_UidualtimesH = [eulerform_UidualtimesH(Q, d, - twist) for twist in 1:r-1]
        Eulerchar_H = [eulerform_H(Q, d, - twist) for twist in 1:r-1]
        Eulerchar_UitimesH = [eulerform_UitimesH(Q, d, - twist) for twist in 0:r-1]

        @info "Teleman vanishings for m = $r"
        @info Teleman_UidualUjH
        @info Teleman_UitimesH
        @info Teleman_UidualtimesH
        @info Teleman_H
        @info "Euler characteristics vanishings for m = $r"
        @info Eulerchar_UidualUjH
        @info Eulerchar_UitimesH
        @info Eulerchar_UidualtimesH
        @info Eulerchar_H
        
        # file_output = open("semiorthogonal-embeddings-output.txt", "a")
        # write(file_output, "Quiver: $Q, d: $d, Stability parameter: $theta, r = $(gcd(theta))\n")
        # write(file_output, "$Teleman_vanishings\n")
        # write(file_output, "$Eulerch_vanishings\n")
        # close(file_output)

    end
end

main()