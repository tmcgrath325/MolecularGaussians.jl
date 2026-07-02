"""
    AtomType(smarts::AbstractString; wrap=true)

A named atom class defined by the SMARTS atom expression `smarts` (the contents
of a bracketed atom, without the enclosing `[...]`).

By default the expression is wrapped as `\$(smarts)` so that it can be combined
with others via [`smarts_or`](@ref) and [`smarts_andnot`](@ref) without
precedence surprises. Pass `wrap=false` to store `smarts` verbatim; the
combining functions use this to avoid re-wrapping already-grouped expressions.
"""
struct AtomType
    smarts::String
    function AtomType(smarts::AbstractString; wrap = true)
        return wrap ? new("\$(" * smarts  * ")") : new(smarts)
    end
end

"""
    FeatureDef(smarts::AbstractString, family::Symbol, weights::AbstractVector{<:Real})

A pharmacophore feature definition: atoms matching the SMARTS pattern `smarts`
form a feature of pharmacophore `family` (e.g. `:Donor`, `:Acceptor`). `weights`
gives a per-matched-atom weighting used when the matched atoms are combined into
a single Gaussian feature.
"""
struct FeatureDef
    smarts::String
    family::Symbol
    weights::Vector{Float64}
    function FeatureDef(smarts::AbstractString, family::Symbol, weights::AbstractVector{<:Real})
        new(smarts, family, weights)
    end
end

"""
    FamilyDef(atomtypes, features, families)

A parsed set of pharmacophore definitions, as returned by
[`parse_feature_definitions`](@ref).

- `atomtypes::Dict{Symbol, AtomType}`: named atom classes referenced by feature
  SMARTS patterns.
- `features::Dict{Symbol, FeatureDef}`: named feature definitions.
- `families::Dict{Symbol, Vector{Symbol}}`: the feature names belonging to each
  pharmacophore family.
"""
struct FamilyDef
    atomtypes::Dict{Symbol, AtomType}
    features::Dict{Symbol, FeatureDef}
    families::Dict{Symbol, Vector{Symbol}}
end

"""
    smarts_or(x::AtomType, y) -> AtomType

Combine two atom types with SMARTS disjunction: the result matches an atom
satisfying `x` or `y` (the `,` operator in a SMARTS atom expression). `y` may be
an `AtomType` or a raw SMARTS `String`.
"""
smarts_or(x::AtomType, y::AtomType) = AtomType(x.smarts * "," * y.smarts; wrap=false)
smarts_or(x::AtomType, y::AbstractString) = smarts_or(x, AtomType(y))

"""
    smarts_andnot(x::AtomType, y) -> AtomType

Restrict `x` to atoms that do not also satisfy `y` (the SMARTS expression
`!y;x`). `y` may be an `AtomType` or a raw SMARTS `String`.
"""
smarts_andnot(x::AtomType, y::AtomType) = AtomType("!" * y.smarts * ";" * x.smarts; wrap=false)
smarts_andnot(x::AtomType, y::AbstractString) = smarts_andnot(x, AtomType(y))

Base.@deprecate Base.:+(x::AtomType, y::AtomType) smarts_or(x, y) false
Base.@deprecate Base.:+(x::AtomType, y::String) smarts_or(x, y) false
Base.@deprecate Base.:-(x::AtomType, y::AtomType) smarts_andnot(x, y) false
Base.@deprecate Base.:-(x::AtomType, y::String) smarts_andnot(x, y) false

# `public` is a soft keyword only on Julia 1.11+; guard so the module still
# loads on the 1.10 LTS, where the declaration is simply absent.
@static if VERSION >= v"1.11"
    eval(Meta.parse("public smarts_or, smarts_andnot, FamilyDef"))
end


smarts(f::AtomType) = "[" * f.smarts * "]"
smarts(f::FeatureDef) = f.smarts

Base.show(io::IO, f::AtomType) = print(io,
    summary(f),
    " with SMARTS string: $(smarts(f))",
)

Base.show(io::IO, f::FeatureDef) = print(io,
    summary(f),
    " of the :$(f.family) family with SMARTS string: $(smarts(f))",
)

Base.show(io::IO, ::MIME"text/plain", f::FeatureDef) = print(io,
    summary(f),
    " of the :$(f.family) family\n",
    "with SMARTS string: $(smarts(f))",
)
