struct AtomType
    smarts::String
    function AtomType(smarts::String; wrap = true)
        return wrap ? new("\$(" * smarts  * ")") : new(smarts)
    end
end

struct FeatureDef
    smarts::String
    family::Symbol
    weights::Vector{Float64}
    function FeatureDef(smarts::String, family::Symbol, weights::Vector{Float64})
        new(smarts, family, weights)
    end
end

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
smarts_or(x::AtomType, y::String) = smarts_or(x, AtomType(y))

"""
    smarts_andnot(x::AtomType, y) -> AtomType

Restrict `x` to atoms that do not also satisfy `y` (the SMARTS expression
`!y;x`). `y` may be an `AtomType` or a raw SMARTS `String`.
"""
smarts_andnot(x::AtomType, y::AtomType) = AtomType("!" * y.smarts * ";" * x.smarts; wrap=false)
smarts_andnot(x::AtomType, y::String) = smarts_andnot(x, AtomType(y))

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

Base.show(io::IO, f::AtomType) = println(io,
    summary(f),
    " with SMARTS string: $(smarts(f))",
)

Base.show(io::IO, f::FeatureDef) = println(io,
    summary(f),
    " of the :$(f.family) family\n",
    "with SMARTS string: $(smarts(f))",
)
