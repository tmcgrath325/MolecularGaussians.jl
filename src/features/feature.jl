using MolecularGraph: smartstomol

const TEST_MOL = MolecularGraph.sdftomol(joinpath(dirname(@__FILE__), "../../assets/data/E1050_3d.sdf"))


struct AtomType
    smarts::String
    function AtomType(smarts::String, wrap = true)
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

Base.:+(x::AtomType, y::AtomType) = AtomType(x.smarts * "," * y.smarts, false)
Base.:+(x::AtomType, y::String) = x + AtomType(y)

Base.:-(x::AtomType, y::AtomType) = AtomType("!" * y.smarts * ";" * x.smarts, false)
Base.:-(x::AtomType, y::String) = x - AtomType(y)


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


function test_smarts(smarts::String)
    try 
        query = smartstomol(smarts)
        matches = MolecularGraph.substructmatches(TEST_MOL, query)
        collect(matches)
    catch
        return false
    end
    return true
end
