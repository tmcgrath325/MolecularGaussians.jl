using MolecularGraph: smartstomol


struct AtomType
    smarts::String
    function AtomType(str::String)
        smarts = str
        if first(smarts) == '['
            @assert last(smarts) == ']'
            smarts = smarts[2:end-1]
        end
        if smarts[1:2] == "\$("
            @assert last(smarts) == ')'
        else
            smarts = "\$(" * smarts * ")"
        end
        return new(smarts)
    end
end

struct FeatureDef
    smarts::String
    family::Symbol
    weights::Vector{Float64}
end

function test_smarts(smarts::String)
    mol = try smartstomol("[" * smarts * "]")
    catch e
        return false
    end
    return true
end

function Base.:+(x::AtomType, y::AtomType)
    xsmarts = strip(x.smarts, ('[',']'))
    ysmarts = strip(y.smarts, ('[',']'))
    AtomType("[" * xsmarts * "," * ysmarts * "]")
end
Base.:+(x::AtomType, y::String) = x + AtomType(y)


function Base.:-(x::AtomType, y::AtomType)
    xsmarts = strip(x.smarts, ('[',']'))
    ysmarts = strip(y.smarts, ('[',']'))
    return AtomType("[!" * ysmarts * ";" * xsmarts * "]")
end
Base.:+(x::AtomType, y::String) = x - AtomType(y)

Base.show(io::IO, a::AtomType) = println(io,
    summary(a),
    " with SMARTS string [$(a.smarts)]",
)