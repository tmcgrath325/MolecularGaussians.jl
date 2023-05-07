mutable struct FeatureDefParser
    lineiter::Base.Iterators.Stateful{Base.EachLine{IOStream}}
    line::String
    words::Vector{String}
    atomtypes::Dict{Symbol,AtomType}
    features::Dict{Symbol,FeatureDef}
    families::Dict{Symbol,Vector{Symbol}}
    function FeatureDefParser(path::String)
        lineiter = Iterators.Stateful(eachline(path))
        line = ""
        words = SubString{String}[]
        atomtypes = Dict{Symbol,AtomType}()
        features = Dict{Symbol,FeatureDef}()
        families = Dict{Symbol,Vector{Symbol}}()
        return new(lineiter, line, words, atomtypes, features, families)
    end
end

parse_feature_definitions(path::String = joinpath(dirname(@__FILE__), "../../assets/const/FeatureDefinitions.fdef")) = parse_feature_definitions!(FeatureDefParser(path))

function parse_feature_definitions!(parser::FeatureDefParser)
    while(!isempty(parser.lineiter))
        read_line!(parser)
        featuretype = lowercase(parser.words[1])
        if featuretype == "#" || featuretype == ""
            continue
        elseif featuretype == "atomtype"
            parse_atomtype!(parser)
        elseif featuretype == "definefeature"
            parse_featuredef!(parser)
        else
            throw(ErrorException("bad input line for feature: $(parser.words[1])"))
        end
    end
    return FamilyDef(parser.atomtypes, parser.features, parser.families)
end

function parse_atomtype!(parser::FeatureDefParser)
    atomname = Symbol(parser.words[2])
    negater = first(parser.words[2]) == '!'
    if negater
        atomname = Symbol(parser.words[2][2:end])
        @assert haskey(parser.atomtypes, atomname)
    end
    smarts = substitute_atom_types(parser.words[3], parser.atomtypes)
    while last(parser.line) == '\\'
        smarts = chop(smarts)
        read_line!(parser)
        @assert parser.words[1] == parser.line
        smarts = smarts * substitute_atom_types(parser.line, parser.atomtypes)
    end
    if haskey(parser.atomtypes, atomname)
        parser.atomtypes[atomname] = negater ? parser.atomtypes[atomname] - smarts : parser.atomtypes[atomname] + smarts
    else
        push!(parser.atomtypes, atomname => AtomType(smarts))
    end
end

function parse_featuredef!(parser::FeatureDefParser)
    featurename = Symbol(parser.words[2])
    family = nothing
    weights = nothing
    smarts = substitute_atom_types(parser.words[3], parser.atomtypes)
    while last(parser.line) == '\\'
        smarts = chop(smarts)
        read_line!(parser)
        @assert parser.words[1] == parser.line
        smarts = smarts * substitute_atom_types(parser.line, parser.atomtypes)
    end
    read_line!(parser)
    while lowercase(parser.words[1]) != "endfeature"
        if lowercase(parser.words[1]) == "family"
            family = Symbol(parser.words[2])
        elseif lowercase(parser.words[1]) == "weights"
            weights = [parse(Float64, string(w)) for w in split(parser.words[2], ',')]
        else
            throw(ErrorException("bad input line for feature: $(parser.words[1])"))
        end
        read_line!(parser)
    end
    push!(parser.features, featurename => FeatureDef(smarts, family, weights))
    if haskey(parser.families, family)
        push!(parser.families[family], featurename)
    else
        push!(parser.families, family => [featurename,])
    end
end

function read_line!(parser::FeatureDefParser)
    parser.line = popfirst!(parser.lineiter)
    parser.words = [string(s) for s in split(strip(parser.line), ' ')]
end

is_open_brace(x) = x == '{'
is_closed_brace(x) = x == '}'

function substitute_atom_types(trimmedline::String, atomtypes::Dict{Symbol,AtomType})
    smarts = trimmedline
    next_open_brace = findfirst(is_open_brace, smarts)
    next_closed_brace = 0
    while !isnothing(next_open_brace)
        next_closed_brace = findnext(is_closed_brace, smarts, next_open_brace)
        atomtypestr = smarts[next_open_brace+1:next_closed_brace-1]
        removelen = length(atomtypestr) + 2
        atomtypekey = Symbol(smarts[next_open_brace+1:next_closed_brace-1])
        insertsmarts = atomtypes[atomtypekey].smarts
        insertlen = length(insertsmarts)
        smarts = smarts[1:next_open_brace-1] * insertsmarts * smarts[next_closed_brace+1:end]
        lastidx = next_closed_brace + (insertlen - removelen)
        next_open_brace = findnext(is_open_brace, smarts, lastidx)
    end
    return smarts
end