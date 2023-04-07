function parse_feature_definitions(path::String = joinpath(dirname(@__FILE__), "../../data/feature_definitions.txt"))
    atomtypes = Dict{Symbol,AtomType}()
    featuredefs = Dict{Symbol,FeatureDef}()
    continuing = false
    inatom = false
    atomname = ""
    negater = false
    infeature = false
    featurename = ""
    smarts = ""
    family = ""
    weights = Float64[]
    for line in eachline(path)
        isempty(line) && continue
        line = strip(line)
        first(line) == '#' && continue
        words = split(line, ' ')
        if inatom
            @assert continuing
            if last(line) == '\\'
                smarts = smarts * substitute_atom_types(string(line[1:end-1]), atomtypes)
            else
                smarts = smarts * substitute_atom_types(string(line), atomtypes)
                k = Symbol(atomname)
                if haskey(atomtypes, k)
                    atomtypes[k] = negater ? atomtypes[k] - smarts : atomtypes[k] + smarts
                else
                    push!(atomtypes, k => AtomType(smarts))
                end
                continuing = false
                inatom = false
                smarts = ""
            end
        elseif infeature
            if continuing
                if last(line) == '\\'
                    smarts = smarts * substitute_atom_types(string(line[1:end-1]), atomtypes)
                else
                    smarts = smarts * substitute_atom_types(string(line), atomtypes)
                    push!(atomtypes, Symbol(atomname) => AtomType(smarts))
                    continuing = false
                    smarts = ""
                end
            else
                words = split(line, ' ')
                if words[1] == "Family"
                    family = words[2]
                elseif words[1] == "Weights"
                    @show words[2]
                    weights = [parse(Float64, string(w)) for w in split(words[2], ',')]
                elseif words[1] == "EndFeature"
                    push!(featuredefs, Symbol(featurename) => FeatureDef(smarts, Symbol(family), weights))
                    smarts = ""
                    family = ""
                    weights = Float64[]
                    infeature = false
                else
                    throw(ErrorException("bad input line for feature: $(words[1])"))
                end
            end
        else
            words = split(line, ' ')
            if words[1] == "AtomType" || words[1] == "Atomtype"
                atomname = words[2]
                if first(atomname) == '!'
                    negater = true
                    atomname = atomname[2:end]
                end
                if last(line) == '\\'
                    smarts = smarts * substitute_atom_types(string(words[3]), atomtypes)
                    inatom = true
                    continuing = true
                else
                    k = Symbol(atomname)
                    smarts = substitute_atom_types(string(line), atomtypes)
                    if haskey(atomtypes, k)
                        @show k, smarts
                        atomtypes[k] = negater ? atomtypes[k] - smarts : atomtypes[k] + smarts
                    else
                        push!(atomtypes, k => AtomType(smarts))
                    end
                end
            elseif words[1] == "DefineFeature"
                featurename = words[2]
                smarts = substitute_atom_types(string(words[3]), atomtypes)
                infeature = true
                if last(line) == '\\'
                    continuing = true
                    smarts = string(words[3][1:end-1])
                else
                    smarts = string(words[3])
                end
            else
                throw(ErrorException("bad input line for feature: $(words[1])"))
            end
        end
    end
    return atomtypes, featuredefs
end

is_open_brace(x) = x == '{'
is_closed_brace(x) = x == '}'

function substitute_atom_types(trimmedline::String, atomtypes::Dict{Symbol,AtomType})
    smarts = trimmedline
    next_open_brace = findfirst(is_open_brace, smarts)
    next_closed_brace = 0
    while !isnothing(next_open_brace)
        next_closed_brace = findnext(is_closed_brace, smarts, next_open_brace)
        @show next_open_brace, next_closed_brace
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