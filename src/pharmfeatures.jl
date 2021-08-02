function pharmstr2pair(str) 
    substr = split(str)
    key = Symbol(substr[length(substr)])
    val = Set{Int}()
    for i = 2:length(substr)-1
        push!(val, parse(Int,substr[i]))
    end     
    return Pair{Symbol,Set{Int}}(key, val)
end

function pharmfeatures(mol::UndirectedGraph, nodes::Set{Int}=nodeset(mol))
    strarray = try attributes(mol, :PUBCHEM_PHARMACOPHORE_FEATURES)
    catch e
        ""
    end
    feats = Dict{Symbol,Vector{Set{Int}}}()
    for i=2:length(strarray)
        feat = pharmstr2pair(strarray[i])
        if feat.second ∩ nodes == feat.second
            if haskey(feats, feat.first)
                push!(feats[feat.first], feat.second)
            else
                push!(feats, Pair(feat.first, [feat.second]))
            end
            # rings are not specified as aromatic in PubChem SDF files, so check here for aromaticity
            if feat.first==:rings && sum(isaromatic(mol)[collect(feat.second)]) == length(feat.second)
                feat = Pair(:aromaticrings, feat.second)
                if haskey(feats, feat.first)
                    push!(feats[feat.first], feat.second)
                else
                    push!(feats, Pair(feat.first, [feat.second]))
                end
            end
        end
    end
    return feats
end

pharmfeatures(submol::SubgraphView) = pharmfeatures(submol.graph, nodeset(submol))

pharmfeatures!(mol::UndirectedGraph) = attributes(mol)[:pharmfeatures] = pharmfeatures(mol)

function pharmfeatures!(mol::SubgraphView)
    pmol = mol.graph
    ppf = get!(attributes(pmol), :pharmfeatures, pharmfeatures!(pmol))
    return filter(pr -> pr.first ∈ nodeset(mol), ppf)
end
