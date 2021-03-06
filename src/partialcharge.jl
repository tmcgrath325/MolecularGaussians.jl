# Reads partial charge information for molecules loaded from the PubChem database
function partialcharges(mol::UndirectedGraph, nodes::Set{Int}=nodeset(mol))
    str2pair(str) = Pair{Int,Float64}(parse(Int, str[1]),parse(Float64, str[2]))
    pcstrarray = try attributes(mol, :PUBCHEM_MMFF94_PARTIAL_CHARGES)
    catch e
        " "
    end
    pcs = Dict{Int,Float64}()
    for i=2:length(pcstrarray)
        p = split(pcstrarray[i]) |> str2pair
        if p.first ∈ nodes
            push!(pcs, p)
        end
    end
    # atoms with no charge need to be added
    for idx in nodes
        get!(pcs, idx, 0)
    end
    return pcs
end

partialcharges(submol::SubgraphView) = partialcharges(submol.graph, nodeset(submol))

"""
    partialcharge!(mol)

Add an attribute `:partialcharge` to `mol`, a dictionary with the atom indices as keys and the MMFF94 
partial charges as values.
"""
partialcharges!(mol::UndirectedGraph) = attributes(mol)[:partialcharges] = partialcharges(mol)

function partialcharges!(mol::SubgraphView)
    pmol = mol.graph
    ppc = get!(attributes(pmol), :partialcharges, partialcharges!(pmol))
    return filter(pr -> pr.first ∈ nodeset(mol), ppc)
end
