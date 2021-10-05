rocs_align(molmoving::Union{UndirectedGraph, SubgraphView}, molfixed::Union{UndirectedGraph, SubgraphView}, σfun=ones, ϕfun=ones; kwargs...
    ) = rocs_align(MolGMM(molmoving, σfun, ϕfun), MolGMM(molfixed, σfun, ϕfun); kwargs...)