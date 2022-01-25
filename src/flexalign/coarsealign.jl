function closest_conformers(gmmx, gmmy; alignfun::Union{typeof(rocs_align),typeof(tiv_gogma_align)} = rocs_align, kwargs...)
    xconformers = conformers(gmmx; kwargs...)
    yconformers = conformers(gmmy; kwargs...)

    bestx = xconformers[1]
    besty = yconformers[1]
    bestolap = -overlap(bestx, besty)
    for xconf in xconformers
        for yconf in yconformers
            res = alignfun(xconf, yconf)
            min = typeof(alignfun) == typeof(rocs_align) ? res.minimum : res.lowerbound
            if min < bestolap
                bestx = res.x
                besty = res.y
                bestolap = min
            end
        end
    end
    return bestx, besty, bestolap
end