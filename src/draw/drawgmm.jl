# using .GaussianMixtureAlignment: equal_volume_radius
# import GaussianMixtureAlignment.plot, GaussianMixtureAlignment.plot!
# export plot, plot!

const FEATURE_COLORS = 
    Dict(:donor         => RGB{N0f8}(1,   0,   1  ),  # magenta
         :acceptor      => RGB{N0f8}(0,   1,   0  ),  # green
         :cation        => RGB{N0f8}(1,   0,   0  ),  # red
         :anion         => RGB{N0f8}(0,   0,   1  ),  # blue
         :hydrophobe    => RGB{N0f8}(0,   1,   1  ),  # cyan
         :rings         => RGB{N0f8}(1,   0.5, 1  ),  # orange
         :aromaticring  => RGB{N0f8}(1,   0.25,0  ),  # brown
         :volume        => RGB{N0f8}(0.5, 0.5, 0.5),  # grey
    )

# function GaussianMixtureAlignment.plot!(m::PharmacophoreGMM; colors=FEATURE_COLORS, kwargs...)
#     stdcolors = Makie.wong_colors()
#     colorcount = 1
#     for (k,gmm) in m.gmms
#         color = stdcolors[colorcount]
#         if haskey(colors, k)
#             color = colors[k]
#         else
#             colorcount += 1
#         end
#         plot!(gmm; color=color, kwargs...)
#     end
# end

# function GaussianMixtureAlignment.plot!(ms::AbstractVector{<:PharmacophoreGMM}; kwargs...)
#     for m in ms
#         GaussianMixtureAlignment.plot!(m; kwargs...)
#     end
# end