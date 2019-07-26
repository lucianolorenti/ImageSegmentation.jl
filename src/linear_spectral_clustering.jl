function refine_seeds(img::Array{CT,2}, seeds::Vector) where CT
    return seeds
end

function map_features(img::Array{Lab,2}, Cc::Float64, Cs::Float64)
    feature_map = zeros((size(img)..., 10))
    @inbounds for IA in CartesianIndices(img)
        feature_map[IA, 1] = Cc * cos(pi/2) * img.l
        feature_map[IA, 2] = Cc * sin(pi/2) * img.l
        feature_map[IA, 3] = 255 * Cc * cos(pi/2) * img.a
        feature_map[IA, 4] = 255 * Cc * sin(pi/2) * img.a
        feature_map[IA, 5] = 255 * Cc * cos(pi/2) * img.b
        feature_map[IA, 6] = 255 * Cc * sin(pi/2) * img.b
        feature_map[IA, 7] = Cs * cos(pi/2) * IA[2]
        feature_map[IA, 8] = Cs * sin(pi/2) * IA[2]
        feature_map[IA, 9] = Cs * cos(pi/2) * IA[1]
        feature_map[IA, 10]  = Cs * sin(pi/2) * IA[1]
    end
    return feature_map
end

"""
```
```
"""
function LSC(img::Array{CT,2}, K::Integer, range_radius::Real; iter::Int = 50, eps::Real = 0.01) where CT
    rows, cols = size(img)
    seeds = linspace(1, rows * cols, K)
    seeds = refine_seeds(img, seeds)
    L = zeros(size(img))
    D = ones(size(img)) * inf
    mapped_features = map_features(img, 0.5, 0.6)
    for (wm, sc) in (weighted_means, search_center)
        for 
            dist = norm(mapped_features[I, :] - weighted_mean)
            if dist < D[I]
                D[I] = dist
                L[I] = k
            end
        end
    end
    
    


end
