using Statistics
function compute_seeds(mapped_features::Array{Float64, 3}, number_of_seeds::Integer, stepx::Integer, stepy::Integer) where T
    rows, cols, mapped_dimension = size(mapped_features)
    N = rows * cols
    seeds_indices = collect(range(1, N, step=Integer(round(N/number_of_seeds))))
    R = CartesianIndices((rows, cols))
    Ifirst, Ilast = first(R), last(R)
    I1 = CartesianIndex(stepx, stepy)
    seeds = zeros(mapped_dimension, length(seeds_indices))
    seeds_location = CartesianIndices((rows, cols))[seeds_indices]
    for (i, I) in enumerate(seeds_location)
        range = max(Ifirst, I-I1):min(Ilast, I+I1)
        seeds[:, i] = mean(mapped_features[range, :], dims=[1, 2])
    end
    return seeds, seeds_location
end

function map_features(img::Array{Lab{T},2}, Cc::Float64, Cs::Float64) where T
    feature_map = zeros((size(img)..., 10))
    @inbounds for IA in CartesianIndices(img)
        c = img[IA]
        feature_map[IA, 1] = Cc * cos(pi/2) * c.l
        feature_map[IA, 2] = Cc * sin(pi/2) * c.l
        feature_map[IA, 3] = 255 * Cc * cos(pi/2) * c.a
        feature_map[IA, 4] = 255 * Cc * sin(pi/2) * c.a
        feature_map[IA, 5] = 255 * Cc * cos(pi/2) * c.b
        feature_map[IA, 6] = 255 * Cc * sin(pi/2) * c.b
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
function LSC(img::Array{CT,2}, nseeds::Integer;  color_importance=0.7, space_importance=0.5, iter::Int = 50) where CT
    img = convert.(Lab, img)
    rows, cols = size(img)
    mapped_features = map_features(img, color_importance, space_importance)    
    sigma = dropdims(mean(mapped_features, dims=[1, 2]), dims=(1, 2))
    W = zeros(size(img))
    @inbounds for I in CartesianIndices(img)
        W[I] = mapped_features[I, :]' * sigma
        mapped_features[I, :] ./= W[I]
    end
    ColNum=sqrt(float(nseeds*rows/cols));
	RowNum=nseeds/ColNum;
	StepX=ceil(Integer, rows/RowNum) 
    StepY=ceil(Integer, cols/ColNum) 
    seeds, seeds_location = compute_seeds(mapped_features, nseeds, StepX, StepY)
    L = zeros(Integer, size(img))  
    clusterSize = zeros(nseeds)
    WSum = zeros(nseeds)
    R = CartesianIndices(img)
    Ifirst, Ilast = first(R), last(R)
    I1 = CartesianIndex(StepX, StepY)
    D =  ones(size(img))
    for it = 1:iter
        @info "Iteration $it"
        fill!(D, Inf)
        for (i, seed) in enumerate(seeds_location)
            for J in max(Ifirst, seed-I1):min(Ilast, seed+I1)
                dist = norm(mapped_features[J, :] - seeds[:, i])
                if dist < D[J]
                    D[J] = dist
                    L[J] = i
                end
            end
            seeds[:, i] .= 0
            clusterSize[i] = 0
            WSum[i] = 0
        end
        for I in CartesianIndices(img)
            label = L[I]        
            seeds[:, label] += mapped_features[I, :] * W[I]
            clusterSize[label] += 1
            WSum[label] += W[I]
            seeds_location[label] += I
        end
        for (i, seed) in enumerate(seeds_location)
            WSum[i]= (WSum[i]==0) ? 1 : WSum[i]
            clusterSize[i]= ( clusterSize[i] == 0 ) ? 1 : clusterSize[i]
            seeds[:, i] ./= WSum[i]
            new_seed_x = round(Integer, seeds_location[i][1] / clusterSize[i])
            new_seed_y = round(Integer, seeds_location[i][2] / clusterSize[i])
            seeds_location[i] = CartesianIndex(new_seed_x, new_seed_y)
        end
    end
    return L
end
