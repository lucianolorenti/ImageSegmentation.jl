using Statistics
function compute_seeds(mapped_features::Array{Float64, 3}, number_of_seeds::Integer) where T
    rows, cols, mapped_dimension = size(mapped_features)
    N = rows * cols
    stepx = round(Integer, cols/sqrt(number_of_seeds))
    stepy = round(Integer, rows/sqrt(number_of_seeds))
    seeds_location = map(
        CartesianIndex,
        collect(
            Iterators.product(
                1:stepy:rows,
                1:stepx:cols)))
    number_of_seeds = length(seeds_location)
    R = CartesianIndices((rows, cols))
    Ifirst, Ilast = first(R), last(R)
    I1 = CartesianIndex(stepy, stepx)
    seeds = zeros(mapped_dimension, number_of_seeds) 
    for (i, I) in enumerate(seeds_location)
        range = max(Ifirst, I-I1):min(Ilast, I+I1)
        seeds[:, i] = mean(mapped_features[range, :], dims=[1, 2])
    end
    return seeds, seeds_location, stepy, stepx
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


#include<iostream>
#include<queue>
#include<vector>
#include<algorithm>
#include<float.h>
using namespace std;

class Superpixel

public:
	int Label;
	int Size;
	vector<int> Neighbor;
	Superpixel(int L=0,int S=0):Label(L),Size(S)end
	vector<int> xLoc;
	vector<int> yLoc;
	friend bool operator==(Superpixel& S,int L);
	friend bool operator==(int L,Superpixel& S);
end;

bool operator==(Superpixel& S,int L)

	return S.Label==L;
end

bool operator==(int L,Superpixel& S)

	return S.Label==L;
end


void EnforceConnectivity(mapped_features,
                         W::Matrix{Float64},
		float** L1,
		float** L2,
		float** a1,
		float** a2,
		float** b1,
		float** b2,
		float** x1,
		float** x2,
		float** y1,
		float** y2,
		double** W,
		unsigned short int* label,
		threshold::Integer,
	)

    mask = false(size(labels))
    sLabel = -1
    loc = []
    centers = []
    centerW = []
    stray = []
    struct Center
        data
        W
        I
    end
    function update!(c::Center, I, data, W)
        c.I += I
        c.data += data
        c.W += W
    end
	for I in CartesianIndices(labels)
        if mask[J]			
            continue
        end
		sLabel += 1
        Count = 1;
        push!(centers, Center())
        #zeros(10))
        #push!(centerW, 0)
        
        update!(centers[sLabel], I, mapped_features[:, I], W[I])
        # centers[sLabel]	+= mapped_features[:, I] * W[I]
        # centerW[sLabel] += W[I]       
		L = labels[I]
        labels[I] = sLabel
        mask[I] = true
        push!(loc, I)
        while !isempty(loc)
            J = pop!(loc)
            for K in max():min()
                if !mask[K] && labels[K] == L
                    Count += 1
                    push!(loc, K)
                    mask[K] = true
                    labels[K] = sLabel
                    update!(centers[sLabel], K, mapped_features[:, K], W[K])
                end
            end
        end
        push!(Size, Count)
        centers[sLabel] ./= centerW[sLabel]
	end
	sLabel += 1
    Count = 0
    
    Sarray = []
    Pointer = []
    for i=1:sLabel
        if Size[i] < threshold
            I = stray[i]
			L = labels[I]
			mask[I] = false
			indexMark = 0
			
			Superpixel S(L,Size[i]);
			S.xLoc.insert(S.xLoc.end(),x);
			S.yLoc.insert(S.yLoc.end(),y);
			while (indexMark<S.xLoc.size())
			
				x=S.xLoc[indexMark];y=S.yLoc[indexMark];
				indexMark++;
				int minX=(x-1<=0)?0:x-1;
				int maxX=(x+1>=nRows-1)?nRows-1:x+1;
				int minY=(y-1<=0)?0:y-1;
				int maxY=(y+1>=nCols-1)?nCols-1:y+1;
				for(int m=minX;m<=maxX;m++)
					for(int n=minY;n<=maxY;n++)
					
						if(mask[m][n]==1&&label[m*nCols+n]==L)
						
							mask[m][n]=0;
							S.xLoc.insert(S.xLoc.end(),m);
							S.yLoc.insert(S.yLoc.end(),n);
						end
						else if(label[m*nCols+n]!=L)
						
							int NewLabel=label[m*nCols+n];
							Pointer=find(S.Neighbor.begin(),S.Neighbor.end(),NewLabel);
							if(Pointer==S.Neighbor.end())
							
								S.Neighbor.insert(S.Neighbor.begin(),NewLabel);
							end
						end
					end

			end
			Sarray.insert(Sarray.end(),S);
		end
	end

	vector<Superpixel>::iterator S;
	vector<int>::iterator I;
	vector<int>::iterator I2;
	S=Sarray.begin();
	while(S!=Sarray.end())
	
		double MinDist=DBL_MAX;
		int Label1=(*S).Label;
		int Label2=-1;
		for(I=(*S).Neighbor.begin();I!=(*S).Neighbor.end();I++)
		
			double D=(centerL1[Label1]-centerL1[*I])*(centerL1[Label1]-centerL1[*I])+
				(centerL2[Label1]-centerL2[*I])*(centerL2[Label1]-centerL2[*I])+
				(centera1[Label1]-centera1[*I])*(centera1[Label1]-centera1[*I])+
				(centera2[Label1]-centera2[*I])*(centera2[Label1]-centera2[*I])+
				(centerb1[Label1]-centerb1[*I])*(centerb1[Label1]-centerb1[*I])+
				(centerb2[Label1]-centerb2[*I])*(centerb2[Label1]-centerb2[*I])+
				(centerx1[Label1]-centerx1[*I])*(centerx1[Label1]-centerx1[*I])+
				(centerx2[Label1]-centerx2[*I])*(centerx2[Label1]-centerx2[*I])+
				(centery1[Label1]-centery1[*I])*(centery1[Label1]-centery1[*I])+
				(centery2[Label1]-centery2[*I])*(centery2[Label1]-centery2[*I]);
			if(D<MinDist)
			
				MinDist=D;
				Label2=(*I);
			end
		end
		double W1=centerW[Label1];
		double W2=centerW[Label2];
		double W=W1+W2;
		centerL1[Label2]=(W2*centerL1[Label2]+W1*centerL1[Label1])/W;
		centerL2[Label2]=(W2*centerL2[Label2]+W1*centerL2[Label1])/W;
		centera1[Label2]=(W2*centera1[Label2]+W1*centera1[Label1])/W;
		centera2[Label2]=(W2*centera2[Label2]+W1*centera2[Label1])/W;
		centerb1[Label2]=(W2*centerb1[Label2]+W1*centerb1[Label1])/W;
		centerb2[Label2]=(W2*centerb2[Label2]+W1*centerb2[Label1])/W;
		centerx1[Label2]=(W2*centerx1[Label2]+W1*centerx1[Label1])/W;
		centerx2[Label2]=(W2*centerx2[Label2]+W1*centerx2[Label1])/W;
		centery1[Label2]=(W2*centery1[Label2]+W1*centery1[Label1])/W;
		centery2[Label2]=(W2*centery2[Label2]+W1*centery2[Label1])/W;
		centerW[Label2]=W;
		for(int i=0;i<(*S).xLoc.size();i++)
		
			int x=(*S).xLoc[i];int y=(*S).yLoc[i];
			label[x*nCols+y]=Label2;
		end
		vector<Superpixel>::iterator Stmp;
		Stmp=find(Sarray.begin(),Sarray.end(),Label2);
		if(Stmp!=Sarray.end())
		
			Size[Label2]=Size[Label1]+Size[Label2];
			if(Size[Label2]>=threshold)
			
				Sarray.erase(Stmp);
				Sarray.erase(S);
			end
			else
			
				(*Stmp).xLoc.insert((*Stmp).xLoc.end(),(*S).xLoc.begin(),(*S).xLoc.end());
				(*Stmp).yLoc.insert((*Stmp).yLoc.end(),(*S).yLoc.begin(),(*S).yLoc.end());
				(*Stmp).Neighbor.insert((*Stmp).Neighbor.end(),(*S).Neighbor.begin(),(*S).Neighbor.end());
				sort((*Stmp).Neighbor.begin(),(*Stmp).Neighbor.end());
				I=unique((*Stmp).Neighbor.begin(),(*Stmp).Neighbor.end());
				(*Stmp).Neighbor.erase(I,(*Stmp).Neighbor.end());
				I=find((*Stmp).Neighbor.begin(),(*Stmp).Neighbor.end(),Label1);
				(*Stmp).Neighbor.erase(I);
				I=find((*Stmp).Neighbor.begin(),(*Stmp).Neighbor.end(),Label2);
				(*Stmp).Neighbor.erase(I);
				Sarray.erase(S);
			end
		end
		else
		
			Sarray.erase(S);
		end
		for(int i=0;i<Sarray.size();i++)
		
			I=find(Sarray[i].Neighbor.begin(),Sarray[i].Neighbor.end(),Label1);
			I2=find(Sarray[i].Neighbor.begin(),Sarray[i].Neighbor.end(),Label2);
			if(I!=Sarray[i].Neighbor.end()&&I2!=Sarray[i].Neighbor.end())
				Sarray[i].Neighbor.erase(I);
			else if(I!=Sarray[i].Neighbor.end()&&I2==Sarray[i].Neighbor.end())
				(*I)=Label2;
		end
		S=Sarray.begin();
	end

	return;
end




function pre_enforce_connectivity(labels::Matrix{Int}; Bond::Integer=20)
	adj = 0
    mask = false(size(labels))
    Ifirst, Ilast = first(labels), last(labels)
    I1 = oneunit(Ifirst)
    loc = []
    for I in CartesianIndices(labels)
        if !mask[I]
            L = labels[I]
            for J in max(Ifirst, I-I1):min(Ilast, I+I1)
                if mask[J] && labels[J] == L
                    adj = labels[J]
                    break
                end
            end
            mask[I] = true
            push!(loc, I)
            indexMarker = 0 
            while indexMarker < length(loc)
                indexMarker += 1
                for J in max(Ifirst, I-I1):min(Ilast, I+I1)
                    if !mask[J] && labels[J]==L
                        mask[J] = true
                        push!(loc, J)
                    end
                end
            end
            if indexMarker<Bond
                for J in loc
                    labels[J] = adj
                end
            end
            loc = []
        end
    end
end

#endif

"""
```
```
"""
function LSC(img::Array{CT,2}, nseeds::Integer;  color_importance=0.7, space_importance=0.5, maxit::Int=50) where CT
    img = convert.(Lab, img)
    rows, cols = size(img)
    mapped_features = map_features(img, color_importance, space_importance)    
    sigma = dropdims(mean(mapped_features, dims=[1, 2]), dims=(1, 2))
    W = zeros(size(img))
    @inbounds for I in CartesianIndices(img)
        W[I] = mapped_features[I, :]' * sigma
        mapped_features[I, :] ./= W[I]
    end
    seeds, seeds_location, stepy, stepx = compute_seeds(mapped_features, nseeds)
    L = zeros(Integer, size(img))  
    nseeds = length(seeds_location)
    clusterSize = zeros(nseeds)
    WSum = zeros(nseeds)
    R = CartesianIndices(img)
    Ifirst, Ilast = first(R), last(R)
    I1 = CartesianIndex(stepy, stepx)
    D =  ones(size(img))
    for it = 1:maxit
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

    num_segments = nseeds
    TM = meantype(CT)
    result = similar(img, Int)
    region_means = Dict{Int, TM}()
    region_pix_count = Dict{Int, Int}()

    cluster_idx = 0
    for cluster in 1:nseeds
        cluster_idx += 1
        region_pix_count[cluster_idx] = clusterSize[cluster]
        pixels = channelview(img[L .== cluster])
        mean_value = dropdims(mean(pixels, dims=[2]), dims=(2))
        if any(isnan.(mean_value))
            mean_value = [0.0, 0.0, 0.0]
        end
        region_means[cluster_idx] = convert(CT, Lab(mean_value...))
    end
    for k in sort(keys(clusterSize))
        println(k, " ", clusterSize[k])
    end
    labels = sort(collect(keys(region_pix_count)))
    seg = SegmentedImage(L, labels, region_means, region_pix_count)
    seg = prune_segments(seg, 
                         i->(segment_pixel_count(seg,i)<50), 
                         (i,j)->(-segment_pixel_count(seg,j)))
    return seg
end
