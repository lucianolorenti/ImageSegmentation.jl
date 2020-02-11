
@testset "Linear Spectral Clustering" begin

    img = fill(Lab(0,0,0), 15, 15)
    img[6:10,4:8] .= Lab(0.5, 0.5, 0.5)
    img[3:7,2:6] .= Lab(0.8, 0.8, 0.8)
    mapped_features = ImageSegmentation.map_features(img, 0.8, 0.2)

    data = zeros(6, 6, 2)
    for I in CartesianIndices((6, 6))
        data[I, :] = [I[1], I[2]]
    end
    seeds = ImageSegmentation.compute_seeds(data, 5, 2)
    @test seeds == [2.0 2.5 3.0 4.0 4.5 5.0;
              2.0 2.5 3.0 4.0 4.5 5.0]

    using ImageSegmentation, Images, Test, TestImages
    img = testimage("mandrill")
    segments = ImageSegmentation.LSC(img, 50)
end