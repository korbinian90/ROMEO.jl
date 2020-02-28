using Statistics#, Images

#TODO open not writable as standard
function readphase(fn; keyargs...)
    phase = niread(fn; keyargs...)
    if any(.!isfinite.(phase))
        println("WARNING: NaNs have been replaced with zeros in the image: $fn")
        phase.raw[.!isfinite.(phase)] .= 0
    end
    minp, maxp = minmaxmiddleslice(phase.raw)
    phase.header.scl_slope = 2pi / (maxp - minp)
    phase.header.scl_inter = -pi - minp * phase.header.scl_slope
    if any(isnan.(phase.raw)) println("WARNING: there are NaNs in the image: $fn") end
    return phase
end

function readmag(fn; keyargs...)
    mag = niread(fn; keyargs...)
    if any(.!isfinite.(mag))
        println("WARNING: NaNs have been replaced with zeros in the image: $fn")
        mag.raw[.!isfinite.(mag)] .= 0
    end
    if mag.header.scl_slope == 0
        _, maxi = minmaxmiddleslice(mag.raw)
        mag.header.scl_slope = 1 / maxi
        mag.header.scl_inter = 0
    end
    if any(isnan.(mag.raw)) println("WARNING: there are NaNs in the image: $fn") end
    return mag
end

function minmaxmiddleslice(image)
    slice = if ndims(image) ≥ 3 (
            ones = repeat([1], ndims(image)-3);
            middle = div(size(image, 3), 2);
            view(image,:,:,middle,ones...))
        else image
        end
     minimum(slice), maximum(slice)
end

savenii(image, name, writedir; kwargs...) = savenii(image, @show joinpath(writedir, name * ".nii"); kwargs...)
savenii(image, filepath; kwargs...) = niwrite(filepath, NIVolume([k[2] for k in kwargs]..., image))
savenii(image::BitArray, filepath; kwargs...) = niwrite(filepath, NIVolume([k[2] for k in kwargs]..., UInt8.(image)))

function createniiforwriting(im, name::AbstractString, writedir::AbstractString; datatype::DataType = Float64, header = NIVolume(zeros(datatype, 1)).header)
    if !occursin(r"\.nii$", name)
        name *= ".nii"
    end
    filepath = joinpath(writedir, name)
    createniiforwriting(im, filepath; datatype = datatype, header = header)
end

function createniiforwriting(im::AbstractArray, filepath::AbstractString; datatype::DataType = eltype(im), header = NIVolume(zeros(datatype, 1)).header)
    nii = createniiforwriting(size(im), filepath; datatype = datatype, header = header)
    nii .= im
end

createniiforwriting(sz::Tuple, filepath::AbstractString; datatype::DataType = Float64, header = NIVolume(zeros(datatype, 1)).header) = createniiforwriting([sz...], filepath; datatype = datatype, header = header)
function createniiforwriting(sz::AbstractVector{Int}, filepath::AbstractString; datatype::DataType = Float64, header = NIVolume(zeros(datatype, 1)).header)
    write_emptynii(sz, filepath, datatype = datatype, header = header)
    niread(filepath, mmap = true, write = true).raw
end

function getHIP(mag, phase; echoes = [1,3])
    e1, e2 = echoes
    compl = zeros(ComplexF64, size(mag)[1:3])
    for iCha in 1:size(mag, 5)
        compl .+= exp.(1.0im .* (phase[:,:,:,e2,iCha] .- phase[:,:,:,e1,iCha])) .* mag[:,:,:,e1,iCha] .* mag[:,:,:,e2,iCha]
    end
    compl
end


function robustmask!(image; maskedvalue = if eltype(image) <: AbstractFloat NaN else 0 end)
    image[.!getrobustmask(image)] .= maskedvalue
    image
end

function getrobustmask(weight)
    notnan = isfinite.(weight)
    mean1 = mean(weight[notnan])
    noisemask = weight .< mean1
    noisemean = mean(weight[noisemask])
    if !isfinite(noisemean) return ones(size(weight)) end

    signalmean = mean(weight[.!noisemask .& notnan])
    noise_σ = std(weight[noisemask])

    mask = weight .> noisemean + noise_σ
    #closing(mask) .& notnan
end

getcomplex(mag::NIVolume, phase::NIVolume) = getcomplex(mag.raw, phase.raw)
getcomplex(fnmag::AbstractString, fnphase::AbstractString) = getcomplex(niread(fnmag), niread(fnphase))

function getcomplex(mag, phase)
    higherdims = ones(Int, length(size(phase)) - 2)
    minp = minimum(phase[:,:,higherdims...])
    maxp = maximum(phase[:,:,higherdims...])

    mag .* exp.((2im * pi / (maxp - minp)) .* phase)
end

function readfromtextheader(filename, searchstring)
    for line in readlines(open(filename, "r"))
        if occursin(searchstring, line)
            # regex to search for "= " or ": " and return the following non-whitespace characters
            return match(r"(?<=(= |: ))(\S+)", line).match
        end
    end
end

# root sum of squares combination
RSS(mag; dim = ndims(mag)) = dropdims(.√sum(mag.^Float32(2); dims = dim); dims = dim)

function getscaledimage(array, div::Number, offset = 0, type::Symbol = :trans)
    array = reshape(array, size(array)[1:2]) # drops singleton dimensions
    scaled = if offset != 0
        (array .- offset) .* (1 / div) .+ 0.5
    else
        array .* (1 / div)
    end
    scaled[isnan.(scaled) .| (scaled .< 0)] .= 0
    scaled[scaled .> 1] .= 1
    if type == :trans
        scaled = reverse(permutedims(scaled, [2 1]); dims = 1)
    else
    end
    Gray.(scaled)
end

function getscaledimage(array, type::Symbol = :trans)
    scaled = robustrescale(array, 0, 1, threshold = true)
    getscaledimage(scaled, 1, 0, type)
end

robustrescale(array, newmin, newmax; threshold = false, mask = trues(size(array)), datatype = Float64) =
    robustrescale!(datatype.(array), newmin, newmax; threshold = threshold, mask = mask)

function robustrescale!(array, newmin, newmax; threshold = false, mask = trues(size(array)))
    array[isnan.(array)] .= minimum(array[.!isnan.(array)])
    q = [0.01, 0.99] # quantiles
    oldq = estimatequantile(array[mask], q)
    oldrange = (oldq[2] - oldq[1]) / (q[2] - q[1])
    oldmin = oldq[1] - q[1] * oldrange
    newrange = newmax - newmin

    array .= (array .- oldmin) .* (newrange / oldrange) .+ newmin

    if threshold
        array[array .< newmin] .= newmin
        array[array .> newmax] .= newmax
    end
    array
end

function estimatequantile(array, p)
    # only use 5 samples of each higher dimension
    # TODO take 100 x 16 consecutive with linear indexing (bound checking!)
    inds = [unique(round.(Int, LinRange(1, sz, 5))) for sz in size(array)[2:end]]
    quantile(array[:, inds...][:], p)
end

function rescale(array, newmin, newmax; datatype = eltype(array))
    rescale!(datatype.(array), newmin, newmax)
end

function rescale!(array, newmin, newmax)
    # TODO sample function to take only every.. elem block (faster)
    oldmin, oldmax = extrema(array)
    factor = (newmax - newmin) / (oldmax - oldmin)
    array .= (array .- oldmin) .* factor .+ newmin
end

mmtovoxel(sizemm, nii::NIVolume) = mmtovoxel(sizemm, nii.header)
mmtovoxel(sizemm, header::NIfTI1Header) = mmtovoxel(sizemm, header.pixdim)
mmtovoxel(sizemm, pixdim) = sizemm ./ pixdim
