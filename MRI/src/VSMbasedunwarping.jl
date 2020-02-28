# VSMBasedUnwarping.jl
# Methods for unwarping geometric distortions in MRI images

using Interpolations

function unwarp(VSM, distorted, dim)
    if dim == 2
        distorted = switchdim(distorted)
        VSM = switchdim(VSM)
    end

    unwarped = unwarp(VSM, distorted)

    if dim == 2
        unwarped = switchdim(unwarped)
    end
    unwarped
end

switchdim(v) = permutedims(v, [2, 1, (3:ndims(v))...])

unwarp(VSM, distorted) = unwarp!(similar(distorted), VSM, distorted)

function unwarp!(unwarped, VSM, distorted)
    xi = 1:size(distorted, 1)
    for J in CartesianIndices(size(distorted)[4:end])
        for I in CartesianIndices(size(distorted)[2:3])
            xtrue = xi .+ VSM[:,I]
            xregrid = (xtrue[1] .<= xi .<= xtrue[end]) # only use x values inside (no extrapolation)
            unwarped[.!xregrid,I,J] .= 0
            unwarped[xregrid,I,J] .= unwarpline(xtrue, distorted[:,I,J], xi[xregrid])
        end
    end
    unwarped
end

function unwarpline(xtrue, distorted, xnew)
    #TODO try better interpolation than linear
    interpolate((xtrue,), distorted, Gridded(Linear()))(xnew)
end

function getVSM(B0, rbw, dim, threshold = 5.0)
    VSM = B0 ./ (2π * rbw)
    thresholdforward(VSM, -0.9, threshold, dim)
end

function thresholdforward(VSM, tmin, tmax, dim)
    if dim == 2
        VSM = switchdim(VSM)
    end

    deltaVSM = VSM[2:end,:,:] .- VSM[1:(end-1),:,:]
    VSMret = copy(VSM)
    nx = size(VSM, 1)
    for I in CartesianIndices(size(VSM)[2:3])
        for x in 1:(nx-1)
            if deltaVSM[x,I] < tmin || deltaVSM[x,I] > tmax
                if deltaVSM[x,I] < tmin
                    diff = tmin - deltaVSM[x,I]
                else
                    diff = tmax - deltaVSM[x,I]
                end
                VSMret[x+1,I] += diff
                if x + 1 < nx
                    deltaVSM[x+1,I] -= diff
                end
            end
        end
    end

    if dim == 2
        VSMret = switchdim(VSMret)
    end
    VSMret
end

# TODO not properly working
threshold(VSM, tmin, tmax, dim) = (thresholdforward(VSM, tmin, tmax, dim) .+ thresholdforward(VSM[end:-1:1,:,:], -tmax, -tmin, dim))[end:-1:1,:,:] ./ 2.0
