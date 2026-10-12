# Run one case through the JIT app and capture the `ItsLiveInput` it hands `write_product`, in the raw
# format `rawio.jl` reads, beside the product the app writes. `writeproduct.jl` rebuilds the input from
# the capture and writes it again, so the two products can be compared.
#
#   julia --project=app -t 6,1 tools/trim/probes/capture_product.jl <sensor> <reference> <secondary> <outdir>
#
# `<outdir>/product.nc` is the app's product; `<outdir>/input/` the capture. Attribute dictionaries are
# written as `<key>\t<type>\t<value>` lines, one file per dictionary.

using AutoriftItsLiveProduction, ItsLiveOffsetProduction
const A = AutoriftItsLiveProduction
const P = ItsLiveOffsetProduction
include(joinpath(@__DIR__, "rawio.jl"))

sensor, ref, sec, out = ARGS
dir = joinpath(out, "input")
mkpath(dir)
isempty(readdir(dir)) || error("$dir is not empty")

function writeattrs(path::String, d::AbstractDict)
    open(path, "w") do io
        for k in sort!(collect(keys(d)))
            v = d[k]
            s = string(v)
            (occursin('\t', s) || occursin('\n', s) || occursin('\t', k)) && error("attribute $k does not fit one line")
            println(io, k, '\t', typeof(v), '\t', s)
        end
    end
end

function P.write_product(path::String, input::P.ItsLiveInput)
    g, info, c = input.georef, input.img_pair_info, input.coefficients
    writemeta(dir, "pair_type", input.pair_type)
    # Free text, one file each, so spaces survive.
    for k in (:detection_method, :coordinates, :autorift_software_version, :parameter_file, :source)
        write(joinpath(dir, string(k) * ".txt"), getfield(input, k))
    end
    writemeta(dir, "dt_seconds", input.dt_seconds === nothing ? "nothing" : input.dt_seconds)
    # `lonlat` transforms from the grid CRS, which `cf_grid_mapping` records as `spatial_epsg`.
    writemeta(dir, "grid_epsg", g.mapping_attrs["spatial_epsg"])
    writemeta(dir, "pixel_size", g.pixel_size_x, g.pixel_size_y)
    writeraw(dir, "x", reshape(g.x, :, 1))
    writeraw(dir, "y", reshape(g.y, :, 1))
    writeattrs(joinpath(dir, "mapping_attrs.tsv"), g.mapping_attrs)
    writemeta(dir, "acquisition_dates", info.acquisition_date_img1, info.acquisition_date_img2)
    writemeta(dir, "missions", info.mission_img1, info.mission_img2)
    writemeta(dir, "satellites", typeof(info.satellite_img1), info.satellite_img1,
              typeof(info.satellite_img2), info.satellite_img2)
    writemeta(dir, "roi_lat_lon", info.roi_valid_percentage, info.latitude, info.longitude)
    writeattrs(joinpath(dir, "extra.tsv"), info.extra)
    for k in (:offset2vx_1, :offset2vx_2, :offset2vy_1, :offset2vy_2, :scale_factor_1, :scale_factor_2)
        writeraw(dir, "c_" * string(k), getfield(c, k))
    end
    c.offset2vr === nothing || writeraw(dir, "c_offset2vr", c.offset2vr)
    c.offset2va === nothing || writeraw(dir, "c_offset2va", c.offset2va)
    writeraw(dir, "ref_vx", input.reference.vx)
    writeraw(dir, "ref_vy", input.reference.vy)
    writeraw(dir, "stable_mask", input.stable_mask)
    input.swath_bias === nothing || error("swath_bias capture is not written yet")
    writeraw(dir, "dx", input.dx)
    writeraw(dir, "dy", input.dy)
    writeraw(dir, "chip_size_x", input.chip_size_x)
    writemeta(dir, "scale_chip_size_y", input.scale_chip_size_y)
    writeraw(dir, "interp_mask", input.interp_mask)
    # The same input with a `roi_valid_percentage` that floors to zero, which selects the uncropped
    # schema: `product_p000.nc`.
    p000 = P.ImagePairInfo(info.acquisition_date_img1, info.acquisition_date_img2, info.mission_img1,
                           info.mission_img2, info.satellite_img1, info.satellite_img2, 0.5,
                           info.latitude, info.longitude, info.extra)
    fields = Any[getfield(input, k) for k in fieldnames(P.ItsLiveInput)]
    fields[findfirst(==(:img_pair_info), fieldnames(P.ItsLiveInput))] = p000
    invoke(P.write_product, Tuple{AbstractString,P.ItsLiveInput}, joinpath(out, "product_p000.nc"),
           P.ItsLiveInput(fields...))
    return invoke(P.write_product, Tuple{AbstractString,P.ItsLiveInput}, path, input)
end

c = A.Job(sensor, [ref], [sec])
A.run_job(c, joinpath(out, "product.nc"))
println("captured to $out")
