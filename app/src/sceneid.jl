# Scenes named by their provider's scene ID alone, resolved to the input the chain reads: derived from the
# ID where the provider's layout fixes it, else looked up in NASA's CMR catalog.
#
#   Landsat              LC08_L1TP_021248_20210507_20210517_02_T1        its landsatlook STAC item
#   Sentinel-2           S2A_MSIL1C_20200624T181931_N0209_R127_T19XDJ_20200624T201221
#                                                                        its B08.jp2 on the Google Cloud mirror
#   Sentinel-1           S1A_IW_SLC__1SDV_20211118T161558_…_04D1E5_D6DE  its zip at ASF
#   Sentinel-1 burst     S1_342239_IW3_20211202T161351_VV_3876-BURST     the burst of an SLC zip CMR names
#   NISAR                NISAR_L1_PR_RSLC_008_155_D_055_…_001            its .h5 at ASF, as CMR lists it
#
# `scene_source` takes the ID of one scene. A burst ID is ESA's relative burst ID, subswath, time,
# polarization and the SLC's checksum; which SLC holds the burst, and at what index, only CMR's burst record
# says, so a burst pair's IDs are resolved before `scene_source`, both scenes' in one query, to bursts of
# SLC zips. A NISAR granule's directory is its collection's, whose version the ID does not carry.

const LANDSAT_STAC_ITEMS = "https://landsatlook.usgs.gov/stac-server/collections/landsat-c2l1/items/"
const ASF_SLC = "https://datapool.asf.alaska.edu/SLC/"
const CMR_GRANULES = "https://cmr.earthdata.nasa.gov/search/granules.umm_json?provider=ASF&page_size=2000"

const LANDSAT_NAME = raw"L[TEOC]0[45789]_L1[A-Z]{2}_\d{6}_\d{8}_\d{8}_\d{2}_[A-Z0-9]{2}"
const LANDSAT_ID = Regex("^" * LANDSAT_NAME * "\$")
const LANDSAT_BAND_FILE = Regex("^(" * LANDSAT_NAME * ")_B(\\d+)\\.TIF\$", "i")
const S2_ID = r"^S2[ABC]_MSIL1C_\d{8}T\d{6}_N\d{4}_R\d{3}_T\d{2}[A-Z]{3}_\d{8}T\d{6}$"
const S1_SLC_ID = r"^S1[A-D]_IW_SLC__1S[SD][HV]_\d{8}T\d{6}_\d{8}T\d{6}_\d{6}_[0-9A-F]{6}_[0-9A-F]{4}$"
const S1_BURST_ID = r"^S1_\d{6}_IW[1-3]_\d{8}T\d{6}_(HH|HV|VH|VV)_[0-9A-F]{4}-BURST$"
const NISAR_ID = r"^NISAR_L[12]_[A-Z]{2}_(RSLC|GSLC)_[0-9A-Z_]+$"

# A Sentinel-1 SLC's ID, which CMR writes with an `-SLC` suffix, without it.
s1_slc_name(id::AbstractString) = String(chopsuffix(id, "-SLC"))
is_s1_slc_id(id::AbstractString) = occursin(S1_SLC_ID, s1_slc_name(id))

"""
    s1_slc_url(name) -> String

Where ASF serves Sentinel-1 SLC `name`'s zip: under `SLC/SA` for Sentinel-1A, `SB` for 1B, and so on.
"""
s1_slc_url(name::String) = ASF_SLC * "S" * name[3:3] * "/" * name * ".zip"

"""
    resolve_burst_ids(reference, secondary) -> (String, String)

A burst pair's two scenes with their IDs replaced: burst IDs, comma-separated, by the bursts of their SLC
zips (`<zip>#IW<n>/<POL>/<index>`), looked up in one query for both scenes; SLC IDs by their zips. A scene
that names no ID is returned as given.
"""
function resolve_burst_ids(reference::String, secondary::String)
    r, s = _id_parts(reference), _id_parts(secondary)
    bursts = String[p for p in [r; s] if occursin(S1_BURST_ID, p)]
    specs = Dict{String,String}(zip(bursts, isempty(bursts) ? String[] : burst_specs(bursts)))
    resolve(input, parts) = all(p -> occursin(S1_BURST_ID, p), parts) ? join([specs[p] for p in parts], ",") :
                            all(is_s1_slc_id, parts) ? join([s1_slc_url(s1_slc_name(p)) for p in parts], ",") :
                            input
    return resolve(reference, r), resolve(secondary, s)
end

# A scene's comma-separated parts, refused when IDs are mixed with paths or URLs or with IDs of another kind.
function _id_parts(input::String)
    parts = String[String(strip(p)) for p in split(input, ','; keepempty = false)]
    nb = count(p -> occursin(S1_BURST_ID, p), parts)
    ns = count(is_s1_slc_id, parts)
    (nb + ns == 0 || nb == length(parts) || ns == length(parts)) || throw(ArgumentError(
        "\"" * input * "\" mixes burst IDs, SLC IDs, paths or URLs; name a scene's parts all one way"))
    return parts
end

# ---------------------------------------------------------------------------
# CMR
# ---------------------------------------------------------------------------
#
# A granule search's UMM-G JSON, only the fields read here. `Mutable` structs, as for STAC items: keys are
# matched by name in any order, and keys without a field are skipped.

mutable struct CmrAttribute
    name::String
    values::Vector{String}
    CmrAttribute() = new("", String[])
end
mutable struct CmrUrl
    url::String
    type::String
    CmrUrl() = new("", "")
end
mutable struct CmrGranule
    granule_ur::String
    attributes::Vector{CmrAttribute}
    urls::Vector{CmrUrl}
    inputs::Vector{String}
    CmrGranule() = new("", CmrAttribute[], CmrUrl[], String[])
end
mutable struct CmrItem
    umm::CmrGranule
    CmrItem() = new(CmrGranule())
end
mutable struct CmrSearch
    items::Vector{CmrItem}
    CmrSearch() = new(CmrItem[])
end

StructTypes.StructType(::Type{CmrAttribute}) = StructTypes.Mutable()
StructTypes.StructType(::Type{CmrUrl}) = StructTypes.Mutable()
StructTypes.StructType(::Type{CmrGranule}) = StructTypes.Mutable()
StructTypes.StructType(::Type{CmrItem}) = StructTypes.Mutable()
StructTypes.StructType(::Type{CmrSearch}) = StructTypes.Mutable()
StructTypes.names(::Type{CmrAttribute}) = ((:name, :Name), (:values, :Values))
StructTypes.names(::Type{CmrUrl}) = ((:url, :URL), (:type, :Type))
StructTypes.names(::Type{CmrGranule}) =
    ((:granule_ur, :GranuleUR), (:attributes, :AdditionalAttributes), (:urls, :RelatedUrls), (:inputs, :InputGranules))

_cmr_error(g::CmrGranule, what::String) = ArgumentError("CMR's record of " * g.granule_ur * " " * what)

"""
    cmr_granules(doc, ids) -> Vector{CmrGranule}

The granules `ids` name, in that order, out of CMR granule search response `doc`. Throws naming any ID the
response does not hold.
"""
function cmr_granules(doc::String, ids::Vector{String})
    found = Dict{String,CmrGranule}(item.umm.granule_ur => item.umm for item in JSON3.read(doc, CmrSearch).items)
    missing_ids = filter(id -> !haskey(found, id), ids)
    isempty(missing_ids) || throw(ArgumentError("ASF's catalog (CMR) lists no granule " * join(missing_ids, ", ")))
    return [found[id] for id in ids]
end

# The granules `ids` name, from one CMR search. IDs the patterns above accept need no escaping in a query.
cmr_lookup(ids::Vector{String}) = cmr_granules(read_text(CMR_GRANULES * join(["&granule_ur=" * id for id in unique(ids)])), ids)

function cmr_attribute(g::CmrGranule, name::String)
    for a in g.attributes
        a.name == name && length(a.values) == 1 && return a.values[1]
    end
    throw(_cmr_error(g, "has no single " * name))
end

"""
    burst_spec(g) -> String

The burst CMR burst record `g` describes, as a burst of its SLC's zip (`<zip>#IW<n>/<POL>/<index>`, the
index counted from zero, as CMR counts it).
"""
function burst_spec(g::CmrGranule)
    length(g.inputs) == 1 || throw(_cmr_error(g, "names no single SLC"))
    m = match(r"^IW([1-3])$", cmr_attribute(g, "SUBSWATH_NAME"))
    m === nothing && throw(_cmr_error(g, "has no IW subswath"))
    k = tryparse(Int, cmr_attribute(g, "BURST_INDEX"))
    k === nothing && throw(_cmr_error(g, "has no burst index"))
    b = SLCDatasets.ZipBurst(s1_slc_url(s1_slc_name(g.inputs[1])), parse(Int, something(m.captures[1])),
                             cmr_attribute(g, "POLARIZATION"), k + 1)
    return SLCDatasets.zip_burst_spec(b)
end

burst_specs(ids::Vector{String}) = [burst_spec(g) for g in cmr_lookup(ids)]

"""
    granule_h5(g) -> String

The one `https://` link to a `.h5` among CMR record `g`'s data links.
"""
function granule_h5(g::CmrGranule)
    hits = [u.url for u in g.urls if u.type == "GET DATA" && startswith(u.url, "https://") && endswith(u.url, ".h5")]
    length(hits) == 1 || throw(_cmr_error(g, "lists " * string(length(hits)) * " https:// .h5 links, not one"))
    return hits[1]
end

nisar_url(id::String) = granule_h5(only(cmr_lookup([id])))
