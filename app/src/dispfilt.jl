# The displacement filter (`DISP_FILT`, `AutoRIFT.GardnerFilter`) each sensor is correlated with.

# Each filter is the loosest whose added points have no more than 1.5 times the stable-surface error of the
# points (3.0, 0.15) keeps. Only `mad_scale` (`MadScalar`) and `agree_tolerance` (`FracSearch`) are set per
# sensor; the window, iterations and minimum agreeing fraction are AutoRIFT's defaults.
const DISPFILT_STRICT = AutoRIFT.GardnerFilter(; mad_scale = 3.0, agree_tolerance = 0.15)
const DISPFILT_OPTICAL = AutoRIFT.GardnerFilter(; mad_scale = 4.0, agree_tolerance = 0.15)

const DISPFILT_BY_PLATFORM = (
    "L4" => DISPFILT_OPTICAL, "L5" => DISPFILT_OPTICAL, "L7" => DISPFILT_OPTICAL,
    "L8" => DISPFILT_OPTICAL, "L9" => DISPFILT_OPTICAL, "AST" => DISPFILT_OPTICAL,
    "NISAR-L1" => DISPFILT_OPTICAL, "NISAR-L2" => DISPFILT_OPTICAL,
    "S2" => DISPFILT_STRICT, "S1-SLC" => DISPFILT_STRICT, "S1-BURST" => DISPFILT_STRICT,
)

"""
    dispfilt(platform) -> AutoRIFT.GardnerFilter

The displacement filter for `platform`. A platform without its own entry gets the optical filter, and a
warning saying so is written to stderr.
"""
function dispfilt(platform::String)
    for (p, f) in DISPFILT_BY_PLATFORM
        p == platform && return f
    end
    say("WARNING: no displacement filter is set for " * platform * "; using the optical one")
    return DISPFILT_OPTICAL
end
