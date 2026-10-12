# Trim probe: EzXML (libxml2) on a Sentinel-1 annotation and an orbit file, with the XPath queries and
# node walks SLCDatasets uses.
#
#   xml <annotation.xml> <orbit.EOF>

import EzXML

say(s::String) = (write(Core.stdout, s, "\n"); nothing)

function message(err)
    m = err isa ErrorException ? err.msg : nothing
    return m isa String ? m : string(typeof(err))
end

function text(root::EzXML.Node, xpath::String)
    n = findfirst(xpath, root)
    n === nothing && error(string("no ", xpath))
    return EzXML.nodecontent(n)
end

function (@main)(args::Vector{String})::Cint
    try
        doc = EzXML.readxml(args[1])
        r = EzXML.root(doc)
        samples = parse(Int, text(r, "imageAnnotation/imageInformation/numberOfSamples"))
        lines = parse(Int, text(r, "imageAnnotation/imageInformation/numberOfLines"))
        bursts = findall("swathTiming/burstList/burst", r)
        t0 = text(bursts[1], "azimuthTime")
        say(string("annotation samples=", samples, " lines=", lines, " bursts=", length(bursts), " first burst=", t0))
        firsts = [parse(Int, s) for s in split(text(bursts[1], "firstValidSample"))]
        say(string("  burst 1 firstValidSample count=", length(firsts), " max=", maximum(firsts)))
        eof = EzXML.root(EzXML.readxml(args[2]))
        osvs = findall("//OSV", eof)
        x = 0.0
        for osv in osvs, e in EzXML.eachelement(osv)
            EzXML.nodename(e) == "X" && (x += parse(Float64, EzXML.nodecontent(e)))
        end
        say(string("orbit state vectors=", length(osvs), " sum X=", x))
        return Cint(0)
    catch err
        say(string("error: ", message(err)))
        return Cint(1)
    end
end
