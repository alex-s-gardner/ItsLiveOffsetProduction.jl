# The juliac trim verifier's errors in a build log, counted by the package whose code holds each
# unresolved call (the first stack frame outside Base and Core), with one example statement each.
#
#   julia tools/trim/verifier_report.jl build.log [max_rows]

function report(io::IO, log::AbstractString; max_rows::Int = 25)
    text = replace(read(log, String), r"\e\[[0-9;]*m" => "")
    blocks = split(text, r"\n(?=Verifier (?:error|warning) #)")
    counts = Dict{Tuple{String,String},Int}()
    example = Dict{Tuple{String,String},String}()
    nerr = 0
    for b in blocks
        startswith(b, "Verifier error") || continue
        nerr += 1
        frames = [(m[1], m[2], m[3]) for m in eachmatch(r"\n\s*\[\d+\] (.*)\n\s+@ (\S+) (\S+)", b)]
        site = something(findfirst(f -> !(f[2] in ("Base", "Core", "Base.Broadcast")), frames), 1)
        key = isempty(frames) ? ("?", "?") :
              (frames[site][2], string(first(split(frames[site][1], '(')), " ", basename(frames[site][3])))
        counts[key] = get(counts, key, 0) + 1
        statement = first(split(b, '\n'))
        get!(example, key, first(replace(statement, r"^Verifier error #\d+: " => ""), 160))
    end
    summary = match(r"Trim verify finished with (\d+) errors?, (\d+) warnings?", text)
    if summary === nothing
        # A clean verification prints no summary line.
        println(io, "verifier: $nerr errors")
    else
        println(io, "verifier: ", summary[1], " errors, ", summary[2], " warnings")
    end
    for (key, n) in first(sort!(collect(counts); by = last, rev = true), max_rows)
        println(io, lpad(n, 5), "  ", key[1], "  ", key[2])
        println(io, "         ", example[key])
    end
end

report(stdout, ARGS[1]; max_rows = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 25)
