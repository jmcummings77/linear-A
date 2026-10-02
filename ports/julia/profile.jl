# Collect actual Julia stack samples from the warmed computation loop.
include(joinpath(@__DIR__, "runner.jl"))
using Profile

function jsonstring(io, value)
    print(io, '"')
    for character in String(value)
        if character == '"'
            print(io, "\\\"")
        elseif character == '\\'
            print(io, "\\\\")
        elseif character == '\n'
            print(io, "\\n")
        elseif character == '\r'
            print(io, "\\r")
        elseif character == '\t'
            print(io, "\\t")
        elseif character < ' '
            print(io, "\\u", string(UInt32(character); base=16, pad=4))
        else
            print(io, character)
        end
    end
    print(io, '"')
end

function profile_main(args=ARGS)
    length(args) == 5 || throw(ArgumentError("usage: profile.jl OP SIZE ITERATIONS SEED OUTPUT"))
    name, a, b, iterations = preparebenchmark(args[1:4])
    runiterations(name, a, b, max(5, min(iterations, 100)))
    delay = 0.001
    Profile.init(n=10_000_000, delay=delay)
    Profile.clear()
    total = 0.0
    Profile.@profile total = runiterations(name, a, b, iterations)
    data, dictionary = Profile.retrieve(include_meta=false)
    stacks = Vector{String}[]
    current = String[]
    for pointer in data
        if pointer == 0
            if !isempty(current)
                push!(stacks, reverse(current))
                current = String[]
            end
        else
            for frame in get(dictionary, pointer, [])
                frame.from_c && continue
                push!(current, string(frame.func, " (", basename(String(frame.file)), ":", frame.line, ")"))
            end
        end
    end
    isempty(current) || push!(stacks, reverse(current))
    isempty(stacks) && throw(ArgumentError("no profile samples collected; increase ITERATIONS"))
    open(args[5], "w") do io
        print(io, "{\"implementation\":\"julia\",\"kind\":\"sampled\",\"chronological\":false,\"unit\":\"milliseconds\",\"checksum\":", total, ",\"stacks\":[")
        for (i, stack) in enumerate(stacks)
            i > 1 && print(io, ',')
            print(io, '[')
            for (j, frame) in enumerate(stack)
                j > 1 && print(io, ',')
                jsonstring(io, frame)
            end
            print(io, ']')
        end
        print(io, "],\"weights\":[")
        for i in eachindex(stacks)
            i > 1 && print(io, ',')
            print(io, delay * 1000)
        end
        print(io, "]}")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        profile_main()
    catch error
        showerror(stderr, error)
        println(stderr)
        exit(1)
    end
end
