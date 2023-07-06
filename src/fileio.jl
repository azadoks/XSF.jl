using AtomsBase
using Unitful
using UnitfulAtomic
using Printf
using PeriodicTable: PeriodicTable

const LENGTH_UNIT = u"Å"
const FORCE_UNIT = u"Eh_au" / u"Å"
const KNOWN_PERIODIC_KEYWORDS = ("PRIMVEC", "CONVVEC", "PRIMCOORD", "CONVCOORD")

function iterate_xsf(itr::Base.EachLine)
    eof(itr.stream) && return itr.ondone()
    line = readline(itr.stream; keep=itr.keep)
    clean_line = strip(first(split(line, "#")))
    isempty(clean_line) && return iterate_xsf(itr)
    return string(clean_line)
end

function parse_boundary_conditions(line)
    occursin("ATOMS", line) && return [DirichletZero(), DirichletZero(), DirichletZero()]
    occursin("POLYMER", line) && return [Periodic(), DirichletZero(), DirichletZero()]
    occursin("SLAB", line) && return [Periodic(), Periodic(), DirichletZero()]
    occursin("CRYSTAL", line) && return [Periodic(), Periodic(), Periodic()]
    return error("Unknown structure type $(line)")
end

function parse_vec_block(T::Type{<:Real}, lines)
    return map(1:3) do _
        return parse.(T, split(iterate_xsf(lines))) .* LENGTH_UNIT
    end
end

function parse_coord_line(T::Type{<:Real}, line)
    words = split(line)
    number = parse(Int, words[1])
    atomic_symbol = Symbol(PeriodicTable.elements[number].symbol)
    position = parse.(T, words[2:4]) .* LENGTH_UNIT
    if length(words) == 7
        force = parse.(T, words[5:7]) .* FORCE_UNIT
        return Atom(; atomic_symbol, position, force=force)
    else
        return Atom(; atomic_symbol, position)
    end
end

function parse_coord_block(T::Type{<:Real}, lines)
    line = iterate_xsf(lines)
    n_atoms = parse(Int, first(split(line)))
    return map(1:n_atoms) do _
        return parse_coord_line(T, iterate_xsf(lines))
    end
end

function parse_xsf_block(T::Type{<:Real}, keyword, lines)
    if keyword in ("PRIMVEC", "CONVVEC")
        return parse_vec_block(T, lines)
    elseif keyword in ("PRIMCOORD", "CONVCOORD")
        return parse_coord_block(T, lines)
    else
        error("Unknown keyword $(keyword)")
    end
end

function parse_periodic_frame(T::Type{<:Real}, lines, bcs, previous_frame)
    should_parse = true
    i = 0
    io_pos = position(lines.stream)
    blocks = Dict()
    while should_parse && i <= 4
        line = iterate_xsf(lines)
        if isnothing(line)
            should_parse = false
        elseif (keyword = first(split(line))) in KNOWN_PERIODIC_KEYWORDS
            if haskey(blocks, keyword)
                should_parse = false
                seek(lines.stream, io_pos)
            else
                blocks[keyword] = parse_xsf_block(T, keyword, lines)
                io_pos = position(lines.stream)
            end
        else
            should_parse = false
            seek(lines.stream, io_pos)
        end
    end
    @assert haskey(blocks, "PRIMCOORD") "Found no PRIMCOORD block in the current frame"
    if !haskey(blocks, "PRIMVEC")
        if !isnothing(previous_frame)
            blocks["PRIMVEC"] = bounding_box(previous_frame)
        else
            error("Found no PRIMVEC block in the current frame and have no previous frame")
        end
    end
    return atomic_system(blocks["PRIMCOORD"], blocks["PRIMVEC"], bcs)
end

function parse_atoms_frame(T::Type{<:Real}, lines)
    line = iterate_xsf(lines)
    atoms = []
    while !isnothing(line) && !startswith(line, "ATOMS") && !startswith(line, "BEGIN")
        atom = parse_coord_line(T, line)
        push!(atoms, atom)
        line = iterate_xsf(lines)
    end
    return isolated_system(atoms)
end

# Single structure
function load_xsf(T::Type{<:Real}, file::Union{AbstractString,IOStream})
    lines = eachline(file; keep=false)
    line = iterate_xsf(lines)
    if occursin("ANIMSTEPS", line)
        n_frames = parse(Int, last(split(line)))
        line = iterate_xsf(lines)
    else
        n_frames = 1
    end
    bcs = parse_boundary_conditions(line)
    frames = AbstractSystem{3}[]
    for _ in 1:n_frames
        if count(Base.Fix2(isa, AtomsBase.Periodic), bcs) > 0 # CRYSTAL, SLAB, or POLYMER
            previous_frame = isempty(frames) ? nothing : last(frames)
            push!(frames, parse_periodic_frame(T, lines, bcs, previous_frame))
        else  # ATOMS
            push!(frames, parse_atoms_frame(T, lines))
        end
    end
    return length(frames) == 1 ? first(frames) : frames
end
load_xsf(file::Union{AbstractString,IOStream}) = load_xsf(Float64, file)

function write_system_type(io::IO, n_periodic_bcs::Int)
    types = Dict{Int,String}(1 => "POLYMER", 2 => "SLAB", 3 => "CRYSTAL")
    return println(io, types[n_periodic_bcs])
end

function write_atom(io::IO, atom)
    n = atomic_number(atom)
    x, y, z = ustrip(uconvert.(LENGTH_UNIT, position(atom)))
    if haskey(atom, :force)
        fx, fy, fz = ustrip(uconvert.(FORCE_UNIT, get(atom, :force, nothing)))
        @printf io "%3d %20.14f %20.14f %20.14f %20.14f %20.14f %20.14f\n" n x y z fx fy fz
    else
        @printf io "%3d %20.14f %20.14f %20.14f\n" n x y z
    end
end
function write_atoms(io::IO, system::AbstractSystem; header="")
    !isempty(header) && println(io, header)
    return map(Base.Fix1(write_atom, io), system[:])
end

function write_bounding_box(io::IO, system::AbstractSystem; header="")
    !isempty(header) && println(io, header)
    for i in 1:3
        x, y, z = ustrip(uconvert.(LENGTH_UNIT, bounding_box(system)[i]))
        @printf io "%20.14f %20.14f %20.14f\n" x y z
    end
end

function save_xsf(io::IO, system::AbstractSystem{3}; frame="")
    bcs = boundary_conditions(system)
    n_periodic_bcs = count(Base.Fix2(isa, Periodic), bcs)
    if n_periodic_bcs == 0
        write_atoms(io, system; header="ATOMS $(frame)")
    else
        n_atoms = length(system)
        if isempty(frame) || isnothing(frame) || frame == 1
            write_system_type(io, n_periodic_bcs)
        end
        write_bounding_box(io, system; header="PRIMVEC $(frame)")
        write_bounding_box(io, system; header="CONVVEC $(frame)")
        write_atoms(io, system; header="PRIMCOORD $(frame)\n$(n_atoms) 1")
    end
end

function save_xsf(io::IO, frames)
    println(io, "ANIMSTEPS $(length(frames))")
    for i in eachindex(frames)
        save_xsf(io, frames[i]; frame=i)
    end
end

function save_xsf(file::AbstractString, frames)
    open(file, "w") do io
        save_xsf(io, frames)
    end
end
