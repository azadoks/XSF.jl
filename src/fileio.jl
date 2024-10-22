const KNOWN_PERIODIC_KEYWORDS = ("PRIMVEC", "CONVVEC", "PRIMCOORD", "CONVCOORD")

# Count the number of boundary conditions which are AtomsBase.Periodic
count_periodic_bcs(system::AbstractSystem) = count(periodicity(system))

function check_system_properties(system::AbstractSystem)
    system_keys = keys(system)
    for key in system_keys
        if !in(key, (:bounding_box, :periodicity, ))
            @warn "Ignoring unsupported property $(key)"
        end
    end
end

function check_atom_properties(atom)
    atom_keys = keys(atom)
    for key in atom_keys
        if !in(key, (:species, :mass, :position, :force,))
            @warn "Ignoring unsupported atomic property $(key)"
        end
    end
end

function check_atomic_mass(atom)
    if haskey(PeriodicTable.elements, atomic_symbol(atom))
        if mass(atom) != PeriodicTable.elements[atomic_symbol(atom)].atomic_mass
            @warn "Atom mass in XSF cannot be mutated."
        end
    end
end

function check_atomic_symbol(atom)
    if atomic_symbol(atom) != Symbol(element(species(atom)).symbol)
        @warn "Atom atomic_symbol in XSF must agree with atomic number"
    end
end

function check_system(system::AbstractSystem)
    check_system_properties(system)
    check_atom_properties.(system)
    check_atomic_mass.(system)
    check_atomic_symbol.(system)
    return nothing
end

# Custom version of `iterate(::Base.EachLine)` which cleans lines for parsing
# and skips comment / empty lines
function iterate_xsf(itr::Base.EachLine)
    eof(itr.stream) && return itr.ondone()
    line = readline(itr.stream; keep=itr.keep)
    clean_line = strip(first(split(line, "#")))
    isempty(clean_line) && return iterate_xsf(itr)
    return string(clean_line)
end

# Get boundary conditions in the XSF convention from the system type keyword
# NOTE: the XSF documentation talks about a keyword `MOLECULE`, but it isn't
# mentioned in the ase.io parser nor in any of the examples, so it is not
# implemented here
function parse_boundary_conditions(line)
    occursin("ATOMS", line)   && return (false, false, false)
    occursin("POLYMER", line) && return (true,  false, false)
    occursin("SLAB", line)    && return (true,   true, false)
    occursin("CRYSTAL", line) && return (true,   true,  true)
    return error("Unknown structure type $(line)")
end

# Parse blocks like PRIMVEC and CONVVEC
function parse_vec_block(T::Type{<:Real}, lines)
    return map(1:3) do _
        return parse.(T, split(iterate_xsf(lines))) .* LENGTH_UNIT
    end
end

# Parse an atomic position line, which is formed as either:
# atomic_number x y z
# or
# atomic_number x y z Fx Fy Fz
function parse_coord_line(T::Type{<:Real}, line)
    words = split(line)
    species = ChemicalSpecies(parse(Int, words[1]))
    position = parse.(T, words[2:4]) .* LENGTH_UNIT
    if length(words) == 7
        force = parse.(T, words[5:7]) .* FORCE_UNIT
        return Atom(; species, position, force=force)
    else
        return Atom(; species, position)
    end
end

# Parse an entire block of atomic coordinates in a periodic system,
# i.e. when the number of atoms is provided explicitly
function parse_coord_block(T::Type{<:Real}, lines)
    line = iterate_xsf(lines)
    n_atoms = parse(Int, first(split(line)))
    return map(1:n_atoms) do _
        return parse_coord_line(T, iterate_xsf(lines))
    end
end

# Parse an arbitrary block from a periodic file
# Supported keywords are PRIMVEC, CONVVEC, PRIMCOORD, and CONVCOORD, i.e.
# DATAGRID_[2,3]D blocks are not (yet) supported
function parse_xsf_block(T::Type{<:Real}, keyword, lines)
    if keyword in ("PRIMVEC", "CONVVEC")
        return parse_vec_block(T, lines)
    elseif keyword in ("PRIMCOORD", "CONVCOORD")
        return parse_coord_block(T, lines)
    else
        error("Unknown keyword $(keyword)")
    end
end

# Parse a frame (ANIMSTEP in XSF lingo) from a periodic file
# The first frame will have at least PRIMVEC and PRIMCOORD and possibly CONVVEC
# and CONVCOORD blocks
# Following frames will have at least PRIMCOORD (PRIMVEC and others are optional)
# If a frame doesn't have PRIMVEC, the bounding box (PRIMVEC) from the previous frame needs
# to be brought forward
function parse_periodic_frame(T::Type{<:Real}, lines, pbcs, previous_frame)
    should_parse = true
    i = 0
    io_pos = position(lines.stream)
    blocks = Dict()
    # Maximum 4 blocks (PRIMVEC, PRIMCOORD, CONVVEC, CONVCOORD)
    while should_parse && i <= 4
        line = iterate_xsf(lines)
        if isnothing(line)  # We've reached the end of the file
            should_parse = false
        elseif (keyword = first(split(line))) in KNOWN_PERIODIC_KEYWORDS
            if haskey(blocks, keyword)
                # If we've already seen this keyword, we're in the next block
                should_parse = false
                seek(lines.stream, io_pos)  # Go back to the keyword line
            else
                blocks[keyword] = parse_xsf_block(T, keyword, lines)
                io_pos = position(lines.stream)
            end
        else  # We've reached the end of the frame
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
    return FlexibleSystem(blocks["PRIMCOORD"], blocks["PRIMVEC"], pbcs)
end

# Parse a frame (ANIMSTEP in XSF lingo) from a non-periodic file (ATOMS)
function parse_atoms_frame(T::Type{<:Real}, lines)
    line = iterate_xsf(lines)
    atoms = []
    # Go until the end of the file or until we see another frame (ATOMS) or data grid (BEGIN)
    while !isnothing(line) && !startswith(line, "ATOMS") && !startswith(line, "BEGIN")
        atom = parse_coord_line(T, line)
        push!(atoms, atom)
        line = iterate_xsf(lines)
    end
    return FlexibleSystem(atoms, IsolatedCell(3))
end

# Load all the frames from an XSF file
function load_xsf(T::Type{<:Real}, file::Union{AbstractString,IOStream})
    lines = eachline(file; keep=false)
    # The first line should be "ANIMSTEPS [N]" for a trajectory (in which case, the system
    # type is in the second line) or one of the system type keywords (ATOMS, POLYMER, SLAB,
    # CRYSTAL), in which case the number of frames is 1
    line = iterate_xsf(lines)
    if occursin("ANIMSTEPS", line)
        n_frames = parse(Int, last(split(line)))
        line = iterate_xsf(lines)
    else
        n_frames = 1
    end
    pbcs = parse_boundary_conditions(line)
    frames = AbstractSystem{3}[]
    for _ in 1:n_frames
        # Check how many boundary conditions are Periodic to determine whether we have
        # ATOMS or one of POLYMER, SLAB, CRYSTAL
        if !any(pbcs)  # no periodic BCs
            push!(frames, parse_atoms_frame(T, lines))
        else
            # We need to pass the previous frame because if the unit cell is fixed,
            # it is written only in the first frame and we need to pass it through
            # frame-by-frame
            previous_frame = isempty(frames) ? nothing : last(frames)
            push!(frames, parse_periodic_frame(T, lines, pbcs, previous_frame))
        end
    end
    return frames
end
# Set a default floating-point type of Float64
load_xsf(file::Union{AbstractString,IOStream}) = load_xsf(Float64, file)

function write_system_type(io::IO, n_periodic_bcs::Int)
    types = Dict{Int,String}(1 => "POLYMER", 2 => "SLAB", 3 => "CRYSTAL")
    println(io, types[n_periodic_bcs])
    return nothing
end

function write_atom(io::IO, atom)
    n = atomic_number(atom)
    x, y, z = ustrip.(uconvert.(LENGTH_UNIT, position(atom)))
    if haskey(atom, :force)
        fx, fy, fz = ustrip.(uconvert.(FORCE_UNIT, get(atom, :force, nothing)))
        println(io, "$(n) $(x) $(y) $(z) $(fx) $(fy) $(fz)")
    else
        println(io, "$(n) $(x) $(y) $(z)")
    end
    return nothing
end

# Write the atomic numbers, positions[, forces] of an ATOMS block or [PRIM,CONV]COORD block
function write_atoms(io::IO, system::AbstractSystem; header="")
    !isempty(header) && println(io, strip(header))
    return map(Base.Fix1(write_atom, io), system[:])
end

# Write a bounding box block (optionally with a header, i.e. PRIMVEC or CONVVEC)
function write_bounding_box(io::IO, system::AbstractSystem; header="")
    !isempty(header) && println(io, strip(header))
    for i in 1:3
        x, y, z = ustrip.(uconvert.(LENGTH_UNIT, bounding_box(system)[i]))
        println(io, "$(x) $(y) $(z)")
    end
    return nothing
end

# Write the PRIMVEC, CONVVEC, and PRIMCOORD blocks which make up a frame of a
# periodic trajectory
function write_periodic_frame(io::IO, system::AbstractSystem; frame="")
    n_atoms = length(system)
    write_bounding_box(io, system; header="PRIMVEC $(frame)")
    write_bounding_box(io, system; header="CONVVEC $(frame)")
    write_atoms(io, system; header="PRIMCOORD $(frame)\n$(n_atoms) 1")
    return nothing
end

# Write an ATOMS frame or periodic frame depending on the periodicity of the system
# This function is used for writing the frames of animated files (not for single
# structures).
function write_frame(io::IO, system::AbstractSystem; frame="")
    if count_periodic_bcs(system) == 0
        write_atoms(io, system; header="ATOMS $(frame)")
    else
        write_periodic_frame(io, system; frame)
    end
    return nothing
end

function save_xsf(io::IO, frames::AbstractVector{<:AbstractSystem})
    # Check the frames and warn about unsupported properties.
    check_system.(frames)
    # Make sure all structures have the same boundary conditions and get those boundary
    # conditions.
    pbcs = only(unique(periodicity.(frames)))
    # Count the number of periodic boundary conditions.
    n_periodic_bcs = count(pbcs)
    # Write the animation line if more than one frame is provided.
    is_animated = length(frames) > 1
    is_animated && println(io, "ANIMSTEPS $(length(frames))")
    # For periodic systems, write the system type line.
    n_periodic_bcs > 0 && write_system_type(io, n_periodic_bcs)
    # Write the structural and force information for each frame.
    for i in eachindex(frames)
        write_frame(io, frames[i], frame=is_animated ? i : "")
    end
    return nothing
end

function save_xsf(file::AbstractString, frames::AbstractVector{<:AbstractSystem})
    open(file, "w") do io
        return save_xsf(io, frames)
    end
end

function save_xsf(file_or_io::Union{IO, AbstractString}, frame::AbstractSystem)
    return save_xsf(file_or_io, [frame])
end
