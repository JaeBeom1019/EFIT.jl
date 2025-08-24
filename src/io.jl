using Printf
import Dates
import IMASdd
import CoordinateConventions

mutable struct GEQDSKFile
    file::String                    # Source file
    time::Float64                   # Time of the equilibrium reconstruction, in seconds
    nw::Int                         # Number of horizontal R grid points
    nh::Int                         # Number of vertical Z grid points
    r::AbstractRange{Float64}       # R grid points
    z::AbstractRange{Float64}       # Z grid points
    rdim::Float64                   # Horizontal dimension in meter of computational box
    zdim::Float64                   # Vertical dimension in meter of computational box
    rleft::Float64                  # Minimum R in meter of rectangular computational box
    zmid::Float64                   # Z of center of computational box in meter
    nbbbs::Int                      # Number of boundary points
    rbbbs::Vector{Float64}          # R of boundary points in meter
    zbbbs::Vector{Float64}          # Z of boundary points in meter
    limitr::Int                     # Number of limiter points
    rlim::Vector{Float64}           # R of surrounding limiter contour in meter
    zlim::Vector{Float64}           # Z of surrounding limiter contour in meter
    rcentr::Float64                 # R in meter of vacuum toroidal magnetic field BCENTR
    bcentr::Float64                 # Vacuum toroidal magnetic field in Tesla at RCENTR
    rmaxis::Float64                 # R of magnetic axis in meter
    zmaxis::Float64                 # Z of magnetic axis in meter
    simag::Float64                  # Poloidal flux at magnetic axis in Weber/rad
    sibry::Float64                  # Poloidal flux at the plasma boundary in Weber/rad
    psi::AbstractRange{Float64}     # Poloidal flux points
    current::Float64                # Plasma current in Ampere
    fpol::Vector{Float64}           # Poloidal current function in m-T, F = RBT on flux grid
    pres::Vector{Float64}           # Plasma pressure in nt / m2 on uniform flux grid
    ffprim::Vector{Float64}         # FF'(psi) in (mT)2 / (Weber/rad) on uniform fl ux grid
    pprime::Vector{Float64}         # P'(psi) in (nt/m2) / (Weber/rad) on uniform flux grid
    qpsi::Vector{Float64}           # q values on uniform flu x grid from axis to boundary
    psirz::Matrix{Float64}          # Poloidal flux in Weber/rad on the rectangular grid points
    rhovn::Vector{Float64}          # Square root of toroidal flux, normalized
end

# Define hash function using fieldnames and getfield
function Base.hash(x::GEQDSKFile, h::UInt)
    for field in fieldnames(GEQDSKFile)
        h = hash(getfield(x, field), h)
    end
    return h
end
# Define == using fieldnames and getfield
Base.:(==)(a::GEQDSKFile, b::GEQDSKFile) = all(field -> getfield(a, field) == getfield(b, field), fieldnames(GEQDSKFile))
# Define isequal using fieldnames and getfield
Base.isequal(a::GEQDSKFile, b::GEQDSKFile) = all(field -> isequal(getfield(a, field), getfield(b, field)), fieldnames(GEQDSKFile))

function Base.show(io::IO, g::GEQDSKFile)
    print(io,"GEQDSKFile: \"",g.file,"\"")
end

Base.broadcastable(g::GEQDSKFile) = (g,)

function file_numbers(f)
    Channel(ctype=Float64) do c
        while ~eof(f)
            line = readline(f)
            pattern = r"[+-]?\d*[\.]?\d+(?:[Ee][+-]?\d+)?"
            for p in eachmatch(pattern,line)
                push!(c,parse(Float64,p.match))
            end
        end
    end
end

function read_array(t,n)
    data = zeros(n)
    for i=1:n
        try
            data[i] = take!(t)
        catch e
            if isa(e, InvalidStateException)
                # InvalidStateException when Channel is closed
            else
                rethrow(e)
            end
        end
    end
    return data
end

function read_array2d(t,n,m)
    data = zeros(n,m)
    for i=1:n, j=1:m
        try
            data[i,j] = take!(t)
        catch e
            if isa(e, InvalidStateException)
                # InvalidStateException when Channel is closed
            else
                rethrow(e)
            end
        end
    end
    return data'
end



function parse_gfile_header(headerline::String; set_time=nothing)
    parts = split(headerline)

    # --- 1. 헤더 라인에서 모든 숫자와 그 원본 문자열을 추출합니다. ---
    all_numeric_parts = []
    for part in parts
        # "ms"나 "s"같은 단위는 제거하고 숫자 부분만 파싱 시도
        numeric_str = replace(part, r"ms$|s$"i => "")
        num = tryparse(Float64, numeric_str)
        if !isnothing(num)
            # (파싱된 숫자, 원본 문자열) 형태로 저장
            push!(all_numeric_parts, (value=num, original=part))
        end
    end

    if length(all_numeric_parts) < 3
        error("헤더에서 필요한 숫자(idum, nw, nh) 3개 이상을 찾지 못했습니다: $headerline")
    end

    # --- 2. 마지막 세 숫자를 idum, nw, nh로 할당합니다. ---
    nh = Int(all_numeric_parts[end].value)
    nw = Int(all_numeric_parts[end-1].value)
    idum = Int(all_numeric_parts[end-2].value)

    # --- 3. 시간 값을 결정합니다. ---
    time_val = 0.0
    if !isnothing(set_time)
        time_val = set_time
    else
        time_found_with_unit = false
        # 1순위: 'ms' 또는 's' 단위가 명시된 값을 찾습니다.
        for part in parts
            if endswith(lowercase(part), "ms")
                num_str = replace(part, r"ms$"i => "")
                parsed_num = tryparse(Float64, num_str)
                if !isnothing(parsed_num)
                    time_val = parsed_num / 1000.0  # ms -> s
                    time_found_with_unit = true
                    break
                end
            elseif endswith(lowercase(part), "s")
                num_str = replace(part, r"s$"i => "")
                parsed_num = tryparse(Float64, num_str)
                if !isnothing(parsed_num)
                    time_val = parsed_num  # s -> s
                    time_found_with_unit = true
                    break
                end
            end
        end

        # 2순위: 단위를 찾지 못했다면, 나머지 숫자 중에서 추론합니다.
        if !time_found_with_unit
            # idum, nw, nh를 제외한 숫자 후보들
            time_candidates = all_numeric_parts[1:end-3]
            if !isempty(time_candidates)
                # 가장 마지막 후보를 시간으로 간주
                last_candidate = time_candidates[end]
                if occursin('.', last_candidate.original)
                    # 원본에 소수점이 있으면 '초(s)'로 간주
                    time_val = last_candidate.value
                else
                    # 정수이면 '밀리초(ms)'로 간주하고 변환
                    time_val = last_candidate.value / 1000.0
                end
            end
        end
    end

    return idum, time_val, nw, nh
end

function readg(gfile; set_time=nothing)

    isfile(gfile) || error("$(gfile) does not exist")

    f = open(gfile)

    desc = readline(f)
    idum, time, nw, nh = parse_gfile_header(desc; set_time)

    token = file_numbers(f)

    rdim   = take!(token)
    zdim   = take!(token)
    rcentr = take!(token)
    rleft  = take!(token)
    zmid   = take!(token)

    rmaxis = take!(token)
    zmaxis = take!(token)
    simag  = take!(token)
    sibry  = take!(token)
    bcentr = take!(token)

    current = take!(token)
    simag  = take!(token)
    xdum    = take!(token)
    rmaxis  = take!(token)
    xdum    = take!(token)

    zmaxis = take!(token)
    xdum   = take!(token)
    sibry = take!(token)
    xdum   = take!(token)
    xdum   = take!(token)

    fpol = read_array(token,nw)
    pres = read_array(token,nw)
    ffprim = read_array(token,nw)
    pprime = read_array(token,nw)
    psirz = read_array2d(token,nh,nw)
    qpsi = read_array(token,nw)

    nbbbs = Int(take!(token))
    limitr = Int(take!(token))

    if nbbbs > 0
        rbbbs = zeros(nbbbs)
        zbbbs = zeros(nbbbs)
        for i=1:nbbbs
            rbbbs[i] = take!(token)
            zbbbs[i] = take!(token)
        end
    else
        rbbbs = [0.0]
        zbbbs = [0.0]
    end

    if limitr > 0
        rlim = zeros(limitr)
        zlim = zeros(limitr)
        for i=1:limitr
            rlim[i] = take!(token)
            zlim[i] = take!(token)
        end
    else
        rlim = [0.0]
        zlim = [0.0]
    end

    local rhovn
    try
        for i ∈ 1:3
            if !eof(f)
                xdum = take!(token)
            end
        end
        if !eof(f)
            rhovn = read_array(token,nw)
        else
            rhovn = zeros(nw)
        end
        if !eof(f)
            xdum = take!(token)
        end
    catch e
        if isa(e, InvalidStateException)
            # InvalidStateException when Channel is closed
            rhovn = zeros(nw)
        else
            rethrow(e)
        end
    end

    close(f)
    close(token)

    r = range(rleft, rleft + rdim, length=nw)
    z = range(zmid - 0.5*zdim, zmid + 0.5*zdim, length=nh)
    psi = range(simag, sibry, length=nw)

    g = GEQDSKFile(
        splitpath(gfile)[end],time,nw,nh,r,z,rdim,zdim,rleft,zmid,nbbbs,rbbbs,zbbbs,
        limitr,rlim,zlim,rcentr,bcentr,rmaxis,zmaxis,simag,sibry,psi,current,fpol,pres,
        ffprim,pprime,qpsi,psirz,rhovn
    )

    return g
end

function write_vector_data_in_chunks(io::IOStream, v::Vector{Float64})
    # Iterate over the vector in chunks of 5 elements
    for chunk in Iterators.partition(v, 5)
        # Format each number with %16.9E
        formatted_numbers = [@sprintf("%16.9E", num) for num in chunk]
        # Join the formatted numbers into a single string without whitespaces
        line = join(formatted_numbers)
        # Write the line to the io stream
        println(io, line)
    end
end

# Function to write the GEQDSKFile struct to a G-file with explicit error handling
function writeg(g::GEQDSKFile, filename::String;
                desc::String="EFIT.jl   $(Dates.format(Dates.today(), "dd/mm/yyyy"))",
                shot::String="#000000",
                time::String= @sprintf("%dms", 1e3 * g.time))
    if isdir(filename)
        @error("Error: A directory with the name '$filename' already exists.")
        return false
    elseif isfile(filename)
        @warn("Warning: A file with the name '$filename' already exists. It will be overwritten.")
    end

    # Note: GEQDSKFile.time is converted to ms (to follow readg's notation)
    # Remove all newline characters in description
    clean_desc = replace(desc, "\n" => "")
    description = join([clean_desc, shot, time], "   ")
    if length(description) > 48
        error("""
            Description too long (length = $(length(description)) > max 48).
            Current description: '$description'
            Please shorten one or more of the following kwargs:
                desc = "$desc"
                shot = "$shot"
                time = "$time"
            """)
    end

    try
        # Open the target file for writing
        open(filename, "w") do f
            # Write header
            @printf(f,"%-48s%4d%4d%4d\n",description, 0, g.nw, g.nh)

            @printf(f,"%16.9E%16.9E%16.9E%16.9E%16.9E\n", g.rdim, g.zdim, g.rcentr, g.rleft, g.zmid)
            @printf(f,"%16.9E%16.9E%16.9E%16.9E%16.9E\n", g.rmaxis, g.zmaxis, g.simag, g.sibry, g.bcentr)
            @printf(f,"%16.9E%16.9E%16.9E%16.9E%16.9E\n", g.current, g.simag, 0.0, g.rmaxis, 0.0)
            @printf(f,"%16.9E%16.9E%16.9E%16.9E%16.9E\n", g.zmaxis, 0.0, g.sibry, 0.0, 0.0)

            write_vector_data_in_chunks(f, g.fpol)
            write_vector_data_in_chunks(f, g.pres)
            write_vector_data_in_chunks(f, g.ffprim)
            write_vector_data_in_chunks(f, g.pprime)

            write_vector_data_in_chunks(f, vec(g.psirz))

            write_vector_data_in_chunks(f, g.qpsi)

            @printf(f,"%d %d\n", g.nbbbs, g.limitr)

            bbbs_rz = Vector(vec(hcat(g.rbbbs, g.zbbbs)'))
            write_vector_data_in_chunks(f, vec(bbbs_rz))

            lim_rz = Vector(vec(hcat(g.rlim, g.zlim)'))
            write_vector_data_in_chunks(f, vec(lim_rz))
        end
        @info "Successfully wrote GEQDSKFile to '$filename'."
        return true
    catch e
        @error """An error occurred while writing to file '$filename':
             $(sprint(showerror, e))"""
        return false
    end
end
