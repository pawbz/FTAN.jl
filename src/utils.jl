function _sample_scheme_colors(name::AbstractString, n::Int)
    if n <= 0
        return String[]
    end

    candidates = (
        Symbol(name),
        Symbol(lowercase(name)),
        Symbol(replace(name, " " => "_")),
        Symbol(lowercase(replace(name, " " => "_"))),
    )

    cs = nothing
    for key in candidates
        if haskey(ColorSchemes.colorschemes, key)
            cs = ColorSchemes.colorschemes[key]
            break
        end
    end
    cs === nothing && (cs = ColorSchemes.colorschemes[:viridis])

    ts = n == 1 ? [0.5] : collect(range(0.0, 1.0, length=n))
    return [
        let c = get(cs, t)
            "rgb($(round(Int, 255 * red(c))),$(round(Int, 255 * green(c))),$(round(Int, 255 * blue(c))))"
        end
        for t in ts
    ]
end

function _periods_from_nyquist(dt::Float64; period_min::Union{Float64, Nothing}=nothing, dT::Float64=0.5, period_max::Float64=60.0)
    dT > 0.0 || throw(ArgumentError("dT must be positive"))
    period_max > 0.0 || throw(ArgumentError("period_max must be positive"))

    period_min = period_min === nothing ? 2.0 * dt : period_min
    period_max >= period_min || throw(ArgumentError("period_max=$(period_max) must be >= Nyquist period $(period_min)"))

    periods = collect(exp10.(range(log10(period_min), log10(period_max), length=400)))
    isempty(periods) && (periods = [period_min])
    return periods
end

function _column_normalise(X)
    Xf = Float64.(X)
    return mapslices(x -> begin
        n = norm(x)
        n > 0 ? x ./ n : x
    end, Xf; dims=1)
end

# ╔═╡ a1000002-0000-0000-0000-000000000001
function _vector_normalise(x)
    xf = Float64.(vec(x))
    n = norm(xf)
    return n > 0 ? xf ./ n : xf
end

# ╔═╡ a1000003-0000-0000-0000-000000000001
function _ncc(a, b)
    a0 = Float64.(vec(a)) .- mean(a)
    b0 = Float64.(vec(b)) .- mean(b)
    return dot(a0, b0) / ((norm(a0) * norm(b0)) + 1e-8)
end

"""
    wavelength_valid_period(period, distance; wavelength_ref_velocity=nothing,
                            wavelength_fraction=nothing)

Return whether `period` passes the optional distance-dependent wavelength filter.
`wavelength_fraction` is the minimum number of wavelengths (at
`wavelength_ref_velocity`) that must fit within `distance` for the period to
be considered valid. When filtering is requested, keep periods satisfying
`wavelength_fraction * wavelength_ref_velocity * period < distance` (i.e.
`distance > wavelength_fraction` wavelengths).
"""
function wavelength_valid_period(period::Real, distance::Real;
                                 wavelength_ref_velocity::Union{Nothing,Real}=nothing,
                                 wavelength_fraction::Union{Nothing,Real}=nothing)
    isnothing(wavelength_ref_velocity) && isnothing(wavelength_fraction) && return true
    (isnothing(wavelength_ref_velocity) || isnothing(wavelength_fraction)) && return false

    period = Float64(period)
    distance = Float64(distance)
    ref_velocity = Float64(wavelength_ref_velocity)
    fraction = Float64(wavelength_fraction)
    all(isfinite, (period, distance, ref_velocity, fraction)) || return false
    (period > 0.0 && distance > 0.0 && ref_velocity > 0.0 && fraction > 0.0) || return false
    return fraction * ref_velocity * period < distance
end

function _wavelength_valid_indices(periods::AbstractVector{<:Real}, distance::Real;
                                   wavelength_ref_velocity::Union{Nothing,Real}=nothing,
                                   wavelength_fraction::Union{Nothing,Real}=nothing)
    return findall(period -> wavelength_valid_period(period, distance;
                                                     wavelength_ref_velocity=wavelength_ref_velocity,
                                                     wavelength_fraction=wavelength_fraction),
                   periods)
end

const EARTH_RADIUS_KM = 6371.0

"""
    haversine_distance_km(lat1, lon1, lat2, lon2) -> Float64

Great-circle distance [km] between two lat/lon points [degrees] using the
haversine formula and a spherical Earth radius of 6371 km.
"""
function haversine_distance_km(lat1::Real, lon1::Real, lat2::Real, lon2::Real)
    φ1, φ2 = deg2rad(lat1), deg2rad(lat2)
    Δφ, Δλ = deg2rad(lat2 - lat1), deg2rad(lon2 - lon1)
    a = sin(Δφ / 2)^2 + cos(φ1) * cos(φ2) * sin(Δλ / 2)^2
    return 2 * EARTH_RADIUS_KM * asin(min(1.0, sqrt(a)))
end

"""
    linear_array_order(station_coords::Dict{String,Tuple{Float64,Float64}}) -> Vector{String}

Order station codes along a roughly-linear array's principal axis: project
centered (lat, lon) positions onto the first principal component of their
2x2 covariance matrix, then sort by that projected coordinate.
"""
function linear_array_order(station_coords::Dict{String,Tuple{Float64,Float64}})
    codes = collect(keys(station_coords))
    isempty(codes) && return String[]
    length(codes) == 1 && return codes

    lats = [station_coords[c][1] for c in codes]
    lons = [station_coords[c][2] for c in codes]
    lat0, lon0 = mean(lats), mean(lons)
    X = hcat(lats .- lat0, lons .- lon0)  # (n x 2)

    cov = X' * X
    evals, evecs = eigen(Symmetric(cov))
    axis = evecs[:, argmax(evals)]  # principal direction (2-vector)

    proj = X * axis
    order = sortperm(proj)
    return codes[order]
end
