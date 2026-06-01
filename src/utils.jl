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
