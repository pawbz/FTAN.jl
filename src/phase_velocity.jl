begin
    DEFAULT_PHASE_BRANCH_NUMBERS = collect(-3:3)
    function _fill_phase_velocity_branches!(phase_velocity_branches::AbstractMatrix{Float64}, frequencies::AbstractVector{Float64}, arrival_times::AbstractVector{Float64}, measured_phases::AbstractVector{Float64}, distance::Float64, phase_branch_numbers::AbstractVector{Int}, phvel_source_phase::Float64; phase_plausible_range::Tuple{Float64,Float64}=(0.5, 20.0))
        fill!(phase_velocity_branches, NaN); c_min, c_max = phase_plausible_range
        for i in eachindex(frequencies)
            t_g = arrival_times[i]; phi = measured_phases[i]
            (isfinite(t_g) && t_g > 0.0 && isfinite(phi)) || continue
            omega = 2π * frequencies[i]
            for (ib, branch) in enumerate(phase_branch_numbers)
                denom = omega * t_g - phi - phvel_source_phase + 2π * branch
                abs(denom) > 1e-6 || continue
                c = omega * distance / denom
                phase_velocity_branches[i, ib] = (c_min <= c <= c_max) ? c : NaN
            end
        end
        return phase_velocity_branches
    end
    function _nearest_finite_branch(branch_values::AbstractVector{Float64}, target::Float64)
        best_j = 0; best_score = Inf
        for j in eachindex(branch_values)
            c = branch_values[j]; isfinite(c) || continue
            score = abs(c - target)
            score < best_score && (best_score = score; best_j = j)
        end
        return best_j
    end
    function _nearest_anchor_branch(branch_values::AbstractVector{Float64}, target::Float64, phase_velocity_range::Tuple{Float64,Float64})
        vmin, vmax = phase_velocity_range; best_j = 0; best_score = Inf
        for j in eachindex(branch_values)
            c = branch_values[j]; (isfinite(c) && vmin <= c <= vmax) || continue
            score = abs(c - target)
            score < best_score && (best_score = score; best_j = j)
        end
        return best_j == 0 ? _nearest_finite_branch(branch_values, target) : best_j
    end
    function compute_group_velocity_from_phase(periods::AbstractVector{Float64}, phase_velocities::AbstractVector{Float64})
        u_pred = fill(NaN, length(periods))
        valid = findall(i -> isfinite(periods[i]) && isfinite(phase_velocities[i]) && periods[i] > 0.0 && phase_velocities[i] > 0.0, eachindex(periods))
        length(valid) >= 2 || return u_pred
        sorted = valid[sortperm(periods[valid])]; n = length(sorted)
        for (pos, i) in enumerate(sorted)
            j1, j2 = pos == 1 ? (sorted[1], sorted[2]) : pos == n ? (sorted[n - 1], sorted[n]) : (sorted[pos - 1], sorted[pos + 1])
            dT = periods[j2] - periods[j1]; abs(dT) > 0.0 || continue
            dcdT = (phase_velocities[j2] - phase_velocities[j1]) / dT; c = phase_velocities[i]
            denom = 1.0 + (periods[i] / c) * dcdT; abs(denom) > 1e-6 || continue
            u_pred[i] = c / denom
        end
        return u_pred
    end
    compute_group_velocity_from_phase(res::MFTResult) = compute_group_velocity_from_phase(res.periods, res.phase_velocities)
    function resolve_phase_velocity_cycles!(phase_velocities::Vector{Float64}, selected_phase_branches::Vector{Int}, phase_suspect::AbstractVector{Bool}, u_predicted_from_phase::Vector{Float64}, phase_velocity_branches::Matrix{Float64}, periods::Vector{Float64}, quality_factors::Vector{Float64}, phase_branch_numbers::Vector{Int}; min_anchor_quality::Float64=3.0, phase_anchor_velocity::Float64=3.3, phase_velocity_range::Tuple{Float64,Float64}=(2.5, 4.5), phase_smoothness_jump::Float64=0.05)
        fill!(phase_velocities, NaN); fill!(selected_phase_branches, 0); fill!(phase_suspect, false); fill!(u_predicted_from_phase, NaN)
        valid = findall(i -> any(isfinite, @view phase_velocity_branches[i, :]), eachindex(periods)); isempty(valid) && return phase_velocities
        anchor_pool = [i for i in valid if isfinite(quality_factors[i]) && quality_factors[i] >= min_anchor_quality]; isempty(anchor_pool) && (anchor_pool = valid)
        anchor = anchor_pool[argmax(periods[anchor_pool])]
        anchor_branch_idx = _nearest_anchor_branch(@view(phase_velocity_branches[anchor, :]), phase_anchor_velocity, phase_velocity_range); anchor_branch_idx == 0 && return phase_velocities
        phase_velocities[anchor] = phase_velocity_branches[anchor, anchor_branch_idx]; selected_phase_branches[anchor] = phase_branch_numbers[anchor_branch_idx]
        sorted = valid[sortperm(periods[valid])]; anchor_pos = findfirst(==(anchor), sorted)
        for pos in (anchor_pos - 1):-1:1
            i = sorted[pos]; prev = sorted[pos + 1]; branch_idx = _nearest_finite_branch(@view(phase_velocity_branches[i, :]), phase_velocities[prev]); branch_idx == 0 && continue
            phase_velocities[i] = phase_velocity_branches[i, branch_idx]; selected_phase_branches[i] = phase_branch_numbers[branch_idx]
            phase_suspect[i] = abs(phase_velocities[i] - phase_velocities[prev]) / max(abs(phase_velocities[prev]), eps(Float64)) > phase_smoothness_jump
        end
        for pos in (anchor_pos + 1):length(sorted)
            i = sorted[pos]; prev = sorted[pos - 1]; branch_idx = _nearest_finite_branch(@view(phase_velocity_branches[i, :]), phase_velocities[prev]); branch_idx == 0 && continue
            phase_velocities[i] = phase_velocity_branches[i, branch_idx]; selected_phase_branches[i] = phase_branch_numbers[branch_idx]
            phase_suspect[i] = abs(phase_velocities[i] - phase_velocities[prev]) / max(abs(phase_velocities[prev]), eps(Float64)) > phase_smoothness_jump
        end
        u_predicted_from_phase .= compute_group_velocity_from_phase(periods, phase_velocities); return phase_velocities
    end
    function finish_phase_velocity_picks!(phase_velocities::Vector{Float64}, phase_velocity_branches::Matrix{Float64}, measured_phases::Vector{Float64}, selected_phase_branches::Vector{Int}, phase_suspect::AbstractVector{Bool}, u_predicted_from_phase::Vector{Float64}, periods::Vector{Float64}, frequencies::Vector{Float64}, arrival_times::Vector{Float64}, distance::Float64, quality_factors::Vector{Float64}, phase_branch_numbers::Vector{Int}, phvel_source_phase::Float64; min_anchor_quality::Float64=3.0, phase_anchor_velocity::Float64=3.3, phase_velocity_range::Tuple{Float64,Float64}=(2.5, 4.5), phase_smoothness_jump::Float64=0.05, phase_plausible_range::Tuple{Float64,Float64}=(0.5, 20.0))
        _fill_phase_velocity_branches!(phase_velocity_branches, frequencies, arrival_times, measured_phases, distance, phase_branch_numbers, phvel_source_phase; phase_plausible_range=phase_plausible_range)
        return resolve_phase_velocity_cycles!(phase_velocities, selected_phase_branches, phase_suspect, u_predicted_from_phase, phase_velocity_branches, periods, quality_factors, phase_branch_numbers; min_anchor_quality=min_anchor_quality, phase_anchor_velocity=phase_anchor_velocity, phase_velocity_range=phase_velocity_range, phase_smoothness_jump=phase_smoothness_jump)
    end
    
    function _interp_phase_velocity_prior(period::Float64, prior_periods, prior_velocities; fallback::Float64=3.3)
        (isnothing(prior_periods) || isnothing(prior_velocities)) && return fallback
        length(prior_periods) == length(prior_velocities) || return fallback
        valid = findall(i -> isfinite(Float64(prior_periods[i])) &&
                             isfinite(Float64(prior_velocities[i])) &&
                             Float64(prior_periods[i]) > 0.0 &&
                             Float64(prior_velocities[i]) > 0.0,
                        eachindex(prior_periods))
        isempty(valid) && return fallback
        ps = Float64[prior_periods[i] for i in valid]
        vs = Float64[prior_velocities[i] for i in valid]
        order = sortperm(ps)
        ps = ps[order]; vs = vs[order]
        period <= ps[1] && return vs[1]
        period >= ps[end] && return vs[end]
        hi = searchsortedfirst(ps, period)
        lo = hi - 1
        w = (period - ps[lo]) / (ps[hi] - ps[lo])
        return (1.0 - w) * vs[lo] + w * vs[hi]
    end

    function _accept_phtovel_velocity(c::Float64,
                                      phase_velocity_range::Tuple{Float64,Float64},
                                      phase_plausible_range::Tuple{Float64,Float64})
        isfinite(c) || return false
        phase_plausible_range[1] <= c <= phase_plausible_range[2] || return false
        phase_velocity_range[1] <= c <= phase_velocity_range[2] || return false
        return true
    end

    # Single propagation step: given an already-resolved period at idx_prev,
    # resolve idx_curr's cycle number/velocity from the group-velocity-predicted
    # continuation (pyaftan-style). Direction-agnostic — idx_prev may be at a
    # longer or shorter period than idx_curr. Mutates
    # phase_velocities/selected_phase_branches/phase_suspect[idx_curr] in place.
    function _phtovel_propagate_step!(phase_velocities::Vector{Float64},
                                      selected_phase_branches::Vector{Int},
                                      phase_suspect::AbstractVector{Bool},
                                      measured_phases::Vector{Float64},
                                      frequencies::Vector{Float64},
                                      arrival_times::Vector{Float64},
                                      dist::Float64,
                                      phvel_source_phase::Float64,
                                      idx_curr::Int,
                                      idx_prev::Int;
                                      phase_velocity_range::Tuple{Float64,Float64}=(2.5, 4.5),
                                      phase_smoothness_jump::Float64=0.05,
                                      phase_plausible_range::Tuple{Float64,Float64}=(0.5, 20.0))
        t_g_curr = arrival_times[idx_curr]
        t_g_prev = arrival_times[idx_prev]
        phi_curr = measured_phases[idx_curr]
        omega_curr = 2π * frequencies[idx_curr]
        omega_prev = 2π * frequencies[idx_prev]

        v_prev = phase_velocities[idx_prev]
        if !isfinite(v_prev)
            phase_suspect[idx_curr] = true
            return phase_velocities
        end

        # Group velocities from arrival times
        u_curr = dist / t_g_curr
        u_prev = dist / t_g_prev

        denom_pred = (1.0 / u_curr + 1.0 / u_prev) * (omega_curr - omega_prev) / 2.0 + omega_prev / v_prev
        if abs(denom_pred) <= 1e-6
            phase_suspect[idx_curr] = true
            return phase_velocities
        end
        vpred_curr = omega_curr / denom_pred

        # Unwrap phase using predicted velocity
        phpred_curr = omega_curr * (t_g_curr - dist / vpred_curr)
        k_curr = round(Int, (phpred_curr - phi_curr) / (2π))

        denom_curr = t_g_curr - (phi_curr + 2π * k_curr + phvel_source_phase) / omega_curr
        if abs(denom_curr) <= 1e-6
            phase_suspect[idx_curr] = true
            return phase_velocities
        end
        v_curr = dist / denom_curr
        selected_phase_branches[idx_curr] = k_curr
        jumpy = abs(v_curr - v_prev) / max(abs(v_prev), eps(Float64)) > phase_smoothness_jump
        if _accept_phtovel_velocity(v_curr, phase_velocity_range, phase_plausible_range)
            phase_velocities[idx_curr] = v_curr
            phase_suspect[idx_curr] = jumpy
        else
            phase_suspect[idx_curr] = true
        end
        return phase_velocities
    end

    # Shared inward-propagation loop: given an already-resolved starting point
    # (idx_long, k_long, v_long) at the longest valid period, walk toward
    # shorter periods predicting each step's phase velocity from the previous
    # step's group-velocity relationship (pyaftan-style). Mutates
    # phase_velocities/selected_phase_branches/phase_suspect in place for all
    # of `sorted` except idx_long itself (assumed already set by the caller).
    function _phtovel_propagate_inward!(phase_velocities::Vector{Float64},
                                        selected_phase_branches::Vector{Int},
                                        phase_suspect::AbstractVector{Bool},
                                        measured_phases::Vector{Float64},
                                        frequencies::Vector{Float64},
                                        arrival_times::Vector{Float64},
                                        dist::Float64,
                                        phvel_source_phase::Float64,
                                        sorted::Vector{Int},
                                        idx_long::Int;
                                        phase_velocity_range::Tuple{Float64,Float64}=(2.5, 4.5),
                                        phase_smoothness_jump::Float64=0.05,
                                        phase_plausible_range::Tuple{Float64,Float64}=(0.5, 20.0))
        long_pos = findfirst(==(idx_long), sorted)
        for pos in (long_pos - 1):-1:1
            _phtovel_propagate_step!(phase_velocities, selected_phase_branches,
                phase_suspect, measured_phases, frequencies, arrival_times, dist,
                phvel_source_phase, sorted[pos], sorted[pos + 1];
                phase_velocity_range=phase_velocity_range,
                phase_smoothness_jump=phase_smoothness_jump,
                phase_plausible_range=phase_plausible_range)
        end
        return phase_velocities
    end

    # pyaftan-style phase unwrapping using group-velocity-based prediction.
    function _phtovel_unwrap_phase_to_velocity!(phase_velocities::Vector{Float64},
                                                selected_phase_branches::Vector{Int},
                                                phase_suspect::AbstractVector{Bool},
                                                measured_phases::Vector{Float64},
                                                periods::Vector{Float64},
                                                frequencies::Vector{Float64},
                                                arrival_times::Vector{Float64},
                                                distance::Float64,
                                                phvel_source_phase::Float64;
                                                phase_anchor_velocity::Float64=3.3,
                                                phase_velocity_range::Tuple{Float64,Float64}=(2.5, 4.5),
                                                phase_smoothness_jump::Float64=0.05,
                                                phase_plausible_range::Tuple{Float64,Float64}=(0.5, 20.0),
                                                phtovel_prior_periods=nothing,
                                                phtovel_prior_velocities=nothing)
        fill!(phase_velocities, NaN)
        fill!(selected_phase_branches, 0)
        fill!(phase_suspect, false)
        dist = Float64(distance)

        # Find valid periods with both group time and measured phase
        valid = findall(i -> isfinite(periods[i]) && isfinite(arrival_times[i]) && isfinite(measured_phases[i]) && arrival_times[i] > 0.0 && periods[i] > 0.0, eachindex(periods))
        isempty(valid) && return phase_velocities

        # Sort by period (ascending)
        sorted = valid[sortperm(periods[valid])]

        # Start from longest period and unwrap downward.
        idx_long = sorted[end]
        t_g_long = arrival_times[idx_long]
        phi_long = measured_phases[idx_long]
        omega_long = 2π * frequencies[idx_long]
        vpred_long = _interp_phase_velocity_prior(periods[idx_long],
            phtovel_prior_periods, phtovel_prior_velocities;
            fallback=phase_anchor_velocity)
        phpred_long = omega_long * (t_g_long - dist / vpred_long)
        k_long = round(Int, (phpred_long - phi_long) / (2π))
        denom_long = t_g_long - (phi_long + 2π * k_long + phvel_source_phase) / omega_long
        if abs(denom_long) <= 1e-6
            phase_suspect[idx_long] = true
            return phase_velocities
        end
        v_long = dist / denom_long
        selected_phase_branches[idx_long] = k_long
        if _accept_phtovel_velocity(v_long, phase_velocity_range, phase_plausible_range)
            phase_velocities[idx_long] = v_long
        else
            phase_suspect[idx_long] = true
            return phase_velocities
        end

        # Propagate from longest to shortest period (descending index order)
        return _phtovel_propagate_inward!(phase_velocities, selected_phase_branches,
            phase_suspect, measured_phases, frequencies, arrival_times, dist,
            phvel_source_phase, sorted, idx_long;
            phase_velocity_range=phase_velocity_range,
            phase_smoothness_jump=phase_smoothness_jump,
            phase_plausible_range=phase_plausible_range)
    end
    
    function finish_phase_velocity_picks_phtovel!(phase_velocities::Vector{Float64},
                                                  measured_phases::Vector{Float64},
                                                  selected_phase_branches::Vector{Int},
                                                  phase_suspect::AbstractVector{Bool},
                                                  u_predicted_from_phase::Vector{Float64},
                                                  periods::Vector{Float64},
                                                  frequencies::Vector{Float64},
                                                  arrival_times::Vector{Float64},
                                                  distance::Float64,
                                                  phvel_source_phase::Float64;
                                                  phase_anchor_velocity::Float64=3.3,
                                                  phase_velocity_range::Tuple{Float64,Float64}=(2.5, 4.5),
                                                  phase_smoothness_jump::Float64=0.05,
                                                  phase_plausible_range::Tuple{Float64,Float64}=(0.5, 20.0),
                                                  phtovel_prior_periods=nothing,
                                                  phtovel_prior_velocities=nothing)
        _phtovel_unwrap_phase_to_velocity!(phase_velocities, selected_phase_branches,
                                           phase_suspect, measured_phases, periods,
                                           frequencies, arrival_times, distance,
                                           phvel_source_phase;
                                           phase_anchor_velocity=phase_anchor_velocity,
                                           phase_velocity_range=phase_velocity_range,
                                           phase_smoothness_jump=phase_smoothness_jump,
                                           phase_plausible_range=phase_plausible_range,
                                           phtovel_prior_periods=phtovel_prior_periods,
                                           phtovel_prior_velocities=phtovel_prior_velocities)
        u_predicted_from_phase .= compute_group_velocity_from_phase(periods, phase_velocities)

        return phase_velocities
    end

end

function _parse_pair_label(label::AbstractString)
    parts = split(label, "-")
    length(parts) == 2 || throw(ArgumentError("Cannot parse pair label: $(label)"))
    return (String(parts[1]), String(parts[2]))
end

"""
    pdsurftomo_dispersion_rows(per_pair_velocities::Dict{String,Vector{Float64}},
        mft_results_by_pair::Dict{String,MFTResult},
        station_coords::Dict{String,Tuple{Float64,Float64}};
        wavelength_ref_velocity=nothing, wavelength_fraction=nothing)
        -> Vector{NamedTuple}

Build flat dispersion rows `(period, lat1, lon1, lat2, lon2, velocity)` for
`write_pdsurftomo_dispersion`/pDSurfTomo, from any per-pair velocity vectors
aligned to each pair's `MFTResult.periods`. `lat1,lon1`/`lat2,lon2` follow the
pair label's station order (`"STA1-STA2"`). NaN velocities are skipped. If both
`wavelength_ref_velocity` and `wavelength_fraction` are given, periods failing
`wavelength_valid_period` (using that pair's `MFTResult.distance`) are skipped.
"""
function pdsurftomo_dispersion_rows(per_pair_velocities::Dict{String,Vector{Float64}},
        mft_results_by_pair::Dict{String,MFTResult},
        station_coords::Dict{String,Tuple{Float64,Float64}};
        wavelength_ref_velocity::Union{Nothing,Real}=nothing,
        wavelength_fraction::Union{Nothing,Real}=nothing)
    rows = NamedTuple[]
    for (pair_label, velocities) in per_pair_velocities
        sta, stb = _parse_pair_label(pair_label)
        (haskey(station_coords, sta) && haskey(station_coords, stb)) || continue
        lat1, lon1 = station_coords[sta]
        lat2, lon2 = station_coords[stb]
        res = mft_results_by_pair[pair_label]
        for ip in eachindex(res.periods)
            period = Float64(res.periods[ip])
            v = velocities[ip]
            isfinite(v) && v > 0 || continue
            wavelength_valid_period(period, res.distance;
                wavelength_ref_velocity=wavelength_ref_velocity,
                wavelength_fraction=wavelength_fraction) || continue
            push!(rows, (; period, lat1, lon1, lat2, lon2, velocity=v))
        end
    end
    sort!(rows, by=r -> (r.period, r.lat1, r.lon1))
    return rows
end

"""
    write_pdsurftomo_dispersion(path::AbstractString, rows) -> String

Write flat dispersion rows (`period, lat1, lon1, lat2, lon2, velocity`) to a
whitespace-separated text file with no header — the input format
`pDSurfTomo_v1.jl` reads directly (Section 1, "Flat dispersion file").
"""
function write_pdsurftomo_dispersion(path::AbstractString, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        for row in rows
            @printf(io, "%.8g %.8f %.8f %.8f %.8f %.8f\n",
                Float64(row.period), Float64(row.lat1), Float64(row.lon1),
                Float64(row.lat2), Float64(row.lon2), Float64(row.velocity))
        end
    end
    return path
end
