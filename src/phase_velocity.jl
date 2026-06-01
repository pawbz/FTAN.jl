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
        for pos in (length(sorted) - 1):-1:1
            idx_curr = sorted[pos]
            idx_prev = sorted[pos + 1]
            
            t_g_curr = arrival_times[idx_curr]
            t_g_prev = arrival_times[idx_prev]
            phi_curr = measured_phases[idx_curr]
            omega_curr = 2π * frequencies[idx_curr]
            omega_prev = 2π * frequencies[idx_prev]
            
            v_prev = phase_velocities[idx_prev]
            isfinite(v_prev) || (phase_suspect[idx_curr] = true; continue)
            
            # Group velocities from arrival times
            u_curr = dist / t_g_curr
            u_prev = dist / t_g_prev
            
            denom_pred = (1.0 / u_curr + 1.0 / u_prev) * (omega_curr - omega_prev) / 2.0 + omega_prev / v_prev
            if abs(denom_pred) <= 1e-6
                phase_suspect[idx_curr] = true
                continue
            end
            vpred_curr = omega_curr / denom_pred
            
            # Unwrap phase using predicted velocity
            phpred_curr = omega_curr * (t_g_curr - dist / vpred_curr)
            k_curr = round(Int, (phpred_curr - phi_curr) / (2π))
            
            denom_curr = t_g_curr - (phi_curr + 2π * k_curr + phvel_source_phase) / omega_curr
            if abs(denom_curr) <= 1e-6
                phase_suspect[idx_curr] = true
                continue
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
        end
        
        return phase_velocities
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

begin
    _wrap_pi(x::Real) = mod(Float64(x) + π, 2π) - π

    function _circular_mean_pair(phi1::Real, phi2::Real, w1::Real, w2::Real)
        z = Float64(w1) * cis(Float64(phi1)) + Float64(w2) * cis(Float64(phi2))
        abs(z) <= eps(Float64) && return NaN
        return angle(z)
    end

    function _sigma_phi_from_quality(quality::Real; floor::Float64=0.05, cap::Float64=π)
        q = Float64(quality)
        isfinite(q) && q > 0 || return cap
        return clamp(1.0 / q, floor, cap)
    end

    function folded_path_azimuth_deg(lat1::Real, lon1::Real, lat2::Real, lon2::Real)
        vals = Float64.((lat1, lon1, lat2, lon2))
        all(isfinite, vals) || return NaN
        φ1, φ2 = deg2rad(vals[1]), deg2rad(vals[3])
        Δλ = deg2rad(vals[4] - vals[2])
        y = sin(Δλ) * cos(φ2)
        x = cos(φ1) * sin(φ2) - sin(φ1) * cos(φ2) * cos(Δλ)
        az = mod(rad2deg(atan(y, x)), 360.0)
        return az >= 180.0 ? az - 180.0 : az
    end
end

# ╔═╡ a2000005-0000-0000-0000-000000000001
function phase_regression_observations(branch_results_by_pair;
        branch_policy::Symbol=:symmetric_average,
        sigma_phi=nothing,
        sigma_from_quality::Bool=true,
        sigma_phi_floor::Float64=0.05,
        sigma_phi_cap::Float64=π,
        phvel_source_phase::Float64=π / 4,
        pair_azimuths=Dict{String,Float64}(),
        min_quality::Float64=0.0)
    branch_policy == :symmetric_average ||
        throw(ArgumentError("Only branch_policy=:symmetric_average is implemented"))

    rows = PhaseRegressionObservation[]
    for pair_key in sort(collect(keys(branch_results_by_pair)); by=string)
        pair_label = String(pair_key)
        result = branch_results_by_pair[pair_key]
        result isa BranchAnalysisResult || continue
        az = get(pair_azimuths, pair_label, NaN)
        for ip in eachindex(result.periods)
            c = result.causal_result
            a = result.acausal_result
            period = Float64(result.periods[ip])
            freq = 1.0 / period
            omega = 2π * freq
            t_c, t_a = c.arrival_times[ip], a.arrival_times[ip]
            phi_c, phi_a = c.measured_phases[ip], a.measured_phases[ip]
            q_c, q_a = c.quality_factors[ip], a.quality_factors[ip]
            all(isfinite, (period, t_c, t_a, phi_c, phi_a, q_c, q_a)) || continue
            (period > 0.0 && t_c > 0.0 && t_a > 0.0 && q_c >= min_quality && q_a >= min_quality) || continue

            wc = max(Float64(q_c), eps(Float64))
            wa = max(Float64(q_a), eps(Float64))
            t_peak = (wc * Float64(t_c) + wa * Float64(t_a)) / (wc + wa)
            phi_meas = _circular_mean_pair(phi_c, phi_a, wc, wa)
            isfinite(phi_meas) || continue
            quality = 0.5 * (Float64(q_c) + Float64(q_a))
            sigma = if !isnothing(sigma_phi)
                sigma_phi isa Function ? Float64(sigma_phi(pair_label, period)) : Float64(sigma_phi)
            elseif sigma_from_quality
                _sigma_phi_from_quality(quality; floor=sigma_phi_floor, cap=sigma_phi_cap)
            else
                sigma_phi_cap
            end
            isfinite(sigma) && sigma > 0.0 || continue
            base_phase = omega * t_peak - phi_meas - phvel_source_phase
            push!(rows, PhaseRegressionObservation(
                pair_label, period, freq, omega, Float64(result.distance),
                t_peak, phi_meas, base_phase, sigma, quality, Float64(az)))
        end
    end
    return rows
end

# ╔═╡ a2000006-0000-0000-0000-000000000001
begin
    function _phase_cycles_for_slope(obs::Vector{PhaseRegressionObservation}, slope::Float64)
        [round(Int, (slope * o.distance - o.base_phase) / (2π)) for o in obs]
    end

    _phase_unwrapped(obs::Vector{PhaseRegressionObservation}, cycles::Vector{Int}) =
        [obs[i].base_phase + 2π * cycles[i] for i in eachindex(obs)]

    _phase_weights(obs::Vector{PhaseRegressionObservation}) =
        [1.0 / max(o.sigma_phi^2, eps(Float64)) for o in obs]

    function _weighted_origin_slope(distances, phases, weights)
        den = sum(weights[i] * distances[i]^2 for i in eachindex(distances))
        den > 0 || return NaN
        return sum(weights[i] * distances[i] * phases[i] for i in eachindex(distances)) / den
    end

    function _robust_scale(xs)
        vals = [abs(Float64(x)) for x in xs if isfinite(x)]
        isempty(vals) && return NaN
        med = median(vals)
        mad = median(abs.(vals .- med))
        return max(1.4826 * mad, median(vals), eps(Float64))
    end

    function _candidate_slopes(obs::Vector{PhaseRegressionObservation}, period::Float64;
            velocity_range::Tuple{Float64,Float64},
            cycle_range::UnitRange{Int},
            prior_c::Union{Nothing,Real}=nothing,
            prior_width_fraction::Float64=0.20,
            max_candidates::Int=400)
        vmin, vmax = velocity_range
        omega = 2π / period
        slopes = Float64[omega / vmax, omega / vmin]
        if !isnothing(prior_c) && isfinite(Float64(prior_c)) && Float64(prior_c) > 0
            cp = Float64(prior_c)
            append!(slopes, omega ./ [cp, cp * (1 - prior_width_fraction), cp * (1 + prior_width_fraction)])
        end
        for o in obs
            o.distance > 0 || continue
            for n in cycle_range
                phi = o.base_phase + 2π * n
                phi > 0 || continue
                slope = phi / o.distance
                c = omega / slope
                vmin <= c <= vmax && push!(slopes, slope)
            end
        end
        slopes = sort(unique(round.(filter(isfinite, slopes); digits=10)))
        if length(slopes) <= max_candidates
            return slopes
        end
        step = max(1, floor(Int, length(slopes) / max_candidates))
        return slopes[1:step:end][1:min(max_candidates, length(slopes[1:step:end]))]
    end
end

# ╔═╡ a2000007-0000-0000-0000-000000000001
function fit_phase_velocity_period(observations::AbstractVector{PhaseRegressionObservation},
        period::Real;
        prior_c::Union{Nothing,Real}=nothing,
        velocity_range::Tuple{Float64,Float64}=(2.5, 4.5),
        cycle_range::UnitRange{Int}=-80:80,
        ransac_tolerance::Float64=2.5,
        min_inliers::Int=4,
        max_refit_iterations::Int=4,
        phvel_source_phase::Float64=π / 4,
        azimuth_rejection::Bool=true,
        azimuth_bin_width::Float64=20.0,
        azimuth_threshold::Float64=2.5,
        prior_penalty::Float64=0.05)
    obs = [o for o in observations if isfinite(o.period) && isapprox(o.period, Float64(period); rtol=1e-8, atol=1e-8)]
    isempty(obs) && return PhaseVelocityPeriodFit(Float64(period), 1.0 / Float64(period), 2π / Float64(period),
        NaN, NaN, NaN, NaN, PhaseRegressionObservation[], Int[], Float64[], Float64[], Bool[], Bool[], Float64[],
        isnothing(prior_c) ? NaN : Float64(prior_c), -Inf)

    omega = 2π / Float64(period)
    weights = _phase_weights(obs)
    distances = [o.distance for o in obs]
    candidates = _candidate_slopes(obs, Float64(period);
        velocity_range=velocity_range, cycle_range=cycle_range, prior_c=prior_c)

    best = nothing
    for slope0 in candidates
        slope = slope0
        cycles = _phase_cycles_for_slope(obs, slope)
        phases = _phase_unwrapped(obs, cycles)
        residuals = phases .- slope .* distances
        inlier = [abs(residuals[i]) / obs[i].sigma_phi <= ransac_tolerance for i in eachindex(obs)]
        for _ in 1:max_refit_iterations
            count(inlier) >= 2 || break
            inds = findall(inlier)
            slope_new = _weighted_origin_slope(distances[inds], phases[inds], weights[inds])
            isfinite(slope_new) || break
            slope = slope_new
            cycles = _phase_cycles_for_slope(obs, slope)
            phases = _phase_unwrapped(obs, cycles)
            residuals = phases .- slope .* distances
            inlier = [abs(residuals[i]) / obs[i].sigma_phi <= ransac_tolerance for i in eachindex(obs)]
        end
        support = sum((weights[i] for i in eachindex(obs) if inlier[i]); init=0.0)
        rss = sum((weights[i] * residuals[i]^2 for i in eachindex(obs) if inlier[i]); init=0.0)
        prior_cost = if isnothing(prior_c) || !isfinite(Float64(prior_c))
            0.0
        else
            c0 = omega / slope
            prior_penalty * abs(c0 - Float64(prior_c)) / max(abs(Float64(prior_c)), eps(Float64))
        end
        score = support - rss - prior_cost
        if isnothing(best) || score > best.score
            best = (; score, slope, cycles, phases, residuals, inlier)
        end
    end

    if isnothing(best)
        best = (; score=-Inf, slope=NaN, cycles=zeros(Int, length(obs)),
            phases=fill(NaN, length(obs)), residuals=fill(NaN, length(obs)),
            inlier=falses(length(obs)))
    end

    azimuth_outlier = falses(length(obs))
    inlier = copy(best.inlier)
    slope = best.slope
    cycles = copy(best.cycles)
    phases = copy(best.phases)
    residuals = copy(best.residuals)

    if azimuth_rejection && count(inlier) >= min_inliers
        finite_az = [i for i in eachindex(obs) if inlier[i] && isfinite(obs[i].azimuth_deg)]
        if !isempty(finite_az)
            scale = _robust_scale(residuals[finite_az] ./ [obs[i].sigma_phi for i in finite_az])
            if isfinite(scale) && scale > 0
                bins = Dict{Int,Vector{Int}}()
                for i in finite_az
                    b = floor(Int, obs[i].azimuth_deg / azimuth_bin_width)
                    push!(get!(bins, b, Int[]), i)
                end
                for inds in values(bins)
                    length(inds) < 2 && continue
                    zmed = median(residuals[inds] ./ [obs[i].sigma_phi for i in inds])
                    if abs(zmed) > azimuth_threshold * scale
                        azimuth_outlier[inds] .= true
                    end
                end
                inlier .&= .!azimuth_outlier
                if count(inlier) >= 2
                    inds = findall(inlier)
                    slope = _weighted_origin_slope(distances[inds], phases[inds], weights[inds])
                    cycles = _phase_cycles_for_slope(obs, slope)
                    phases = _phase_unwrapped(obs, cycles)
                    residuals = phases .- slope .* distances
                    inlier = [inlier[i] && abs(residuals[i]) / obs[i].sigma_phi <= ransac_tolerance for i in eachindex(obs)]
                end
            end
        end
    end

    n_in = count(inlier)
    cfit = isfinite(slope) && slope > 0 ? omega / slope : NaN
    if n_in >= max(min_inliers, 2) && isfinite(slope)
        inds = findall(inlier)
        den = sum(weights[i] * distances[i]^2 for i in inds)
        dof = max(n_in - 1, 1)
        wrss = sum(weights[i] * residuals[i]^2 for i in inds)
        sigma_slope = den > 0 ? sqrt((wrss / dof) / den) : NaN
        sigma_c = isfinite(cfit) && isfinite(sigma_slope) ? abs(omega / slope^2) * sigma_slope : NaN
    else
        sigma_slope = NaN
        sigma_c = NaN
        cfit = NaN
    end

    return PhaseVelocityPeriodFit(Float64(period), 1.0 / Float64(period), omega,
        cfit, sigma_c, slope, sigma_slope, obs, cycles, phases, residuals,
        collect(inlier), collect(azimuth_outlier), weights,
        isnothing(prior_c) ? NaN : Float64(prior_c), best.score)
end

# ╔═╡ a2000008-0000-0000-0000-000000000001
function fit_phase_velocity_curve(observations::AbstractVector{PhaseRegressionObservation};
        period_order::Symbol=:descending,
        warm_start::Bool=true,
        kwargs...)
    periods = sort(unique(o.period for o in observations), rev=(period_order == :descending))
    fits = PhaseVelocityPeriodFit[]
    prior = nothing
    for period in periods
        fit = fit_phase_velocity_period(observations, period; prior_c=warm_start ? prior : nothing, kwargs...)
        push!(fits, fit)
        if isfinite(fit.phase_velocity) && fit.phase_velocity > 0
            prior = fit.phase_velocity
        end
    end
    return PhaseVelocityCurveFit(
        Float64[f.period for f in fits],
        Float64[f.phase_velocity for f in fits],
        Float64[f.sigma_c for f in fits],
        fits,
    )
end

