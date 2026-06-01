function zero_lag_correlation(filtered_traces_causal::AbstractMatrix{<:Real},
                              filtered_traces_acausal::AbstractMatrix{<:Real})
    nfreq = size(filtered_traces_causal, 2)
    correlations = Float64[]
    
    for i in 1:nfreq
        fc = filtered_traces_causal[:, i]
        fa = filtered_traces_acausal[:, i]
        
        # Normalize
        fc_norm = fc .- mean(fc)
        fa_norm = fa .- mean(fa)
        
        # Compute correlation
        numerator = sum(fc_norm .* fa_norm)
        denominator = sqrt(sum(fc_norm.^2) * sum(fa_norm.^2))
        
        if denominator > 0
            corr = numerator / denominator
            push!(correlations, clamp(corr, -1.0, 1.0))  # Clamp to [-1, 1] range
        else
            push!(correlations, NaN)
        end
    end
    
    return correlations
end

begin
	"""
	    analyze_causal_acausal_branches(w_c, w_ac, dt, dist, periods; ...)
	        -> BranchAnalysisResult
	
	Array-first single-state branch analysis. A temporary `MFTFilterBank` is built
	from the fixed `dt`, raw waveform length, and requested periods, then both
	branches are processed together through the banked batch path.
	"""
	function analyze_causal_acausal_branches(w_c::AbstractVector{<:Real},
	                                         w_ac::AbstractVector{<:Real},
	                                         dt::Real,
	                                         dist::Real,
	                                         periods::Vector{Float64};
	                                         velocity_range::Tuple{Float64,Float64}=(2.0, 6.0),
	                                         bandwidth_factor::Float64=sqrt(2.0 / 25.0),
	                                         zero_pad_factor::Int=4,
	                                         upsample_factor::Real=2.0,
	                                         precision::Type{<:AbstractFloat}=Float32,
	                                         storage_mode::Symbol=:picks_only,
	                                         kwargs...)
	    bank = MFTFilterBank(Float64(dt), length(w_c), periods;
	                         bandwidth_factor=bandwidth_factor,
	                         zero_pad_factor=zero_pad_factor,
	                         upsample_factor=upsample_factor,
	                         velocity_range=velocity_range,
	                         precision=precision,
	                         storage_mode=storage_mode,
	                         N_initial=2)
	    return analyze_causal_acausal_branches(w_c, w_ac, Float64(dist), bank; kwargs...)
	end
	
	"""
	    analyze_causal_acausal_branches(W_c, W_ac, dt, dists, periods; ...)
	        -> BranchBatchAnalysisResult
	
	Array-first N-D branch analysis. `W_c` and `W_ac` have shape
	`(nt × trailing dims...)`; `dists` is scalar or has shape `trailing dims`.
	All causal and acausal waveforms are processed in one banked batch.
	"""
	function analyze_causal_acausal_branches(W_c::AbstractArray{<:Real},
	                                         W_ac::AbstractArray{<:Real},
	                                         dt::Real,
	                                         dists::Union{Real,AbstractArray{<:Real}},
	                                         periods::Vector{Float64};
	                                         velocity_range::Tuple{Float64,Float64}=(2.0, 6.0),
	                                         bandwidth_factor::Float64=sqrt(2.0 / 25.0),
	                                         zero_pad_factor::Int=4,
	                                         upsample_factor::Real=2.0,
	                                         precision::Type{<:AbstractFloat}=Float32,
	                                         storage_mode::Symbol=:picks_only,
	                                         kwargs...)
	    ndims(W_c) >= 2 || throw(ArgumentError("W_c must have shape (nt × trailing dims...)"))
	    bank = MFTFilterBank(Float64(dt), size(W_c, 1), periods;
	                         bandwidth_factor=bandwidth_factor,
	                         zero_pad_factor=zero_pad_factor,
	                         upsample_factor=upsample_factor,
	                         velocity_range=velocity_range,
	                         precision=precision,
	                         storage_mode=storage_mode,
	                         N_initial=2 * prod(size(W_c)[2:end]))
	    return analyze_causal_acausal_branches(W_c, W_ac, dists, bank; kwargs...)
	end
	

end

function high_correlation_indices(result::BranchAnalysisResult; correlation_threshold::Float64=0.85)
    return findall(c -> isfinite(c) && c >= correlation_threshold, result.branch_correlation)
end

# ╔═╡ 4f47e6c2-0f17-4f71-a7da-08f651176632
"""
    high_correlation_indices(result::BranchBatchAnalysisResult; correlation_threshold=0.85)

Return a boolean mask [period, state] selecting finite correlation values above threshold.
"""
function high_correlation_indices(result::BranchBatchAnalysisResult; correlation_threshold::Float64=0.85)
    return map(c -> isfinite(c) && c >= correlation_threshold, result.branch_correlation)
end

function _matched_peak_average_velocities(causal_peaks::Vector{Tuple{Float64,Float64}},
                                          acausal_peaks::Vector{Tuple{Float64,Float64}},
                                          distance::Float64;
                                          velocity_tolerance_fraction::Float64=0.10)
    velocity_tolerance_fraction >= 0.0 || throw(ArgumentError("velocity_tolerance_fraction must be >= 0"))

    causal_vels = Float64[]
    acausal_vels = Float64[]

    for (t, _) in causal_peaks
        if isfinite(t) && t > 0.0
            v = distance / t
            isfinite(v) && v > 0.0 && push!(causal_vels, v)
        end
    end
    for (t, _) in acausal_peaks
        if isfinite(t) && t > 0.0
            v = distance / t
            isfinite(v) && v > 0.0 && push!(acausal_vels, v)
        end
    end

    isempty(causal_vels) && return Float64[]
    isempty(acausal_vels) && return Float64[]

    used_acausal = falses(length(acausal_vels))
    matched_avg = Float64[]

    for vc in causal_vels
        best_j = 0
        best_rel = Inf
        for (j, va) in enumerate(acausal_vels)
            used_acausal[j] && continue
            denom = max((abs(vc) + abs(va)) / 2.0, eps(Float64))
            rel = abs(vc - va) / denom
            if rel <= velocity_tolerance_fraction && rel < best_rel
                best_rel = rel
                best_j = j
            end
        end

        if best_j > 0
            used_acausal[best_j] = true
            push!(matched_avg, 0.5 * (vc + acausal_vels[best_j]))
        end
    end

    sort!(matched_avg)
    return matched_avg
end

# ╔═╡ fb000010-0000-0000-0000-000000000001
"""
    analyze_causal_acausal_branches(W_c, W_ac, dists, bank; ...)
        -> BranchBatchAnalysisResult

Array-based, bank-accelerated multi-state analysis. `W_c` and `W_ac` may be
2-D `(nt × nstates)` matrices or N-D arrays `(nt × d1 × d2 × ...)`; each
trailing index is one waveform/state. `dists` is either a scalar distance or an
array with shape `(d1 × d2 × ...)`. All 2·nstates waveforms are FFT'd in a
single batch (one `_run_phase1!` call), then peak-finding and assembly use the
per-state distances.

# Arguments
- `W_c`   : causal waveforms, shape `(nt × trailing dims...)`
- `W_ac`  : acausal waveforms, same shape as `W_c`
- `dists` : distances [km], either scalar or shape `trailing dims`
- `bank`  : pre-computed `MFTFilterBank` (sets periods, bandwidth, etc.)
"""
function analyze_causal_acausal_branches(W_c::AbstractArray{<:Real},
                                         W_ac::AbstractArray{<:Real},
                                         dists::Union{Real,AbstractArray{<:Real}},
                                         bank::MFTFilterBank;
                                         state_labels=nothing,
                                         max_modes::Int=4,
                                         compute_correlation::Bool=true,
                                         compute_phase::Bool=false,
                                         phvel_correction::Float64=0.0,
                                         phvel_source_phase::Float64=phvel_correction,
                                         kwargs...)
    ndims(W_c) >= 2 || throw(ArgumentError("W_c must have shape (nt × trailing dims...)"))
    size(W_ac) == size(W_c) || throw(ArgumentError("W_c and W_ac must have the same shape"))
    nt = size(W_c, 1)
    trailing_dims = size(W_c)[2:end]
    nstates = prod(trailing_dims)
    nstates > 0 || throw(ArgumentError("At least one waveform/state is required"))

    dists_vec = if dists isa Real
        fill(Float64(dists), nstates)
    else
        size(dists) == trailing_dims ||
            throw(ArgumentError("dists shape $(size(dists)) must match waveform trailing dims $(trailing_dims)"))
        vec(Float64.(dists))
    end

    labels = if isnothing(state_labels)
        ["State $(i)" for i in 1:nstates]
    else
        length(state_labels) == nstates ||
            throw(ArgumentError("state_labels length must match number of waveform states"))
        vec(string.(state_labels))
    end

    periods = bank.periods
    T = _bank_float_type(bank)
    W_c_flat = reshape(T.(W_c), nt, nstates)
    W_ac_flat = reshape(T.(W_ac), nt, nstates)
    W_all = Matrix{T}(undef, nt, 2 * nstates)
    W_all[:, 1:nstates] .= W_c_flat
    W_all[:, nstates+1:2*nstates] .= W_ac_flat
    dists_all = vcat(dists_vec, dists_vec)

    _require_batch_size(bank, 2 * nstates)
    correlations_picks = nothing
    results_all = if bank.storage_mode == :picks_only
        correlation_pairs = compute_correlation ?
            [(i, nstates + i) for i in 1:nstates] : nothing
        results, corr = _run_picks_only!(bank, W_all, dists_all, 2 * nstates;
                                         compute_phase=compute_phase,
                                         phvel_source_phase=phvel_source_phase,
                                         correlation_pairs=correlation_pairs,
                                         kwargs...)
        correlations_picks = corr
        results
    else
        _run_phase1!(bank, W_all, 2 * nstates)
        _run_phase2!(bank, dists_all, 2 * nstates;
                     compute_phase=compute_phase,
                     phvel_source_phase=phvel_source_phase,
                     kwargs...)
    end

    state_results = BranchAnalysisResult[]
    corr_matrix   = fill(NaN, length(periods), nstates)
    for i in 1:nstates
        res_c = results_all[i]
        res_a = results_all[nstates + i]
        d     = dists_vec[i]
        modes_c = extract_all_modes(res_c; distance=d, max_modes=max_modes)
        modes_a = extract_all_modes(res_a; distance=d, max_modes=max_modes)
        correlations = fill(NaN, length(periods))
        if compute_correlation
            correlations = bank.storage_mode == :picks_only ?
                correlations_picks[:, i] :
                zero_lag_correlation(res_c.filtered_traces, res_a.filtered_traces)
        end
        push!(state_results,
              BranchAnalysisResult(res_c, res_a, modes_c, modes_a, correlations, periods, d))
        corr_matrix[:, i] = correlations
    end
    return BranchBatchAnalysisResult(state_results, corr_matrix, periods, labels)
end

# ╔═╡ fb000011-0000-0000-0000-000000000001
"""
    analyze_causal_acausal_branches(w_c, w_ac, dist, bank; ...)
        -> BranchAnalysisResult

Array-based, bank-accelerated single-state analysis. Both vectors are FFT'd
together in a 2-column batch.

# Arguments
- `w_c`  : length-`nt` causal waveform vector
- `w_ac` : length-`nt` acausal waveform vector
- `dist` : source–receiver distance [km]
- `bank` : pre-computed `MFTFilterBank`
"""
function analyze_causal_acausal_branches(w_c::AbstractVector{<:Real},
                                         w_ac::AbstractVector{<:Real},
                                         dist::Real,
                                         bank::MFTFilterBank;
                                         max_modes::Int=4,
                                         compute_correlation::Bool=true,
                                         compute_phase::Bool=false,
                                         phvel_correction::Float64=0.0,
                                         phvel_source_phase::Float64=phvel_correction,
                                         kwargs...)
    W_c  = reshape(w_c,  :, 1)
    W_ac = reshape(w_ac, :, 1)
    res = analyze_causal_acausal_branches(W_c, W_ac, dist, bank;
              max_modes=max_modes, compute_correlation=compute_correlation,
              compute_phase=compute_phase, phvel_source_phase=phvel_source_phase,
              state_labels=["state 1"], kwargs...)
    return only(res.state_results)
end

# ╔═╡ fb000012-0000-0000-0000-000000000001
"""
    analyze_self_branches(W, dists, bank; ...)
        -> BranchBatchAnalysisResult

Analyze waveforms once and wrap each result as both causal and acausal branch.
This is useful for synthetic/codebook mixture candidates where the two branches
are intentionally identical and `analyze_causal_acausal_branches(W, W, ...)`
would do twice the FFT/IFFT work for no new information.
"""
function analyze_self_branches(W::AbstractArray{<:Real},
                               dists::Union{Real,AbstractArray{<:Real}},
                               bank::MFTFilterBank;
                               state_labels=nothing,
                               max_modes::Int=4,
                               compute_phase::Bool=false,
                               phvel_correction::Float64=0.0,
                               phvel_source_phase::Float64=phvel_correction,
                               kwargs...)
    ndims(W) >= 2 || throw(ArgumentError("W must have shape (nt × trailing dims...)"))
    nt = size(W, 1)
    trailing_dims = size(W)[2:end]
    nstates = prod(trailing_dims)
    nstates > 0 || throw(ArgumentError("At least one waveform/state is required"))

    dists_vec = if dists isa Real
        fill(Float64(dists), nstates)
    else
        size(dists) == trailing_dims ||
            throw(ArgumentError("dists shape $(size(dists)) must match waveform trailing dims $(trailing_dims)"))
        vec(Float64.(dists))
    end

    labels = if isnothing(state_labels)
        ["State $(i)" for i in 1:nstates]
    else
        length(state_labels) == nstates ||
            throw(ArgumentError("state_labels length must match number of waveform states"))
        vec(string.(state_labels))
    end

    results = perform_mft_analysis_batch!(bank, reshape(W, nt, nstates), dists_vec;
                                          compute_phase=compute_phase,
                                          phvel_source_phase=phvel_source_phase,
                                          kwargs...)

    periods = bank.periods
    state_results = BranchAnalysisResult[]
    corr_matrix = ones(Float64, length(periods), nstates)
    for i in 1:nstates
        res = results[i]
        d = dists_vec[i]
        modes = extract_all_modes(res; distance=d, max_modes=max_modes)
        push!(state_results,
              BranchAnalysisResult(res, res, modes, modes, corr_matrix[:, i], periods, d))
    end
    return BranchBatchAnalysisResult(state_results, corr_matrix, periods, labels)
end

# ╔═╡ f4e12fc4-b775-4c17-9459-920ad3c8bbd6
begin
	"""
	    wavelength_valid_period(period, distance; wavelength_ref_velocity=nothing,
	                            wavelength_fraction=nothing)
	
	Return whether `period` passes the optional distance-dependent wavelength filter.
	When filtering is requested, keep periods satisfying
	`wavelength_ref_velocity * period < wavelength_fraction * distance`.
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
	    return ref_velocity * period < fraction * distance
	end
	
	function _wavelength_valid_indices(periods::AbstractVector{<:Real}, distance::Real;
	                                   wavelength_ref_velocity::Union{Nothing,Real}=nothing,
	                                   wavelength_fraction::Union{Nothing,Real}=nothing)
	    return findall(period -> wavelength_valid_period(period, distance;
	                                                     wavelength_ref_velocity=wavelength_ref_velocity,
	                                                     wavelength_fraction=wavelength_fraction),
	                   periods)
	end
end
