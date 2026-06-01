function perform_mft_analysis(data_in::AbstractVector{<:Real},
                             dt_original::Real,
                             distance::Real,
                             periods::Vector{Float64};
                             velocity_range::Tuple{Float64,Float64}=(2.0, 6.0),
                             bandwidth_factor::Float64=sqrt(2.0 / 25.0),
                             zero_pad_factor::Int=4,
                             upsample_factor::Real=2.0,
                             phvel_correction::Float64=0.0,
                             phvel_source_phase::Float64=phvel_correction,
                             compute_phase::Bool=true,
                             use_phtovel::Bool=false,
                             min_anchor_quality::Float64=3.0,
                             phase_anchor_velocity::Float64=3.3,
                             phase_velocity_range::Tuple{Float64,Float64}=(2.5, 4.5),
                             phase_smoothness_jump::Float64=0.05,
                             phase_plausible_range::Tuple{Float64,Float64}=(0.5, 20.0),
                             phtovel_prior_periods=nothing,
                             phtovel_prior_velocities=nothing,
                             precision::Type{<:AbstractFloat}=Float32,
                             storage_mode::Symbol=:picks_only)

    if @isdefined(MFTFilterBank)
        bank = MFTFilterBank(dt_original, length(data_in), periods;
                             velocity_range=velocity_range,
                             bandwidth_factor=bandwidth_factor,
                             zero_pad_factor=zero_pad_factor,
                             upsample_factor=upsample_factor,
                             precision=precision,
                             storage_mode=storage_mode,
                             N_initial=1)
        return only(perform_mft_analysis_batch!(
            bank, reshape(data_in, :, 1), distance;
            compute_phase=compute_phase,
            use_phtovel=use_phtovel,
            phvel_source_phase=phvel_source_phase,
            min_anchor_quality=min_anchor_quality,
            phase_anchor_velocity=phase_anchor_velocity,
            phase_velocity_range=phase_velocity_range,
            phase_smoothness_jump=phase_smoothness_jump,
            phase_plausible_range=phase_plausible_range,
            phtovel_prior_periods=phtovel_prior_periods,
            phtovel_prior_velocities=phtovel_prior_velocities,
        ))
    end

    upsample = Float64(upsample_factor)
    upsample > 0.0 || throw(ArgumentError("upsample_factor must be positive"))

    # Modest oversampling improves envelope peak localization without the 10× cost.
    data_unpadded = upsample == 1.0 ? Float64.(data_in) : DSP.resample(Float64.(data_in), upsample)
    dt            = Float64(dt_original) / upsample
    dist          = Float64(distance)

    if !(zero_pad_factor in (1, 2, 4, 8))
        throw(ArgumentError("zero_pad_factor must be one of 1, 2, 4, 8"))
    end

    npts_original = length(data_unpadded)
    npts_padded   = npts_original * zero_pad_factor
    data = if zero_pad_factor == 1
        data_unpadded
    else
        vcat(data_unpadded, zeros(Float64, npts_padded - npts_original))
    end

    nfreq       = length(periods)
    frequencies = 1.0 ./ periods

    # Output arrays
    filtered_traces  = zeros(Float64, npts_original, nfreq)
    envelopes        = zeros(Float64, npts_original, nfreq)
    arrival_times    = fill(NaN,     nfreq)
    group_velocities = fill(NaN,     nfreq)
    phase_velocities = fill(NaN,     nfreq)
    phase_branch_numbers = copy(DEFAULT_PHASE_BRANCH_NUMBERS)
    phase_velocity_branches = fill(NaN, nfreq, length(phase_branch_numbers))
    measured_phases  = fill(NaN,     nfreq)
    selected_phase_branches = fill(0, nfreq)
    phase_suspect    = falses(nfreq)
    u_predicted_from_phase = fill(NaN, nfreq)
    amplitudes       = zeros(Float64, nfreq)
    quality_factors  = zeros(Float64, nfreq)
    all_peaks        = [Tuple{Float64,Float64}[] for _ in 1:nfreq]
    time             = collect(range(dt, step=dt, length=npts_original))

    min_vel, max_vel = velocity_range

    for (i, freq) in enumerate(frequencies)
        # --- Complex analytic signal via one-sided Gaussian ---
        z = narrow_band_filter_analytic(data, dt, freq, bandwidth_factor)[1:npts_original]
        filtered_traces[:, i] = real.(z)

        # Envelope = modulus of analytic signal
        envelope = abs.(z)
        envelopes[:, i] = envelope

        # Velocity-bounded search window
        t_min = dist / max_vel
        t_max = dist / min_vel
        search_window = (t_min, t_max)

        # Up to four strongest peaks in the search window; keep the strongest
        # one as the canonical group-velocity pick for backward compatibility.
        peaks = find_group_arrivals(envelope, time, search_window; max_peaks=4)
        all_peaks[i] = peaks
        t_group, amp_peak = isempty(peaks) ? (NaN, 0.0) : first(peaks)
        arrival_times[i] = t_group

        if !isnan(t_group) && t_group > 0.0
            group_velocities[i] = dist / t_group
        end

        amplitudes[i] = isnan(t_group) ? maximum(envelope) : amp_peak

        # Quality factor: peak / mean (SNR proxy, same as sacmft96)
        mean_env = mean(envelope)
        quality_factors[i] = mean_env > 0.0 ? amplitudes[i] / mean_env : 0.0

        # Collect wrapped analytic phase at group arrival for candidate branches.
        if compute_phase && !isnan(t_group) && t_group > 0.0
            i_g = clamp(round(Int, (t_group - time[1]) / dt) + 1, 1, npts_original)
            measured_phases[i] = angle(z[i_g])
        end
    end

    if compute_phase
        finish_phase_velocity_picks!(phase_velocities, phase_velocity_branches,
                                     measured_phases, selected_phase_branches,
                                     phase_suspect, u_predicted_from_phase,
                                     periods, frequencies, arrival_times, dist,
                                     quality_factors, phase_branch_numbers,
                                     phvel_source_phase;
                                     min_anchor_quality=min_anchor_quality,
                                     phase_anchor_velocity=phase_anchor_velocity,
                                     phase_velocity_range=phase_velocity_range,
                                     phase_smoothness_jump=phase_smoothness_jump,
                                     phase_plausible_range=phase_plausible_range)
    end

    return MFTResult(periods, frequencies, group_velocities, phase_velocities,
                     phase_branch_numbers, phase_velocity_branches,
                     measured_phases, selected_phase_branches, phase_suspect,
                     u_predicted_from_phase,
                     amplitudes, filtered_traces, envelopes, arrival_times,
                     dist, quality_factors, all_peaks, time, :full)
end

	function perform_mft_analysis(data::AbstractVector{<:Real}, dt::Real, distance::Real,
	                              period_range, n_periods; kwargs...)
		period_min, period_max = Float64.(period_range)
		periods_analysis = collect(exp10.(range(log10(period_min), log10(period_max), length=n_periods)))
		return perform_mft_analysis(data, dt, distance, periods_analysis; kwargs...)
	end

function _require_batch_size(bank::MFTFilterBank, N::Int)
    N <= bank.N_buf || throw(ArgumentError(
        "MFTFilterBank was preallocated for $(bank.N_buf) waveform columns, " *
        "but this call needs $(N). Allocate a bank with N_initial >= $(N) " *
        "or use non-bang perform_mft_analysis_batch to allocate a call-local bank."
    ))
    return nothing
end

"""
    perform_mft_analysis!(bank, data, distance; compute_phase=true, phvel_correction=0.0)

Single-waveform MFT using a pre-computed `MFTFilterBank`.
"""
function perform_mft_analysis!(bank::MFTFilterBank,
                               data::AbstractVector{<:Real},
                               distance::Real;
                               compute_phase::Bool=true,
                               phvel_correction::Float64=0.0,
                               phvel_source_phase::Float64=phvel_correction,
                               kwargs...)
    W = reshape(data, :, 1)
    return only(perform_mft_analysis_batch!(bank, W, distance;
        compute_phase=compute_phase,
        phvel_source_phase=phvel_source_phase,
        kwargs...))
end

# ╔═╡ bdb2f320-c713-43aa-bb6b-98cca0606266
function _assemble_mft_result(bank::MFTFilterBank, dist::Float64, n::Int,
                               all_peaks_n::Vector{Vector{Tuple{Float64,Float64}}};
                               compute_phase::Bool=false,
                               use_phtovel::Bool=false,
                               phvel_correction::Float64=0.0,
                               phvel_source_phase::Float64=phvel_correction,
                               min_anchor_quality::Float64=3.0,
                               phase_anchor_velocity::Float64=3.3,
                               phase_velocity_range::Tuple{Float64,Float64}=(2.5, 4.5),
                               phase_smoothness_jump::Float64=0.05,
                               phase_plausible_range::Tuple{Float64,Float64}=(0.5, 20.0),
                               phtovel_prior_periods=nothing,
                               phtovel_prior_velocities=nothing)
    nfreq = bank.nfreq
    phase_branch_numbers = copy(DEFAULT_PHASE_BRANCH_NUMBERS)
    nbranches = length(phase_branch_numbers)

    arrivals     = bank.arrivals_buf[:, n]
    bank.storage_mode == :full ||
        throw(ArgumentError("_assemble_mft_result requires a :full MFTFilterBank"))
    filtered     = bank.filtered_buf[:, :, n]
    envelopes    = bank.envelopes_buf[:, :, n]

    group_velocities      = fill(NaN, nfreq)
    amplitudes            = zeros(Float64, nfreq)
    quality_factors       = zeros(Float64, nfreq)
    phase_velocities      = fill(NaN, nfreq)
    phase_velocity_branches = fill(NaN, nfreq, nbranches)
    measured_phases       = fill(NaN, nfreq)
    selected_phase_branches = fill(0, nfreq)
    phase_suspect         = falses(nfreq)
    u_predicted_from_phase = fill(NaN, nfreq)

    for ifreq in 1:nfreq
        t_g = arrivals[ifreq]
        envelope = @view envelopes[:, ifreq]
        peaks = all_peaks_n[ifreq]
        amp_peak = isempty(peaks) ? 0.0 : first(peaks)[2]

        if !isnan(t_g) && t_g > 0.0
            group_velocities[ifreq] = dist / t_g
        end
        amplitudes[ifreq]     = isnan(t_g) ? maximum(envelope) : amp_peak
        mean_env              = mean(envelope)
        quality_factors[ifreq] = mean_env > 0.0 ? amplitudes[ifreq] / mean_env : 0.0

        if compute_phase && !isnan(t_g) && t_g > 0.0
            measured_phases[ifreq] = bank.phases_buf[ifreq, n]
        end
    end

    if compute_phase
        if use_phtovel
            finish_phase_velocity_picks_phtovel!(phase_velocities, measured_phases,
                                                 selected_phase_branches,
                                                 phase_suspect, u_predicted_from_phase,
                                                 bank.periods, bank.frequencies, arrivals,
                                                 dist, phvel_source_phase;
                                                 phase_anchor_velocity=phase_anchor_velocity,
                                                 phase_velocity_range=phase_velocity_range,
                                                 phase_smoothness_jump=phase_smoothness_jump,
                                                 phase_plausible_range=phase_plausible_range,
                                                 phtovel_prior_periods=phtovel_prior_periods,
                                                 phtovel_prior_velocities=phtovel_prior_velocities)
        else
            finish_phase_velocity_picks!(phase_velocities, phase_velocity_branches,
                                         measured_phases, selected_phase_branches,
                                         phase_suspect, u_predicted_from_phase,
                                         bank.periods, bank.frequencies, arrivals, dist,
                                         quality_factors, phase_branch_numbers,
                                         phvel_source_phase;
                                         min_anchor_quality=min_anchor_quality,
                                         phase_anchor_velocity=phase_anchor_velocity,
                                         phase_velocity_range=phase_velocity_range,
                                         phase_smoothness_jump=phase_smoothness_jump,
                                         phase_plausible_range=phase_plausible_range)
        end
    end

    return MFTResult(
        bank.periods, bank.frequencies,
        group_velocities, phase_velocities,
        phase_branch_numbers, phase_velocity_branches,
        measured_phases, selected_phase_branches, phase_suspect,
        u_predicted_from_phase,
        amplitudes, filtered, envelopes, arrivals,
        dist, quality_factors, all_peaks_n, bank.time, :full,
    )
end

# ╔═╡ 9a6cb7e1-6ea0-4dd1-97f3-47f3b8a717de
# --- Phase 1: dist-independent (upsample, FFT, filter, IFFT → fills envelopes_buf/filtered_buf)
function _run_phase1!(bank::MFTFilterBank{T}, W_flat::AbstractMatrix{T}, N::Int) where {T}
    bank.storage_mode == :full ||
        throw(ArgumentError("_run_phase1! requires a :full MFTFilterBank"))
    W_up = bank.upsample_factor == 1.0 ? W_flat :
        DSP.resample(W_flat, bank.upsample_factor; dims=1)   # (npts_original × N)
    @views bank.W_pad_buf[:, 1:N] .= 0
    @views bank.W_pad_buf[1:bank.npts_original, 1:N] .= W_up
    mul!(bank.SPEC_buf, bank.plan_fwd, bank.W_pad_buf)
    for ifreq in 1:bank.nfreq
        @views @. bank.Z_buf[:, 1:N] = bank.H_full[:, ifreq] * bank.SPEC_buf[:, 1:N]
        mul!(bank.Z_time_buf, bank.plan_inv, bank.Z_buf)
        @views @. bank.filtered_buf[1:bank.npts_original, ifreq, 1:N]  =
            real(bank.Z_time_buf[1:bank.npts_original, 1:N])
        @views @. bank.envelopes_buf[1:bank.npts_original, ifreq, 1:N] =
            abs(bank.Z_time_buf[1:bank.npts_original, 1:N])
    end
end

# ╔═╡ 6a5c34ae-f358-40bd-8954-6d7e8fb0fa8f
begin
	# --- Phase 2: dist-dependent (peak finding + assembly → returns Vector{MFTResult})
	function _run_phase2!(bank::MFTFilterBank, dists_flat::AbstractVector{Float64}, N::Int;
	                      compute_phase::Bool=false,
	                      use_phtovel::Bool=false,
	                      phvel_correction::Float64=0.0,
	                      phvel_source_phase::Float64=phvel_correction,
	                      min_anchor_quality::Float64=3.0,
	                      phase_anchor_velocity::Float64=3.3,
	                      phase_velocity_range::Tuple{Float64,Float64}=(2.5, 4.5),
	                      phase_smoothness_jump::Float64=0.05,
	                      phase_plausible_range::Tuple{Float64,Float64}=(0.5, 20.0),
	                      phtovel_prior_periods=nothing,
	                      phtovel_prior_velocities=nothing)
	    bank.arrivals_buf[:, 1:N] .= NaN
	    bank.phases_buf[:, 1:N]   .= NaN
	    min_vel, max_vel = bank.velocity_range
	
	    all_peaks = [[Tuple{Float64,Float64}[] for _ in 1:bank.nfreq] for _ in 1:N]
	    for n in 1:N
	        dist  = dists_flat[n]
	        t_min = dist / max_vel
	        t_max = dist / min_vel
	        for ifreq in 1:bank.nfreq
	            if compute_phase
	                @views @. bank.Z_buf[:, 1:N] = bank.H_full[:, ifreq] * bank.SPEC_buf[:, 1:N]
	                mul!(bank.Z_time_buf, bank.plan_inv, bank.Z_buf)
	            end
	            env   = @view bank.envelopes_buf[:, ifreq, n]
	            peaks = find_group_arrivals(env, bank.time, (t_min, t_max); max_peaks=4)
	            all_peaks[n][ifreq] = peaks
	            t_g   = isempty(peaks) ? NaN : first(peaks)[1]
	            bank.arrivals_buf[ifreq, n] = t_g
	            if compute_phase && isfinite(t_g) && t_g > 0.0
	                i_g = clamp(round(Int, (t_g - bank.time[1]) / bank.dt) + 1, 1, bank.npts_original)
	                bank.phases_buf[ifreq, n] = Float64(angle(bank.Z_time_buf[i_g, n]))
	            end
	        end
	    end
	
	    results = Vector{MFTResult}(undef, N)
	    for n in 1:N
	        results[n] = _assemble_mft_result(bank, dists_flat[n], n, all_peaks[n];
	                                          compute_phase=compute_phase,
	                                          use_phtovel=use_phtovel,
	                                          phvel_source_phase=phvel_source_phase,
	                                          min_anchor_quality=min_anchor_quality,
	                                          phase_anchor_velocity=phase_anchor_velocity,
	                                          phase_velocity_range=phase_velocity_range,
	                                          phase_smoothness_jump=phase_smoothness_jump,
	                                          phase_plausible_range=phase_plausible_range,
	                                          phtovel_prior_periods=phtovel_prior_periods,
	                                          phtovel_prior_velocities=phtovel_prior_velocities)
	    end
	    return results
	end
	
	# Picks-only banks keep only one filtered period in memory at a time.  The
	# dispersion summaries are still Float64 so downstream tables and plots keep the
	# same public shape as full-storage results.
	function _empty_pick_result(bank::MFTFilterBank{T}, dist::Float64) where {T}
	    nfreq = bank.nfreq
	    phase_branch_numbers = copy(DEFAULT_PHASE_BRANCH_NUMBERS)
	    return MFTResult(
	        bank.periods, bank.frequencies,
	        fill(NaN, nfreq), fill(NaN, nfreq),
	        phase_branch_numbers, fill(NaN, nfreq, length(phase_branch_numbers)),
	        fill(NaN, nfreq), fill(0, nfreq), falses(nfreq), fill(NaN, nfreq),
	        zeros(Float64, nfreq), zeros(T, 0, 0), zeros(T, 0, 0),
	        fill(NaN, nfreq), dist, zeros(Float64, nfreq),
	        [Tuple{Float64,Float64}[] for _ in 1:nfreq], bank.time, :picks_only,
	    )
	end
	
	function _finish_phase_velocities!(res::MFTResult, measured_phases::Vector{Float64},
	                                   phvel_source_phase::Float64;
	                                   use_phtovel::Bool=false,
	                                   min_anchor_quality::Float64=3.0,
	                                   phase_anchor_velocity::Float64=3.3,
	                                   phase_velocity_range::Tuple{Float64,Float64}=(2.5, 4.5),
	                                   phase_smoothness_jump::Float64=0.05,
	                                   phase_plausible_range::Tuple{Float64,Float64}=(0.5, 20.0),
	                                   phtovel_prior_periods=nothing,
	                                   phtovel_prior_velocities=nothing)
	    res.measured_phases .= measured_phases
	    if use_phtovel
	        finish_phase_velocity_picks_phtovel!(res.phase_velocities, res.measured_phases,
	                                             res.selected_phase_branches,
	                                             res.phase_suspect, res.u_predicted_from_phase,
	                                             res.periods, res.frequencies, res.arrival_times,
	                                             res.distance, phvel_source_phase;
	                                             phase_anchor_velocity=phase_anchor_velocity,
	                                             phase_velocity_range=phase_velocity_range,
	                                             phase_smoothness_jump=phase_smoothness_jump,
	                                             phase_plausible_range=phase_plausible_range,
	                                             phtovel_prior_periods=phtovel_prior_periods,
	                                             phtovel_prior_velocities=phtovel_prior_velocities)
	    else
	        finish_phase_velocity_picks!(res.phase_velocities, res.phase_velocity_branches,
	                                     res.measured_phases, res.selected_phase_branches,
	                                     res.phase_suspect, res.u_predicted_from_phase,
	                                     res.periods, res.frequencies, res.arrival_times,
	                                     res.distance, res.quality_factors,
	                                     res.phase_branch_numbers, phvel_source_phase;
	                                     min_anchor_quality=min_anchor_quality,
	                                     phase_anchor_velocity=phase_anchor_velocity,
	                                     phase_velocity_range=phase_velocity_range,
	                                     phase_smoothness_jump=phase_smoothness_jump,
	                                     phase_plausible_range=phase_plausible_range)
	    end
	    return res
	end
	
	function _zero_lag_correlation_vectors(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
	    a0 = a .- mean(a)
	    b0 = b .- mean(b)
	    den = sqrt(sum(abs2, a0) * sum(abs2, b0))
	    den > 0 || return NaN
	    return clamp(Float64(sum(a0 .* b0) / den), -1.0, 1.0)
	end
	
	function _prepare_fft_inputs!(bank::MFTFilterBank{T}, W_flat::AbstractMatrix{T}, N::Int) where {T}
	    W_up = bank.upsample_factor == 1.0 ? W_flat :
	        DSP.resample(W_flat, bank.upsample_factor; dims=1)
	    @views bank.W_pad_buf[:, 1:N] .= 0
	    @views bank.W_pad_buf[1:bank.npts_original, 1:N] .= W_up
	    mul!(bank.SPEC_buf, bank.plan_fwd, bank.W_pad_buf)
	    return nothing
	end
	
	function _run_picks_only!(bank::MFTFilterBank{T}, W_flat::AbstractMatrix{T},
	                          dists_flat::AbstractVector{Float64}, N::Int;
	                          compute_phase::Bool=false,
	                          use_phtovel::Bool=false,
	                          phvel_correction::Float64=0.0,
	                          phvel_source_phase::Float64=phvel_correction,
	                          min_anchor_quality::Float64=3.0,
	                          phase_anchor_velocity::Float64=3.3,
	                          phase_velocity_range::Tuple{Float64,Float64}=(2.5, 4.5),
	                          phase_smoothness_jump::Float64=0.05,
	                          phase_plausible_range::Tuple{Float64,Float64}=(0.5, 20.0),
	                          phtovel_prior_periods=nothing,
	                          phtovel_prior_velocities=nothing,
	                          correlation_pairs::Union{Nothing,Vector{Tuple{Int,Int}}}=nothing) where {T}
	    bank.storage_mode == :picks_only ||
	        throw(ArgumentError("_run_picks_only! requires a :picks_only MFTFilterBank"))
	    _prepare_fft_inputs!(bank, W_flat, N)
	
	    min_vel, max_vel = bank.velocity_range
	    results = [_empty_pick_result(bank, dists_flat[n]) for n in 1:N]
	    measured_phases = [fill(NaN, bank.nfreq) for _ in 1:N]
	    correlations = isnothing(correlation_pairs) ? nothing :
	        fill(NaN, bank.nfreq, length(correlation_pairs))
	
	    for ifreq in 1:bank.nfreq
	        @views @. bank.Z_buf[:, 1:N] = bank.H_full[:, ifreq] * bank.SPEC_buf[:, 1:N]
	        mul!(bank.Z_time_buf, bank.plan_inv, bank.Z_buf)
	
	        if !isnothing(correlations)
	            for (ipair, (ia, ib)) in enumerate(correlation_pairs)
	                correlations[ifreq, ipair] = _zero_lag_correlation_vectors(
	                    real.(@view bank.Z_time_buf[1:bank.npts_original, ia]),
	                    real.(@view bank.Z_time_buf[1:bank.npts_original, ib]),
	                )
	            end
	        end
	
	        for n in 1:N
	            res = results[n]
	            envelope = abs.(@view bank.Z_time_buf[1:bank.npts_original, n])
	            dist = res.distance
	            peaks = find_group_arrivals(envelope, bank.time,
	                                        (dist / max_vel, dist / min_vel); max_peaks=4)
	            res.all_peaks[ifreq] = peaks
	            t_g, amp_peak = isempty(peaks) ? (NaN, 0.0) : first(peaks)
	            res.arrival_times[ifreq] = t_g
	            if isfinite(t_g) && t_g > 0.0
	                res.group_velocities[ifreq] = dist / t_g
	            end
	            res.amplitudes[ifreq] = isnan(t_g) ? Float64(maximum(envelope)) : amp_peak
	            mean_env = mean(envelope)
	            res.quality_factors[ifreq] =
	                mean_env > 0.0 ? Float64(res.amplitudes[ifreq] / mean_env) : 0.0
	
	            if compute_phase && isfinite(t_g) && t_g > 0.0
	                i_g = clamp(round(Int, (t_g - bank.time[1]) / bank.dt) + 1, 1, bank.npts_original)
	                measured_phases[n][ifreq] = Float64(angle(bank.Z_time_buf[i_g, n]))
	            end
	        end
	    end
	
	    if compute_phase
	        for n in eachindex(results)
	            _finish_phase_velocities!(results[n], measured_phases[n], phvel_source_phase;
	                                       use_phtovel=use_phtovel,
	                                       min_anchor_quality=min_anchor_quality,
	                                       phase_anchor_velocity=phase_anchor_velocity,
	                                       phase_velocity_range=phase_velocity_range,
	                                       phase_smoothness_jump=phase_smoothness_jump,
	                                       phase_plausible_range=phase_plausible_range,
	                                       phtovel_prior_periods=phtovel_prior_periods,
	                                       phtovel_prior_velocities=phtovel_prior_velocities)
	        end
	    end
	    return results, correlations
	end
end

# ╔═╡ 0a2beec9-fe91-4af8-8646-064ec215c1e6
"""
    perform_mft_analysis_batch!(bank, W, dist; compute_phase=false) -> Vector{MFTResult}
    perform_mft_analysis_batch!(bank, W, dists; compute_phase=false) -> Array{MFTResult}

Run MFT on a batch of waveforms using a pre-computed `MFTFilterBank`. The bang
form never grows or replans `bank`; allocate it with `N_initial` large enough
for the number of waveform columns in the call.

**2-D form** — `W::AbstractMatrix{Float64}` of shape `(nt × N)`, `dist::Float64`:
Returns `Vector{MFTResult}` of length N. One shared distance for all waveforms.

**N-D form** — `W::AbstractArray{Float64}` of shape `(nt × d1 × d2 × ...)`,
`dists::AbstractArray{Float64}` of shape `(d1 × d2 × ...)`:
Returns `Array{MFTResult}` with the same shape `(d1 × d2 × ...)`.
Each position uses its own distance from `dists`. All waveforms are FFT'd in one
batch pass (Phase 1), then peak-finding and assembly use the per-position distance
(Phase 2). `reshape` is zero-copy so no data is duplicated.

Scalar-dist convenience: passing `dist::Float64` to the N-D form broadcasts to
`fill(dist, size(W)[2:end])` — equivalent to the 2-D form.
"""
function perform_mft_analysis_batch!(bank::MFTFilterBank,
                                     W::AbstractMatrix{<:Real},
                                     dist::Real;
                                     compute_phase::Bool=false,
                                     use_phtovel::Bool=false,
                                     phvel_correction::Float64=0.0,
                                     phvel_source_phase::Float64=phvel_correction,
                                     min_anchor_quality::Float64=3.0,
                                     phase_anchor_velocity::Float64=3.3,
                                     phase_velocity_range::Tuple{Float64,Float64}=(2.5, 4.5),
                                     phase_smoothness_jump::Float64=0.05,
                                     phase_plausible_range::Tuple{Float64,Float64}=(0.5, 20.0),
                                     phtovel_prior_periods=nothing,
                                     phtovel_prior_velocities=nothing)
    N = size(W, 2)
    _require_batch_size(bank, N)
    T = _bank_float_type(bank)
    W_flat = T.(W)
    dists = fill(Float64(dist), N)
    if bank.storage_mode == :picks_only
        results, _ = _run_picks_only!(bank, W_flat, dists, N;
                                      compute_phase=compute_phase,
                                      use_phtovel=use_phtovel,
                                      phvel_source_phase=phvel_source_phase,
                                      min_anchor_quality=min_anchor_quality,
                                      phase_anchor_velocity=phase_anchor_velocity,
                                      phase_velocity_range=phase_velocity_range,
                                      phase_smoothness_jump=phase_smoothness_jump,
                                      phase_plausible_range=phase_plausible_range,
                                      phtovel_prior_periods=phtovel_prior_periods,
                                      phtovel_prior_velocities=phtovel_prior_velocities)
        return results
    end
    _run_phase1!(bank, W_flat, N)
    return _run_phase2!(bank, dists, N;
                        compute_phase=compute_phase,
                        use_phtovel=use_phtovel,
                        phvel_source_phase=phvel_source_phase,
                        min_anchor_quality=min_anchor_quality,
                        phase_anchor_velocity=phase_anchor_velocity,
                        phase_velocity_range=phase_velocity_range,
                        phase_smoothness_jump=phase_smoothness_jump,
                        phase_plausible_range=phase_plausible_range,
                        phtovel_prior_periods=phtovel_prior_periods,
                        phtovel_prior_velocities=phtovel_prior_velocities)
end

# ╔═╡ f57cbf2e-5dcb-4e0a-90a5-5f2c363990cd
function perform_mft_analysis_batch!(bank::MFTFilterBank,
                                     W::AbstractArray{<:Real},
                                     dists::AbstractArray{<:Real};
                                     compute_phase::Bool=false,
                                     use_phtovel::Bool=false,
                                     phvel_correction::Float64=0.0,
                                     phvel_source_phase::Float64=phvel_correction,
                                     min_anchor_quality::Float64=3.0,
                                     phase_anchor_velocity::Float64=3.3,
                                     phase_velocity_range::Tuple{Float64,Float64}=(2.5, 4.5),
                                     phase_smoothness_jump::Float64=0.05,
                                     phase_plausible_range::Tuple{Float64,Float64}=(0.5, 20.0),
                                     phtovel_prior_periods=nothing,
                                     phtovel_prior_velocities=nothing)
    dims = size(W)[2:end]
    @assert size(dists) == dims "dists shape $(size(dists)) must match trailing dims of W $(dims)"
    N = prod(dims)
    _require_batch_size(bank, N)
    T = _bank_float_type(bank)
    W_flat = reshape(T.(W), size(W, 1), N)
    dists_flat = vec(Float64.(dists))
    results_flat = if bank.storage_mode == :picks_only
        first(_run_picks_only!(bank, W_flat, dists_flat, N;
                               compute_phase=compute_phase,
                               use_phtovel=use_phtovel,
                               phvel_source_phase=phvel_source_phase,
                               min_anchor_quality=min_anchor_quality,
                               phase_anchor_velocity=phase_anchor_velocity,
                               phase_velocity_range=phase_velocity_range,
                               phase_smoothness_jump=phase_smoothness_jump,
                               phase_plausible_range=phase_plausible_range,
                               phtovel_prior_periods=phtovel_prior_periods,
                               phtovel_prior_velocities=phtovel_prior_velocities))
    else
        _run_phase1!(bank, W_flat, N)
        _run_phase2!(bank, dists_flat, N;
                     compute_phase=compute_phase,
                     use_phtovel=use_phtovel,
                     phvel_source_phase=phvel_source_phase,
                     min_anchor_quality=min_anchor_quality,
                     phase_anchor_velocity=phase_anchor_velocity,
                     phase_velocity_range=phase_velocity_range,
                     phase_smoothness_jump=phase_smoothness_jump,
                     phase_plausible_range=phase_plausible_range,
                     phtovel_prior_periods=phtovel_prior_periods,
                     phtovel_prior_velocities=phtovel_prior_velocities)
    end
    return reshape(results_flat, dims)
end

# ╔═╡ f310a7fd-cab6-4d92-b54e-060e113ab2f1
begin
	# Scalar-dist convenience for N-D arrays (broadcasts one distance to all positions)
	function perform_mft_analysis_batch!(bank::MFTFilterBank,
	                                     W::AbstractArray{<:Real},
	                                     dist::Real;
	                                     kwargs...)
	    dims = size(W)[2:end]
	    return perform_mft_analysis_batch!(bank, W, fill(Float64(dist), dims); kwargs...)
	end
	
	"""
	    perform_mft_analysis_batch(W, dt, dist, periods; ...) -> Vector{MFTResult}
	    perform_mft_analysis_batch(W, dt, dists, periods; ...) -> Array{MFTResult}
	
	Allocating batch MFT entrypoint. A call-local `MFTFilterBank` is sized for the
	waveform columns in `W`, then passed to `perform_mft_analysis_batch!`.
	"""
	function perform_mft_analysis_batch(W::AbstractArray{<:Real},
	                                    dt::Real,
	                                    dists::Union{Real,AbstractArray{<:Real}},
	                                    periods::Vector{Float64};
	                                    velocity_range::Tuple{Float64,Float64}=(2.0, 6.0),
	                                    bandwidth_factor::Float64=sqrt(2.0 / 25.0),
	                                    zero_pad_factor::Int=4,
	                                    upsample_factor::Real=2.0,
	                                    precision::Type{<:AbstractFloat}=Float32,
	                                    storage_mode::Symbol=:picks_only,
	                                    compute_phase::Bool=false,
	                                    use_phtovel::Bool=false,
	                                    phvel_correction::Float64=0.0,
	                                    phvel_source_phase::Float64=phvel_correction,
	                                    min_anchor_quality::Float64=3.0,
	                                    phase_anchor_velocity::Float64=3.3,
	                                    phase_velocity_range::Tuple{Float64,Float64}=(2.5, 4.5),
	                                    phase_smoothness_jump::Float64=0.05,
	                                    phase_plausible_range::Tuple{Float64,Float64}=(0.5, 20.0),
	                                    phtovel_prior_periods=nothing,
	                                    phtovel_prior_velocities=nothing)
	    ndims(W) >= 2 || throw(ArgumentError("W must have shape (nt × trailing dims...)"))
	    N = prod(size(W)[2:end])
	    bank = MFTFilterBank(dt, size(W, 1), periods;
	                         velocity_range=velocity_range,
	                         bandwidth_factor=bandwidth_factor,
	                         zero_pad_factor=zero_pad_factor,
	                         upsample_factor=upsample_factor,
	                         precision=precision,
	                         storage_mode=storage_mode,
	                         N_initial=N)
	    return perform_mft_analysis_batch!(bank, W, dists;
	                                       compute_phase=compute_phase,
	                                       use_phtovel=use_phtovel,
	                                       phvel_source_phase=phvel_source_phase,
	                                       min_anchor_quality=min_anchor_quality,
	                                       phase_anchor_velocity=phase_anchor_velocity,
	                                       phase_velocity_range=phase_velocity_range,
	                                       phase_smoothness_jump=phase_smoothness_jump,
	                                       phase_plausible_range=phase_plausible_range,
	                                       phtovel_prior_periods=phtovel_prior_periods,
	                                       phtovel_prior_velocities=phtovel_prior_velocities)
	end
end
