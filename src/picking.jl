function find_group_arrivals(envelope::AbstractVector{<:Real}, time::AbstractVector{Float64},
                             search_window::Tuple{Float64,Float64};
                             max_peaks::Int=4)
    t_start, t_end = search_window
    idx_window = findall(t -> t_start <= t <= t_end, time)
    isempty(idx_window) && return Tuple{Float64,Float64}[]

    env_win = envelope[idx_window]

    # All strict local maxima in the windowed envelope
    peak_idx, peak_vals = findmaxima(env_win)

    if isempty(peak_idx)
        # Fall back: absolute maximum (flat or monotone envelope in window)
        i_best = argmax(env_win)
        return [(time[idx_window[i_best]], env_win[i_best])]
    end

    order = sortperm(peak_vals, rev=true)
    nkeep = min(max_peaks, length(order))
    peaks = Tuple{Float64,Float64}[]
    for j in 1:nkeep
        idx = peak_idx[order[j]]
        push!(peaks, (time[idx_window[idx]], env_win[idx]))
    end
    return peaks
end

# ╔═╡ 74f5430a-8d8a-41cf-b683-52cc29bb12c0
"""
    find_group_arrival(envelope, time, search_window) -> (t_group, amplitude)

Return the strongest envelope peak in the search window. This preserves the
original group-velocity pick behaviour while `find_group_arrivals` exposes the
additional overtone candidates.
"""
function find_group_arrival(envelope::AbstractVector{<:Real}, time::Vector{Float64},
                            search_window::Tuple{Float64,Float64})
    peaks = find_group_arrivals(envelope, time, search_window; max_peaks=1)
    isempty(peaks) && return (NaN, 0.0)
    return first(peaks)
end

# ╔═╡ 82c85400-3cdb-11f1-9af8-1da7c901570c
"""
    extract_all_peaks_matrix(res::MFTResult, distance::Float64)

Convert sparse all_peaks storage into dense matrices for multipeak analysis.

Returns: `(periods, peak_times_matrix, peak_amps_matrix)`
- `peak_times_matrix[i, j]` — arrival time (s) of j-th peak at i-th period, or NaN if peak not found
- `peak_amps_matrix[i, j]` — amplitude of j-th peak, or NaN if peak not found
"""
function extract_all_peaks_matrix(res::MFTResult, distance::Float64)
    nperiods = length(res.periods)
    
    # Find maximum number of peaks across all periods
    max_peaks = maximum(length(peaks) for peaks in res.all_peaks; init=0)
    max_peaks = max(max_peaks, 1)  # Ensure at least 1
    
    # Initialize matrices with NaN
    peak_times = fill(NaN, nperiods, max_peaks)
    peak_amps = fill(NaN, nperiods, max_peaks)
    
    # Populate matrices from sparse storage
    for (i, peaks) in enumerate(res.all_peaks)
        for (j, (t_arrival, amplitude)) in enumerate(peaks)
            if j <= max_peaks
                peak_times[i, j] = t_arrival
                peak_amps[i, j] = amplitude
            end
        end
    end
    
    return (res.periods, peak_times, peak_amps)
end

# ╔═╡ 82c857c8-3cdb-11f1-8b14-3f44a70d8c9a
"""
    match_peaks_continuous(periods::Vector{Float64}, peak_times::Matrix{Float64},
                          peak_amps::Matrix{Float64}; velocity_range::Tuple{Float64,Float64}=(2.0, 6.0),
                          distance::Float64=1.0, start_peak_idx::Int=1) -> MultimodalDispersion

Extract one continuous dispersion curve via greedy period-to-period peak matching.

**Algorithm:** For each period, select the peak closest (in time) to the previous period's pick.
This minimizes period-to-period velocity jumps, producing a smooth single-mode dispersion curve.

# Arguments
- `periods` : Analysis periods [s]
- `peak_times` : Matrix of arrival times [period × peak], NaN for missing peaks
- `peak_amps` : Matrix of peak amplitudes [period × peak], for quality tracking
- `velocity_range` : Expected group velocity bounds for validation (currently unused; for validation only)
- `distance` : Source-receiver distance [km] for group velocity computation
- `start_peak_idx` : Which peak (column) to start from (default: 1 = strongest peak)

# Returns
- `MultimodalDispersion` — Continuous mode with periods, arrival times, velocities, amplitudes
"""
function match_peaks_continuous(periods::Vector{Float64}, peak_times::Matrix{Float64},
                               peak_amps::Matrix{Float64}; 
                               velocity_range::Tuple{Float64,Float64}=(2.0, 6.0),
                               distance::Float64=1.0,
                               start_peak_idx::Int=1)
    nperiods = length(periods)
    
    arrival_times = fill(NaN, nperiods)
    amplitudes = fill(NaN, nperiods)
    
    # Start with specified peak from first period
    if start_peak_idx <= size(peak_times, 2) && !isnan(peak_times[1, start_peak_idx])
        arrival_times[1] = peak_times[1, start_peak_idx]
        amplitudes[1] = peak_amps[1, start_peak_idx]
    else
        # Fallback: use first non-NaN peak in first period
        for j in 1:size(peak_times, 2)
            if !isnan(peak_times[1, j])
                arrival_times[1] = peak_times[1, j]
                amplitudes[1] = peak_amps[1, j]
                break
            end
        end
    end
    
    # Return NaN if first period has no valid peaks
    if isnan(arrival_times[1])
        group_vels = fill(NaN, nperiods)
        return MultimodalDispersion(periods, arrival_times, group_vels, amplitudes, 0)
    end
    
    # Greedy matching for subsequent periods
    for i in 2:nperiods
        prev_time = arrival_times[i-1]
        
        # Find closest non-NaN peak in current period
        best_j = 0
        best_dist = Inf
        
        for j in 1:size(peak_times, 2)
            t_candidate = peak_times[i, j]
            if !isnan(t_candidate)
                dist = abs(t_candidate - prev_time)
                if dist < best_dist
                    best_dist = dist
                    best_j = j
                end
            end
        end
        
        if best_j > 0
            arrival_times[i] = peak_times[i, best_j]
            amplitudes[i] = peak_amps[i, best_j]
        end
    end
    
    # Compute group velocities: v_g = distance / arrival_time
    group_vels = similar(arrival_times)
    for i in 1:nperiods
        if !isnan(arrival_times[i]) && arrival_times[i] > 0
            group_vels[i] = distance / arrival_times[i]
        else
            group_vels[i] = NaN
        end
    end
    
    return MultimodalDispersion(periods, arrival_times, group_vels, amplitudes, 1)
end

# ╔═╡ 82c86256-3cdb-11f1-befc-c7b40eeb6b06
"""
    extract_all_modes(res::MFTResult; distance::Float64=1.0, max_modes::Int=4,
                     velocity_range::Tuple{Float64,Float64}=(2.0, 6.0)) -> Vector{MultimodalDispersion}

Extract all continuous dispersion modes from MFT result.

**Strategy:** Iteratively extract modes by greedy peak matching, starting from each unused initial peak.
Marks peaks as "used" to prevent duplicate mode extraction.

# Arguments
- `res` : MFTResult from perform_mft_analysis
- `distance` : Source-receiver distance [km]
- `max_modes` : Maximum number of modes to extract (default: 4)
- `velocity_range` : Velocity bounds for validation (default: 2–6 km/s)

# Returns
- `Vector{MultimodalDispersion}` — Extracted modes, sorted by initial peak amplitude (strongest first)
"""
function extract_all_modes(res::MFTResult; distance::Float64=1.0, max_modes::Int=4,
                          velocity_range::Tuple{Float64,Float64}=(2.0, 6.0))
    periods, peak_times, peak_amps = extract_all_peaks_matrix(res, distance)
    
    modes = MultimodalDispersion[]
    used_peaks = falses(size(peak_times))
    
    for mode_idx in 1:max_modes
        # Find next unused peak in first period
        start_peak_idx = 0
        for j in 1:size(peak_times, 2)
            if !used_peaks[1, j] && !isnan(peak_times[1, j])
                start_peak_idx = j
                break
            end
        end
        
        # No more peaks to extract
        if start_peak_idx == 0
            break
        end
        
        # Extract continuous mode starting from this peak
        mode = match_peaks_continuous(periods, peak_times, peak_amps; 
                                       distance=distance, 
                                       velocity_range=velocity_range,
                                       start_peak_idx=start_peak_idx)
        
        # Mark used peaks
        for (i, t_arrival) in enumerate(mode.arrival_times)
            if !isnan(t_arrival)
                # Find which column was used
                for j in 1:size(peak_times, 2)
                    if !isnan(peak_times[i, j]) && abs(peak_times[i, j] - t_arrival) < 1e-6
                        used_peaks[i, j] = true
                        break
                    end
                end
            end
        end
        
        # Assign mode index and append
        mode = MultimodalDispersion(mode.periods, mode.arrival_times, mode.group_velocities, 
                                   mode.peak_amplitudes, length(modes) + 1)
        push!(modes, mode)
    end
    
    return modes
end
