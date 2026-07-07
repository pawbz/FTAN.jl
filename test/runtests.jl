using FTAN
using Test

function _synthetic_waveform(; dt=0.25, n=2048, distance=120.0)
    t = collect((1:n) .* dt)
    data = @. exp(-0.5 * ((t - distance / 3.2) / 7.0)^2) * cos(2π * 0.18 * t)
    data .+= @. 0.35 * exp(-0.5 * ((t - distance / 2.6) / 9.0)^2) * cos(2π * 0.10 * t)
    return Float64.(data), dt, distance
end

@testset "FTAN array-first MFT" begin
    data, dt, distance = _synthetic_waveform()
    periods = collect(range(4.0, 14.0, length=12))

    res = perform_mft_analysis(data, dt, distance, periods;
        velocity_range=(2.0, 5.0),
        bandwidth_factor=0.30,
        compute_phase=true,
        storage_mode=:full)
    @test res isa MFTResult
    @test length(res.periods) == length(periods)
    @test any(isfinite, res.group_velocities)

    bank = MFTFilterBank(dt, length(data), periods;
        velocity_range=(2.0, 5.0),
        bandwidth_factor=0.30,
        storage_mode=:picks_only,
        N_initial=2)
    res_bank = perform_mft_analysis!(bank, data, distance; compute_phase=true)
    @test res_bank isa MFTResult
    @test res_bank.storage_mode == :picks_only

    W = hcat(data, 0.8 .* data)
    batch = perform_mft_analysis_batch!(bank, W, distance; compute_phase=true)
    @test length(batch) == 2

    plot_obj = plot_uc_consistency_comparison(res; pair_label="synthetic")
    @test !isnothing(plot_obj)

    res_phtovel = perform_mft_analysis(data, dt, distance, periods;
        velocity_range=(2.0, 5.0),
        bandwidth_factor=0.30,
        compute_phase=true,
        use_phtovel=true,
        phase_velocity_range=(1.5, 6.0))
    @test res_phtovel isa MFTResult
end

@testset "phvel_source_phase period-dependent correction" begin
    # _interp_phvel_correction: exact-match keys return the stored value;
    # in-between periods interpolate linearly; out-of-range periods clamp.
    d = Dict(5.0 => -π/4, 10.0 => -π/8, 20.0 => -π/16)
    @test FTAN._interp_phvel_correction(5.0, d) ≈ -π/4
    @test FTAN._interp_phvel_correction(10.0, d) ≈ -π/8
    @test FTAN._interp_phvel_correction(20.0, d) ≈ -π/16
    @test FTAN._interp_phvel_correction(7.5, d) ≈ (-π/4 + -π/8) / 2
    @test FTAN._interp_phvel_correction(1.0, d) ≈ -π/4   # clamped below range
    @test FTAN._interp_phvel_correction(100.0, d) ≈ -π/16 # clamped above range

    # Single-entry dict behaves like a constant scalar at every period.
    d1 = Dict(9.0 => -π/6)
    @test FTAN._interp_phvel_correction(1.0, d1) ≈ -π/6
    @test FTAN._interp_phvel_correction(9.0, d1) ≈ -π/6
    @test FTAN._interp_phvel_correction(50.0, d1) ≈ -π/6

    @test_throws ArgumentError FTAN._phvel_correction_per_period(Dict{Float64,Float64}(), [5.0, 10.0])

    periods = [5.0, 7.5, 10.0, 20.0]
    @test FTAN._phvel_correction_per_period(0.3, periods) == fill(0.3, length(periods))
    @test FTAN._phvel_correction_per_period(d, periods) ≈
        [FTAN._interp_phvel_correction(T, d) for T in periods]

    # End-to-end: a dict correction changes phase-velocity picks relative to
    # a 0.0 scalar, for both the legacy branch-matrix and phtovel paths.
    data, dt, distance = _synthetic_waveform()
    periods = collect(range(4.0, 14.0, length=12))
    correction = Dict(4.0 => π/6, 9.0 => π/8, 14.0 => π/10)

    for use_phtovel in (false, true)
        res0 = perform_mft_analysis(data, dt, distance, periods;
            velocity_range=(2.0, 5.0), bandwidth_factor=0.30,
            compute_phase=true, use_phtovel=use_phtovel,
            phase_velocity_range=(1.5, 6.0))
        res_corr = perform_mft_analysis(data, dt, distance, periods;
            velocity_range=(2.0, 5.0), bandwidth_factor=0.30,
            compute_phase=true, use_phtovel=use_phtovel,
            phase_velocity_range=(1.5, 6.0), phvel_source_phase=correction)
        @test res_corr isa MFTResult
        finite = findall(i -> isfinite(res0.phase_velocities[i]) && isfinite(res_corr.phase_velocities[i]),
                          eachindex(periods))
        @test !isempty(finite)
        @test any(i -> !isapprox(res0.phase_velocities[i], res_corr.phase_velocities[i]), finite)
    end
end
