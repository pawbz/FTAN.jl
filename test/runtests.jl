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

    br = analyze_causal_acausal_branches(data, reverse(data), dt, distance, periods;
        velocity_range=(2.0, 5.0),
        bandwidth_factor=0.30,
        compute_phase=true)
    @test br isa BranchAnalysisResult

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
