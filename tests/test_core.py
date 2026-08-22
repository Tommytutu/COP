from __future__ import annotations

import numpy as np

from ahpcop import solve_mnvdm, validate_pcm
from ahpcop.evrim import solve_evrim, solve_evrim_check_first
from ahpcop.metrics import recovery_metrics, violation_score
from ahpcop.pcm import contains_directed_cycle, generate_pcm, llsm_weights
from ahpcop.priority import (
    GurobiSettings,
    Stage1Result,
    alpha_grid,
    solve_alpha_mnvdm,
    solve_mnvem,
    solve_mnvllsm,
    solve_stage1,
)
from ahpcop.sensitivity import _ordered_weak_orders


SETTINGS = GurobiSettings(time_limit=30, threads=1, output_flag=0)


def test_generators_respect_cr_regimes_and_cycle() -> None:
    bounds = {"low": (0.0, 0.05), "moderate": (0.05, 0.10), "high": (0.10, 0.20)}
    for regime, (lo, hi) in bounds.items():
        sample = generate_pcm(4, regime, 0, 20260815)
        assert sample.cr <= hi + 1e-10
        assert sample.cr >= lo - 1e-10 if lo == 0 else sample.cr > lo
        assert np.allclose(sample.matrix * sample.matrix.T, 1.0)
    cyclic = generate_pcm(4, "cyclic", 0, 20260815)
    assert cyclic.cycle is not None
    assert contains_directed_cycle(cyclic.matrix, cyclic.cycle)


def test_recovery_metrics_are_exact_at_ground_truth() -> None:
    sample = generate_pcm(4, "low", 1, 20260815)
    metrics = recovery_metrics(sample.matrix, sample.latent_weights, sample.latent_weights)
    assert metrics["kendall_tau_b"] == 1.0
    assert metrics["best_choice_accuracy"] == 1
    assert metrics["lrmse"] < 1e-12


def test_strong_and_basic_stage1_agree() -> None:
    sample = generate_pcm(4, "moderate", 2, 20260815)
    strong = solve_stage1(sample.matrix, SETTINGS, "strong")
    basic = solve_stage1(sample.matrix, SETTINGS, "basic")
    assert strong.solved and basic.solved
    assert strong.nv2 == basic.nv2


def test_cop_llsm_sanity_when_zero_violation_is_feasible() -> None:
    z = np.array([0.9, 0.3, -0.2, -1.0])
    consistent = np.exp(z[:, None] - z[None, :])
    stage1 = solve_stage1(consistent, SETTINGS, "strong")
    assert stage1.solved and stage1.nv == 0
    result = solve_mnvllsm(consistent, stage1, SETTINGS)
    assert result.solved and result.weights is not None
    assert violation_score(consistent, result.weights) == 0
    assert np.allclose(result.weights, llsm_weights(consistent), atol=2e-4)


def test_evrim_returns_an_order_representable_pcm() -> None:
    sample = generate_pcm(3, "cyclic", 0, 20260815)
    result = solve_evrim(sample.matrix, SETTINGS)
    assert result.solved and result.revised_matrix is not None
    assert result.nv == 0


def test_direct_and_oa_evrim_certificates_agree() -> None:
    sample = generate_pcm(3, "high", 1, 20260815)
    direct = solve_evrim(sample.matrix, SETTINGS, backend="direct")
    oa = solve_evrim(sample.matrix, SETTINGS, backend="oa_callback")
    assert direct.solved and oa.solved
    assert direct.nrp == oa.nrp
    assert abs(float(direct.aoc) - float(oa.aoc)) <= 1e-6
    assert direct.gci <= 0.31 + 1e-7
    assert oa.gci <= 0.31 + 1e-7


def test_stage1_preserves_repeated_intensity_equalities() -> None:
    a = np.array([
        [1, 3, 3, 1],
        [1 / 3, 1, 1, 1 / 3],
        [1 / 3, 1, 1, 1 / 3],
        [1, 3, 3, 1],
    ], dtype=float)
    result = solve_stage1(a, SETTINGS, "indicator")
    assert result.solved and result.nv == 0


def test_stage2_is_not_run_after_uncertified_stage1_timeout() -> None:
    timed_out = Stage1Result(
        status="TIME_LIMIT", solved=False, nv=2.0, nv2=4, runtime=30.0,
        gap=0.1, node_count=100.0, signs=[1], y=np.zeros(3), model_variant="strong",
    )
    a = np.ones((3, 3))
    for result in (solve_mnvllsm(a, timed_out, SETTINGS), solve_mnvem(a, timed_out, SETTINGS)):
        assert not result.solved
        assert result.status == "TIME_LIMIT"
        assert result.weights is None
        assert result.gap == 0.1


def test_public_single_matrix_mnvdm_api() -> None:
    a = np.array([
        [1, 2, 4, 9],
        [1 / 2, 1, 3, 7],
        [1 / 4, 1 / 3, 1, 5],
        [1 / 9, 1 / 7, 1 / 5, 1],
    ])
    assert np.array_equal(validate_pcm(a), a)
    result = solve_mnvdm(a, method="LLSM", settings=SETTINGS)
    assert result.solved
    assert result.weights is not None
    assert result.nv_star == 0
    assert np.isclose(result.weights.sum(), 1.0)


def test_public_api_rejects_nonreciprocal_matrix() -> None:
    bad = np.array([[1, 2, 3], [0.4, 1, 2], [1 / 3, 1 / 2, 1]])
    try:
        validate_pcm(bad)
    except ValueError as error:
        assert "reciprocal" in str(error)
    else:
        raise AssertionError("validate_pcm accepted a nonreciprocal matrix")


def test_four_by_four_weak_order_enumerator_has_4683_cases() -> None:
    assert sum(1 for _ in _ordered_weak_orders(6)) == 4683


def test_alpha_grid_has_101_values_for_every_n() -> None:
    values = alpha_grid(0.01)
    assert len(values) == 101
    assert values[0] == 1.0
    assert values[-1] == 0.0
    assert all(np.isclose(values[index] - values[index + 1], 0.01) for index in range(100))


def test_alpha_mnvdm_records_both_objective_components() -> None:
    sample = generate_pcm(3, "cyclic", 0, 20260815)
    settings = GurobiSettings(time_limit=30, threads=1, output_flag=0)
    result = solve_alpha_mnvdm(sample.matrix, 0.5, settings)
    assert result.has_solution and result.certified
    assert result.weights is not None and np.isclose(result.weights.sum(), 1.0)
    assert result.nvr is not None and result.gci_deviation is not None
    assert result.num_general_constraints == 0
    assert result.min_strict_slack is None or result.min_strict_slack >= -settings.feasibility_tol
    assert result.max_equality_residual is None or result.max_equality_residual <= settings.feasibility_tol
    expected = 0.5 * result.nvr + 0.5 * result.gci_deviation
    assert np.isclose(result.weighted_objective, expected, atol=1e-7)


def test_alpha_one_endpoint_has_fixed_state_certificate() -> None:
    sample = generate_pcm(3, "cyclic", 1, 20260815)
    settings = GurobiSettings(time_limit=30, threads=1, output_flag=0)
    interior = solve_alpha_mnvdm(sample.matrix, 0.99, settings)
    endpoint = solve_alpha_mnvdm(
        sample.matrix,
        1.0,
        settings,
        warm_start_y=interior.y,
        warm_start_states=interior.relation_states,
    )
    assert endpoint.has_solution and endpoint.certified
    assert endpoint.endpoint_fixed_state_status == "OPTIMAL"
    assert not endpoint.endpoint_fallback_used
    assert endpoint.min_strict_slack is None or endpoint.min_strict_slack >= -settings.feasibility_tol
    assert endpoint.max_equality_residual is None or endpoint.max_equality_residual <= settings.feasibility_tol
    assert np.isclose(endpoint.solver_objective, endpoint.nvr, atol=1e-7)


def test_check_first_evrim_stops_on_certified_gci_feasible_screening_solution() -> None:
    a = np.array([
        [1.0, 2.0, 4.0],
        [0.5, 1.0, 2.0],
        [0.25, 0.5, 1.0],
    ])
    settings = GurobiSettings(time_limit=30, threads=1, output_flag=0)
    result = solve_evrim_check_first(a, settings)
    assert result.solved and result.early_stop
    assert result.screening_solved
    assert result.bnc_runtime == 0.0
    assert result.gci is not None and result.gci <= result.gci_threshold + settings.gci_tolerance
