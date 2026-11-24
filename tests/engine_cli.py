"""End-to-end checks of the released executable and its CSV files."""

import csv
import math
from pathlib import Path
import subprocess
import sys
import tempfile


BINARY = Path(sys.argv[1]).resolve()
EXAMPLE = Path(sys.argv[2]).resolve()
CASE = sys.argv[3]


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def prepare_case(directory, engine_edits=(), simulation_edits=()):
    directory.mkdir()
    for kind, edits in (("engine", engine_edits), ("simulation", simulation_edits)):
        filename = f"{kind}_parameters.yaml"
        contents = (EXAMPLE / filename).read_text()
        for before, after in edits:
            require(before in contents, f"Example no longer contains {before!r}")
            contents = contents.replace(before, after)
        (directory / filename).write_text(contents)


def execute(directory):
    return subprocess.run(
        [str(BINARY)], cwd=directory, capture_output=True, text=True, timeout=90
    )


def read_trace(path, start, end):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        columns = reader.fieldnames
        require(columns and "Y_CH4[-]" in columns, f"Missing species columns in {path.name}")
        rows = []
        species = [name for name in columns if name.startswith("Y_")]
        for row in reader:
            require(None not in row and all(row.values()), f"Malformed data in {path.name}")
            values = {key: float(value) for key, value in row.items()}
            require(all(math.isfinite(value) for value in values.values()), f"Nonfinite value in {path.name}")
            for name in ("Pressure[Pa]", "Temperature[K]", "Volume[m3]", "Mass[kg]"):
                require(values[name] > 0.0, f"Nonpositive {name} in {path.name}")
            fractions = [values[key] for key in species]
            require(min(fractions) >= -1e-10 and abs(sum(fractions) - 1.0) < 3e-6,
                    f"Invalid species mass fractions in {path.name}")
            rows.append(values)
    angles = [row["CrankAngle[deg]"] for row in rows]
    require(len(angles) >= 2 and angles[0] == start and angles[-1] == end,
            f"{path.name} is missing the start or end state")
    require(all(left < right for left, right in zip(angles, angles[1:])),
            f"{path.name} angles are not increasing")
    return rows


def check_balances(metrics):
    intake = metrics["Intake Enthalpy"]
    outlet = metrics["Outlet Enthalpy"]
    heat = metrics["Wall Heat Into Gas"]
    work = metrics["Indicated Work"]
    stored = metrics["Cylinder Internal Energy Change"]
    residual = heat + intake - outlet - work - stored
    reported = metrics["Energy Balance Error"]
    require(abs(reported - residual) < 1e-8, "Printed energy terms do not reproduce the residual")
    require(abs(intake - metrics["Total Energy Accounted"] - residual) < 1e-8,
            "Total energy accounted does not match the first-law terms")
    require(abs(residual) < 5e-6, f"First-law balance failed with residual {residual * 1000:.6g} J")
    pv_residual = residual + work - metrics["Pressure-Volume Work"]
    require(abs(pv_residual - metrics["Pressure-Volume Balance Error"]) < 1e-8,
            "Pressure-volume residual is inconsistent with the independent work integral")
    require(abs(pv_residual) < 5e-5, f"Pressure-volume residual too large: {pv_residual * 1000:.6g} J")
    mass_error = metrics["Intake Mass"] - metrics["Outlet Mass"] - metrics["Cylinder Mass Change"]
    require(abs(mass_error) < 1e-12 and abs(metrics["Mass Balance Error"]) < 1e-12,
            f"Cycle mass balance failed with residual {mass_error:.6g} kg")
    require(abs(metrics["IMEP"] - metrics["GMEP"] - metrics["PMEP"]) < 1e-8,
            "Gross and pumping work do not sum to indicated work")
    require(abs(metrics["BMEP"] - metrics["IMEP"] + metrics["FMEP"]) < 1e-8,
            "Friction work is inconsistent with brake and indicated work")


def run_cycle(root, name, engine_edits=(), simulation_edits=(), start=0.0, end=720.0):
    directory = root / name
    prepare_case(directory, engine_edits, simulation_edits)
    result = execute(directory)
    require(result.returncode == 0, f"{name} exited {result.returncode}: {result.stderr}\n{result.stdout[-1500:]}")
    with (directory / "engine_performance.csv").open(newline="") as stream:
        metrics = {row["Parameter"]: float(row["Value"]) for row in csv.DictReader(stream)}
    require(all(math.isfinite(value) for value in metrics.values()), f"{name} produced nonfinite metrics")
    check_balances(metrics)
    motored = read_trace(directory / "motored_results.csv", start, end)
    fired = read_trace(directory / "engine_results.csv", start, end)
    require(abs(metrics["Cylinder Mass Change"] - (fired[-1]["Mass[kg]"] - fired[0]["Mass[kg]"])) < 2e-9,
            "Cycle mass inventory is inconsistent with the beginning and end states")
    print(f"{name}: IMEP={metrics['IMEP']:.6g} bar, "
          f"first-law residual={1000 * metrics['Energy Balance Error']:.6g} J, "
          f"independent p-V residual={1000 * metrics['Pressure-Volume Balance Error']:.6g} J")
    return metrics, motored, fired


def convergence(root):
    results = [
        run_cycle(root, f"step_{step}", simulation_edits=(
            ("crank_angle_step: 0.25", f"crank_angle_step: {step}"),
            ("adaptive_stepping: true", "adaptive_stepping: false"),
        ))[0]
        for step in (0.5, 0.25, 0.125)
    ]
    errors = [abs(metrics["Pressure-Volume Balance Error"]) for metrics in results]
    require(errors[1] < errors[0] and errors[2] < errors[1] and errors[2] < 5e-6,
            f"Independent pressure-volume work did not converge: {errors} kJ")
    for key, tolerance in (("IMEP", 0.02), ("T3 (Post-Combustion Temperature)", 2.0)):
        coarse, medium, fine = (metrics[key] for metrics in results)
        require(abs(medium - fine) < tolerance and abs(coarse - fine) < 2.0 * tolerance,
                f"{key} changed too much on refinement: {coarse}, {medium}, {fine}")


def operating_range(root):
    cases = (
        ("baseline", (), ()),
        ("low_speed_throttle", (("2000.0", "800.0"), ("intake_pressure: 101325.0", "intake_pressure: 75000.0")), ()),
        ("high_speed_backpressure", (("2000.0", "4500.0"), ("exhaust_pressure: 101325.0", "exhaust_pressure: 115000.0")), ()),
        ("frozen_expansion", (), (("equilibrium_expansion: true", "equilibrium_expansion: false"),)),
    )
    for name, engine_edits, simulation_edits in cases:
        metrics, motored, fired = run_cycle(root, name, engine_edits, simulation_edits)
        require(0.0 < metrics["Brake Power"] and 0.0 < metrics["Cycle Efficiency"] < 100.0,
                f"Unphysical fired-cycle performance in {name}")
        require(max(row["Temperature[K]"] for row in fired) > 2 * max(row["Temperature[K]"] for row in motored),
                f"Spark failed to heat the fired cycle in {name}")
        for trace in (motored, fired):
            require(all(200.0 < row["Temperature[K]"] < 4000.0 and
                        1e3 < row["Pressure[Pa]"] < 1e8 for row in trace),
                    f"Temperature or pressure escaped the test envelope in {name}")


def output_interval(root):
    dense, dense_motored, dense_fired = run_cycle(root, "dense_output")
    sparse, sparse_motored, sparse_fired = run_cycle(root, "sparse_output", simulation_edits=(
        ("output_interval: 0.5", "output_interval: 17.0"),
    ))
    require(len(sparse_fired) < len(dense_fired) and len(sparse_motored) < len(dense_motored),
            "Output interval did not change CSV sampling")
    require(dense == sparse, "Changing the CSV output interval changed the simulation")
    require(dense_fired[-1] == sparse_fired[-1], "Changing the output interval changed the final state")


def closed_cylinder(root):
    air_only = (("enabled: true", "enabled: false"), ("CH4:1, O2:2, N2:7.52", "O2:1, N2:3.76"))
    metrics, _, _ = run_cycle(root, "motored_air", engine_edits=air_only)
    require(metrics["Total LHV Energy Available"] == metrics["Indicated Efficiency"] == 0.0,
            "Spark-off air cycle reported fuel energy or efficiency")
    require(metrics["Wall Heat Into Gas"] > 0.0, "Net wall heating was reported with the wrong sign")

    metrics, _, fired = run_cycle(
        root, "closed_adiabatic",
        engine_edits=air_only + (("max_lift: 0.008", "max_lift: 0.0"),
                                 ("close_angle: 180.0", "close_angle: 0.0")),
        simulation_edits=(("woschni_c1: 3.26", "woschni_c1: 0.0"),
                          ("start_angle: 0.0", "start_angle: 65.0"),
                          ("end_angle: 720.0", "end_angle: 425.0")),
        start=65.0, end=425.0,
    )
    for key in ("Intake Mass", "Outlet Mass", "Wall Heat Into Gas"):
        require(metrics[key] == 0.0, f"Closed adiabatic cylinder produced a nonzero {key}")
    require(abs(metrics["Indicated Work"]) < 1e-7, "Reversible closed loop produced net work")
    for key in ("Pressure[Pa]", "Temperature[K]", "Volume[m3]", "Mass[kg]"):
        require(math.isclose(fired[0][key], fired[-1][key], rel_tol=5e-6),
                f"Closed loop did not return to its initial {key}")
    require(abs(fired[0]["Pressure[Pa]"] - 101325.0) < 0.1,
            "Nonzero starting angle changed the specified initial pressure")


def check_failure(directory, diagnostic):
    result = execute(directory)
    require(result.returncode == 1 and result.stderr.startswith("Error:") and diagnostic in result.stderr,
            f"{directory.name} should report {diagnostic!r}, got {result.returncode}: {result.stderr}")


def invalid_inputs(root):
    missing = root / "missing_config"
    missing.mkdir()
    check_failure(missing, "engine_parameters.yaml")
    malformed = root / "malformed_yaml"
    prepare_case(malformed, (("geometry:", "geometry: ["),))
    check_failure(malformed, "Failed to load engine parameters")
    cases = (
        ("zero_rpm", (("2000.0", "0.0"),), (), "engine_speed"),
        ("bad_geometry", (("compression_ratio: 10.0", "compression_ratio: 1.0"),), (), "compression_ratio"),
        ("unknown_species", (("CH4:1", "NOT_A_SPECIES:1"),), (), "NOT_A_SPECIES"),
        ("missing_mechanism", (("gri30.yaml", "missing_mechanism.yaml"),), (), "missing_mechanism.yaml"),
        ("zero_spark_duration", (("duration: 2.0", "duration: 0"),), (), "spark duration"),
        ("zero_valve_duration", (("close_angle: 180.0", "close_angle: 0.0"),), (), "valve close_angle"),
        ("zero_step", (), (("crank_angle_step: 0.25", "crank_angle_step: 0.0"),), "crank_angle_step"),
        ("nonfinite_step", (), (("crank_angle_step: 0.25", "crank_angle_step: .nan"),), "finite"),
        ("reversed_range", (), (("end_angle: 720.0", "end_angle: -1.0"),), "end_angle"),
        ("unresolvable_step", (), (("crank_angle_step: 0.25", "crank_angle_step: 1e-300"),
                                   ("adaptive_stepping: true", "adaptive_stepping: false"),
                                   ("start_angle: 0.0", "start_angle: 65.0")), "too small to advance"),
    )
    for name, engine_edits, simulation_edits, diagnostic in cases:
        directory = root / name
        prepare_case(directory, engine_edits, simulation_edits)
        check_failure(directory, diagnostic)
    blocked = root / "blocked_output"
    prepare_case(blocked, simulation_edits=(("end_angle: 720.0", "end_angle: 3.0"),))
    (blocked / "motored_results.csv").mkdir()
    check_failure(blocked, "Cannot open output file: motored_results.csv")


if __name__ == "__main__":
    cases = {
        "convergence": convergence,
        "operating_range": operating_range,
        "output_interval": output_interval,
        "closed_cylinder": closed_cylinder,
        "invalid_inputs": invalid_inputs,
    }
    with tempfile.TemporaryDirectory(prefix=f"dtpower_{CASE}_") as directory:
        cases[CASE](Path(directory))
