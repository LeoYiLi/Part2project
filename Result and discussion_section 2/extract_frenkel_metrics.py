#!/usr/bin/env python3
"""Extract Section 2 Frenkel-pair metrics from CIF, LOG, and TRAJ files."""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from ase.io import read


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
OUT_DIR = SCRIPT_DIR / "analysis_outputs"

SUPER_CELL = np.array([5, 5, 3], dtype=float)
PRISTINE_CIF = ROOT / "YBa2Cu3O7_553.cif"

OCCUPANCY_CUTOFF_A = 1.20
CU_O_CUTOFF_A = 2.60
SUSTAINED_TIME_PS = 0.05
DEFAULT_DT_PS = 0.0005


@dataclass(frozen=True)
class CaseInput:
    case_id: str
    cif: Path
    log: Path | None
    trajs: tuple[Path, ...]


CASE_INPUTS = (
    CaseInput(
        "F1",
        ROOT / "YBa2Cu3O7_553_Frenkel.cif",
        ROOT / "Frankel_md_298k.log",
        (ROOT / "Frankel_md_77K.traj",),
    ),
    CaseInput(
        "F2",
        ROOT / "YBa2Cu3O7_553_Frenkel2.cif",
        ROOT / "Frankel2_md_298k.log",
        (),
    ),
    CaseInput(
        "F3",
        ROOT / "YBa2Cu3O7_553_Frenkel3.cif",
        ROOT / "YBa2Cu3O7_553_Frenkel3.log",
        (
            ROOT / "YBa2Cu3O7_553_Frenkel3_113.traj",
            ROOT / "YBa2Cu3O7_553_Frenkel3_331.traj",
        ),
    ),
    CaseInput(
        "F4",
        ROOT / "YBa2Cu3O7_553_Frenkel4.cif",
        ROOT / "YBa2Cu3O7_553_Frenkel4.log",
        (),
    ),
)


def mic_delta_frac(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    delta = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    return delta - np.round(delta)


def mic_distance_frac(a: np.ndarray, b: np.ndarray, cell: np.ndarray) -> float:
    return float(np.linalg.norm(mic_delta_frac(a, b) @ cell))


def distances_to_site(site_frac: np.ndarray, frac_positions: np.ndarray, cell: np.ndarray) -> np.ndarray:
    delta = np.asarray(frac_positions, dtype=float) - np.asarray(site_frac, dtype=float)
    delta -= np.round(delta)
    return np.linalg.norm(delta @ cell, axis=1)


def pairwise_mic_distances(frac_a: np.ndarray, frac_b: np.ndarray, cell: np.ndarray) -> np.ndarray:
    delta = frac_a[:, None, :] - frac_b[None, :, :]
    delta -= np.round(delta)
    return np.linalg.norm(delta @ cell, axis=2)


def unit_cell_index(frac: np.ndarray) -> np.ndarray:
    scaled = np.asarray(frac) * SUPER_CELL
    return np.floor(scaled + 1e-8).astype(int)


def wrapped_cell_offset(cell_i: np.ndarray, cell_v: np.ndarray) -> np.ndarray:
    raw = np.asarray(cell_i, dtype=int) - np.asarray(cell_v, dtype=int)
    sizes = SUPER_CELL.astype(int)
    out = raw.copy()
    for i, size in enumerate(sizes):
        if out[i] > size / 2:
            out[i] -= size
        elif out[i] < -size / 2:
            out[i] += size
    return out


def fmt_vec(vec: np.ndarray, ndigits: int = 5) -> str:
    return "(" + ",".join(f"{float(x):.{ndigits}f}" for x in vec) + ")"


def fmt_int_vec(vec: np.ndarray) -> str:
    return "(" + ",".join(str(int(x)) for x in vec) + ")"


def primitive_site_type(frac: np.ndarray) -> str:
    prim = (np.asarray(frac) * SUPER_CELL) % 1.0
    z = float(prim[2])
    if min(abs(z), abs(z - 1.0)) < 0.06:
        return "chain"
    if min(abs(z - 0.15918), abs(z - 0.84082)) < 0.06:
        return "apical"
    if min(abs(z - 0.37831), abs(z - 0.62169)) < 0.06:
        return "plane"
    return "unknown"


def read_log_rows(path: Path | None) -> list[dict[str, float]]:
    if path is None or not path.exists():
        return []
    rows: list[dict[str, float]] = []
    for line in path.read_text(errors="ignore").splitlines():
        parts = line.split()
        if len(parts) < 5:
            continue
        try:
            time_ps, etot, epot, ekin, temp = map(float, parts[:5])
        except ValueError:
            continue
        rows.append(
            {
                "time_ps": time_ps,
                "Etot_eV": etot,
                "Epot_eV": epot,
                "Ekin_eV": ekin,
                "T_K": temp,
            }
        )
    return rows


def unique_log_times(rows: list[dict[str, float]]) -> np.ndarray:
    seen = set()
    times = []
    for row in rows:
        t = row["time_ps"]
        if t in seen:
            continue
        seen.add(t)
        times.append(t)
    return np.asarray(times, dtype=float)


def detect_defect_sites(pristine, defect) -> dict[str, object]:
    cell = pristine.cell.array
    pristine_symbols = np.asarray(pristine.get_chemical_symbols())
    defect_symbols = np.asarray(defect.get_chemical_symbols())
    pristine_o_idx = np.where(pristine_symbols == "O")[0]
    defect_o_idx = np.where(defect_symbols == "O")[0]
    pristine_o_frac = pristine.get_scaled_positions(wrap=True)[pristine_o_idx]
    defect_o_frac = defect.get_scaled_positions(wrap=True)[defect_o_idx]

    distances = pairwise_mic_distances(pristine_o_frac, defect_o_frac, cell)
    nearest_defect_distance = distances.min(axis=1)
    nearest_pristine_distance = distances.min(axis=0)

    vacancy_local = int(np.argmax(nearest_defect_distance))
    interstitial_local = int(np.argmax(nearest_pristine_distance))

    vacancy_atom_index = int(pristine_o_idx[vacancy_local])
    interstitial_atom_index = int(defect_o_idx[interstitial_local])
    vacancy_frac = pristine_o_frac[vacancy_local]
    interstitial_frac = defect_o_frac[interstitial_local]

    vacancy_cell = unit_cell_index(vacancy_frac)
    interstitial_cell = unit_cell_index(interstitial_frac)
    cell_offset = wrapped_cell_offset(interstitial_cell, vacancy_cell)

    pristine_frac = pristine.get_scaled_positions(wrap=True)
    pristine_cu_idx = np.where(pristine_symbols == "Cu")[0]
    cu_distances = distances_to_site(vacancy_frac, pristine_frac[pristine_cu_idx], cell)
    vacancy_cu_neighbors = int(np.sum(cu_distances <= CU_O_CUTOFF_A))

    return {
        "vacancy_atom_index_pristine": vacancy_atom_index,
        "interstitial_atom_index_defect": interstitial_atom_index,
        "vacancy_frac": vacancy_frac,
        "interstitial_frac": interstitial_frac,
        "vacancy_primitive_frac": (vacancy_frac * SUPER_CELL) % 1.0,
        "interstitial_primitive_frac": (interstitial_frac * SUPER_CELL) % 1.0,
        "vacancy_cell": vacancy_cell,
        "interstitial_cell": interstitial_cell,
        "cell_offset": cell_offset,
        "same_parent_unit_cell": bool(np.all(vacancy_cell == interstitial_cell)),
        "initial_oi_ov_distance_A": mic_distance_frac(interstitial_frac, vacancy_frac, cell),
        "vacancy_nearest_defect_o_distance_A": float(nearest_defect_distance[vacancy_local]),
        "interstitial_nearest_pristine_o_distance_A": float(nearest_pristine_distance[interstitial_local]),
        "vacancy_type": primitive_site_type(vacancy_frac),
        "vacancy_cu_neighbors_pristine": vacancy_cu_neighbors,
    }


def first_sustained_time(times: np.ndarray, flags: np.ndarray, sustained_ps: float) -> float | None:
    if len(flags) == 0:
        return None
    if len(times) > 1:
        dt = float(np.median(np.diff(times)))
        if dt <= 0:
            dt = DEFAULT_DT_PS
    else:
        dt = DEFAULT_DT_PS
    window = max(1, int(math.ceil(sustained_ps / dt)))
    if len(flags) < window:
        return None
    run = 0
    for idx, flag in enumerate(flags):
        run = run + 1 if bool(flag) else 0
        if run >= window:
            return float(times[idx - window + 1])
    return None


def infer_times(rows: list[dict[str, float]], n_frames: int) -> np.ndarray:
    log_times = unique_log_times(rows)
    if len(log_times) >= n_frames:
        return log_times[:n_frames]
    return np.arange(n_frames, dtype=float) * DEFAULT_DT_PS


def analyze_trajectory(case_id: str, case_info: dict[str, object], traj_path: Path, log_rows: list[dict[str, float]]) -> dict[str, object]:
    frames = read(traj_path, index=":")
    if not isinstance(frames, list):
        frames = [frames]
    n_frames = len(frames)
    times = infer_times(log_rows, n_frames)

    vacancy_frac = np.asarray(case_info["vacancy_frac"], dtype=float)
    interstitial_frac = np.asarray(case_info["interstitial_frac"], dtype=float)
    cell = frames[0].cell.array
    symbols = np.asarray(frames[0].get_chemical_symbols())
    o_idx = np.where(symbols == "O")[0]
    cu_idx = np.where(symbols == "Cu")[0]

    first_frac = frames[0].get_scaled_positions(wrap=True)
    first_o_frac = first_frac[o_idx]
    first_oi_site_distances = distances_to_site(interstitial_frac, first_o_frac, cell)
    tracked_oi_atom_index = int(o_idx[int(np.argmin(first_oi_site_distances))])
    tracked_oi_start_site_distance = float(first_oi_site_distances.min())

    ts_rows = []
    tracked_distances = []
    vacancy_site_distances = []
    coord_counts = []
    occupied_flags = []
    coord_recovered_flags = []

    coord_requirement = max(1, min(int(case_info["vacancy_cu_neighbors_pristine"]), 2))

    for frame_idx, atoms in enumerate(frames):
        frac = atoms.get_scaled_positions(wrap=True)
        o_frac = frac[o_idx]
        vac_site_dists = distances_to_site(vacancy_frac, o_frac, cell)
        nearest_o_local = int(np.argmin(vac_site_dists))
        nearest_o_atom_index = int(o_idx[nearest_o_local])
        nearest_o_frac = o_frac[nearest_o_local]
        nearest_o_vac_distance = float(vac_site_dists[nearest_o_local])

        tracked_distance = mic_distance_frac(frac[tracked_oi_atom_index], vacancy_frac, cell)

        if len(cu_idx):
            cu_distances_to_occupied_o = distances_to_site(nearest_o_frac, frac[cu_idx], cell)
            occupied_o_cu_coordination = int(np.sum(cu_distances_to_occupied_o <= CU_O_CUTOFF_A))
        else:
            occupied_o_cu_coordination = 0

        vacancy_occupied = nearest_o_vac_distance <= OCCUPANCY_CUTOFF_A
        coord_recovered = vacancy_occupied and occupied_o_cu_coordination >= coord_requirement

        tracked_distances.append(tracked_distance)
        vacancy_site_distances.append(nearest_o_vac_distance)
        coord_counts.append(occupied_o_cu_coordination)
        occupied_flags.append(vacancy_occupied)
        coord_recovered_flags.append(coord_recovered)

        ts_rows.append(
            {
                "case_id": case_id,
                "traj_file": traj_path.name,
                "frame": frame_idx,
                "time_ps": f"{times[frame_idx]:.6f}",
                "tracked_oi_atom_index": tracked_oi_atom_index,
                "tracked_oi_to_vacancy_distance_A": f"{tracked_distance:.6f}",
                "nearest_O_to_vacancy_site_atom_index": nearest_o_atom_index,
                "nearest_O_to_vacancy_site_distance_A": f"{nearest_o_vac_distance:.6f}",
                "vacancy_site_occupied": int(vacancy_occupied),
                "occupied_O_Cu_coordination": occupied_o_cu_coordination,
                "coordination_recovered": int(coord_recovered),
            }
        )

    tracked = np.asarray(tracked_distances, dtype=float)
    vac_site = np.asarray(vacancy_site_distances, dtype=float)
    coord = np.asarray(coord_counts, dtype=int)
    occupied = np.asarray(occupied_flags, dtype=bool)
    coord_recovered = np.asarray(coord_recovered_flags, dtype=bool)

    occupancy_time = first_sustained_time(times, occupied, SUSTAINED_TIME_PS)
    recombination_time = first_sustained_time(times, coord_recovered, SUSTAINED_TIME_PS)

    if n_frames == 1:
        recombined = "unclear"
        outcome = "single-frame only"
    elif recombination_time is not None:
        recombined = "yes"
        outcome = "recombined"
    elif occupancy_time is not None:
        recombined = "unclear"
        outcome = "vacancy occupied, coordination unclear"
    else:
        recombined = "no"
        outcome = "separated/metastable over sampled window"

    ts_path = OUT_DIR / f"{case_id}_{traj_path.stem}_timeseries.csv"
    write_csv(ts_path, ts_rows)

    return {
        "case_id": case_id,
        "traj_file": traj_path.name,
        "n_frames": n_frames,
        "n_atoms": len(frames[0]),
        "time_start_ps": f"{float(times[0]):.6f}" if len(times) else "",
        "time_end_ps": f"{float(times[-1]):.6f}" if len(times) else "",
        "tracked_oi_atom_index": tracked_oi_atom_index,
        "tracked_oi_start_site_distance_A": f"{tracked_oi_start_site_distance:.6f}",
        "tracked_oi_initial_to_vacancy_distance_A": f"{tracked[0]:.6f}",
        "tracked_oi_min_to_vacancy_distance_A": f"{tracked.min():.6f}",
        "tracked_oi_final_to_vacancy_distance_A": f"{tracked[-1]:.6f}",
        "vacancy_site_initial_nearest_O_distance_A": f"{vac_site[0]:.6f}",
        "vacancy_site_min_nearest_O_distance_A": f"{vac_site.min():.6f}",
        "vacancy_site_final_nearest_O_distance_A": f"{vac_site[-1]:.6f}",
        "max_occupied_O_Cu_coordination": int(coord.max()) if len(coord) else 0,
        "vacancy_occupied_fraction": f"{occupied.mean():.6f}",
        "coordination_recovered_fraction": f"{coord_recovered.mean():.6f}",
        "first_sustained_occupancy_time_ps": "" if occupancy_time is None else f"{occupancy_time:.6f}",
        "first_sustained_recombination_time_ps": "" if recombination_time is None else f"{recombination_time:.6f}",
        "recombined": recombined,
        "outcome_class": outcome,
        "timeseries_csv": ts_path.name,
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def summarize_log(case_id: str, log_path: Path | None, rows: list[dict[str, float]]) -> dict[str, object]:
    if not rows:
        return {
            "case_id": case_id,
            "log_file": "" if log_path is None else log_path.name,
            "status": "missing",
        }
    temps = np.asarray([r["T_K"] for r in rows], dtype=float)
    epot = np.asarray([r["Epot_eV"] for r in rows], dtype=float)
    times = np.asarray([r["time_ps"] for r in rows], dtype=float)
    unique_times = unique_log_times(rows)
    half = max(1, len(rows) // 2)
    return {
        "case_id": case_id,
        "log_file": "" if log_path is None else log_path.name,
        "status": "ok",
        "n_rows": len(rows),
        "n_unique_times": len(unique_times),
        "time_start_ps": f"{times[0]:.6f}",
        "time_end_ps": f"{times[-1]:.6f}",
        "T_mean_all_K": f"{temps.mean():.3f}",
        "T_mean_second_half_K": f"{temps[-half:].mean():.3f}",
        "T_final_K": f"{temps[-1]:.3f}",
        "Epot_initial_eV": f"{epot[0]:.6f}",
        "Epot_final_eV": f"{epot[-1]:.6f}",
    }


def main() -> None:
    OUT_DIR.mkdir(exist_ok=True)
    pristine = read(PRISTINE_CIF)

    inventory_rows = []
    log_rows_out = []
    event_rows = []

    for case in CASE_INPUTS:
        if not case.cif.exists():
            continue
        defect = read(case.cif)
        info = detect_defect_sites(pristine, defect)
        log_rows = read_log_rows(case.log)
        log_rows_out.append(summarize_log(case.case_id, case.log, log_rows))

        inventory_rows.append(
            {
                "case_id": case.case_id,
                "cif_file": case.cif.name,
                "vacancy_atom_index_pristine": info["vacancy_atom_index_pristine"],
                "vacancy_type": info["vacancy_type"],
                "vacancy_frac": fmt_vec(np.asarray(info["vacancy_frac"])),
                "vacancy_primitive_frac": fmt_vec(np.asarray(info["vacancy_primitive_frac"])),
                "interstitial_atom_index_defect": info["interstitial_atom_index_defect"],
                "interstitial_frac": fmt_vec(np.asarray(info["interstitial_frac"])),
                "interstitial_primitive_frac": fmt_vec(np.asarray(info["interstitial_primitive_frac"])),
                "vacancy_cell": fmt_int_vec(np.asarray(info["vacancy_cell"])),
                "interstitial_cell": fmt_int_vec(np.asarray(info["interstitial_cell"])),
                "same_parent_unit_cell": "yes" if info["same_parent_unit_cell"] else "no",
                "cell_offset": fmt_int_vec(np.asarray(info["cell_offset"])),
                "initial_oi_ov_distance_A": f"{float(info['initial_oi_ov_distance_A']):.6f}",
                "vacancy_cu_neighbors_pristine": info["vacancy_cu_neighbors_pristine"],
                "available_log": case.log.name if case.log and case.log.exists() else "",
                "available_traj": ";".join(p.name for p in case.trajs if p.exists()),
            }
        )

        for traj_path in case.trajs:
            if traj_path.exists():
                event_rows.append(analyze_trajectory(case.case_id, info, traj_path, log_rows))

    for case in CASE_INPUTS:
        if not case.trajs:
            event_rows.append(
                {
                    "case_id": case.case_id,
                    "traj_file": "",
                    "n_frames": 0,
                    "n_atoms": 0,
                    "time_start_ps": "",
                    "time_end_ps": "",
                    "tracked_oi_atom_index": "",
                    "tracked_oi_start_site_distance_A": "",
                    "tracked_oi_initial_to_vacancy_distance_A": "",
                    "tracked_oi_min_to_vacancy_distance_A": "",
                    "tracked_oi_final_to_vacancy_distance_A": "",
                    "vacancy_site_initial_nearest_O_distance_A": "",
                    "vacancy_site_min_nearest_O_distance_A": "",
                    "vacancy_site_final_nearest_O_distance_A": "",
                    "max_occupied_O_Cu_coordination": "",
                    "vacancy_occupied_fraction": "",
                    "coordination_recovered_fraction": "",
                    "first_sustained_occupancy_time_ps": "",
                    "first_sustained_recombination_time_ps": "",
                    "recombined": "not assessed",
                    "outcome_class": "missing trajectory",
                    "timeseries_csv": "",
                }
            )

    write_csv(OUT_DIR / "frenkel_static_case_inventory.csv", inventory_rows)
    write_csv(OUT_DIR / "frenkel_log_summary.csv", log_rows_out)
    write_csv(OUT_DIR / "frenkel_event_summary.csv", event_rows)
    print(f"Wrote outputs to {OUT_DIR}")


if __name__ == "__main__":
    main()
