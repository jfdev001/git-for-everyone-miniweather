#!/usr/bin/env python

from collections import namedtuple
from glob import glob
from os import listdir
from os.path import join, realpath, basename, getmtime
from pathlib import Path
from re import findall, search

import matplotlib.pyplot as plt

FORTRAN_PROJECT: Path = Path(realpath(__file__)).parents[2]
BUILD_OUTPUT_DIR: str = join(
    FORTRAN_PROJECT, "build", "build_output")


def main():
    # get the subdirectories corresponding to the experiments
    subdirs: list[str] = listdir(BUILD_OUTPUT_DIR)
    experiment_subdirs: list[str] = [
        join(BUILD_OUTPUT_DIR, dir) for dir in subdirs if "nx_" in dir]

    # extract parameters based on these subdirectories
    # NOTE: you could adapt this based on your naming conventions
    SimParams = namedtuple(
        "SimParams", (
            "nx",
            "nz",
            "sim_time",
            "out_freq"
        )
    )
    sim_params_set: list[SimParams] = []
    subdir_name_regex = r"nx_(\d*)_nz_(\d*)_sim_time_(\d*)_out_freq_(\d*)"
    for subdir in experiment_subdirs:
        match_sim_params: list[tuple[str]] = findall(
            subdir_name_regex, subdir)
        assert len(match_sim_params) == 1, "must find only one match"
        sim_params: SimParams = SimParams(
            *tuple(map(lambda p: int(p), match_sim_params[0])))
        sim_params_set.append(sim_params)

    # extract measurements from output directories and build object holding
    # experiment simulation parameters and measurements
    # NOTE: You could adapt this to your measurements
    Measurement = namedtuple(
        "Measurement", (
            "nthreads",
            "cpu_time"
        )
    )

    Experiment = namedtuple(
        "Experiment", (
            "sim_params",
            "results"
        )
    )

    experiments: list[Experiment] = []
    experimental_output_dir_regex = r"output_nthread_(\d*)"
    for subdir, sim_params in zip(experiment_subdirs, sim_params_set):

        experimental_output_dirs = [
            join(subdir, d) for d in listdir(subdir) if "output_" in d]
        results: list[Measurement] = []
        for experimental_output_dir in experimental_output_dirs:
            experimental_output_dir_basename = basename(
                experimental_output_dir)

            match_nthreads = findall(
                experimental_output_dir_regex,
                experimental_output_dir_basename)
            assert len(match_nthreads) == 1, "must find only one match"
            nthreads = int(match_nthreads[0])

            log_files = get_mtime_sorted_log_files(experimental_output_dir)
            most_recent_log_file = log_files[-1]

            cpu_time = get_cpu_time(most_recent_log_file)

            measurement = Measurement(nthreads, cpu_time)
            results.append(measurement)

        experiment = Experiment(sim_params, results)
        experiments.append(experiment)

    # -- plot performance --
    # NOTE: you could adapt this to plot efficiency or other things
    fig, axs = plt.subplots(figsize=(10, 8), nrows=2)
    fig.suptitle("Miniweather OpenMP Performance Results")

    # ----- Plot CPU Time scaling -----
    for exp in experiments:
        # Sort results by nthreads for consistent plotting
        sorted_results = sorted(exp.results, key=lambda m: m.nthreads)
        nthreads = [m.nthreads for m in sorted_results]
        cpu_time = [m.cpu_time for m in sorted_results]

        label = f"nx={exp.sim_params.nx}, nz={exp.sim_params.nz}"
        axs[0].plot(nthreads, cpu_time, marker='o', label=label)

    axs[0].set_ylabel("CPU Time [s]")
    axs[0].set_title("CPU Time Scaling")
    axs[0].grid(True, linestyle='--', alpha=0.5)
    axs[0].legend()
    axs[0].set_xticks([1, 2, 4])  # adjust based on your thread counts

    # --- plot speedup ---
    for exp in experiments:
        sorted_results = sorted(exp.results, key=lambda m: m.nthreads)
        nthreads = [m.nthreads for m in sorted_results]
        cpu_time = [m.cpu_time for m in sorted_results]
        # speedup = single-thread time / parallel time
        speedup = [cpu_time[0]/t for t in cpu_time]

        label = f"nx={exp.sim_params.nx}, nz={exp.sim_params.nz}"
        axs[1].plot(nthreads, speedup, marker='o', label=label)

    # ideal speedup
    ideal_threads = sorted(
        set([m.nthreads for exp in experiments for m in exp.results]))
    ideal_speedup = ideal_threads  # linear scaling: speedup = nthreads
    axs[1].plot(ideal_threads, ideal_speedup,
                'k--', label="Ideal linear scaling")

    axs[1].set_xlabel("Number of Threads")
    axs[1].set_ylabel("Speedup")
    axs[1].set_title("Parallel Speedup")
    axs[1].grid(True, linestyle='--', alpha=0.5)
    axs[1].legend()
    axs[1].set_xticks([1, 2, 4])

    plt.show()
    return


def get_mtime_sorted_log_files(experimental_output_dir: str) -> list[str]:
    log_files: list[str] = glob(
        join(experimental_output_dir, "**", "*.out"), recursive=True)
    sorted_log_files: list[str] = sorted(log_files, key=lambda f: getmtime(f))
    return sorted_log_files


def get_cpu_time(log_file: str) -> float:
    cpu_time: float
    cpu_time_regex = r"CPU Time:\s*([0-9.]+)"
    with open(log_file, "r") as f:
        for line in f:
            match = search(cpu_time_regex, line)
            if match:
                cpu_time = float(match.group(1))
                break  # stop after first match
    return cpu_time


if __name__ == "__main__":
    main()
