"""Compare compiled logdiffexp kernels with mpmath on selected inputs"""

import argparse
import csv
import ctypes
import hashlib
import json
import math
import pathlib
import platform
import shlex
import shutil
import statistics
import subprocess
import sys
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np


SOURCE = r'''
#include <xsf/log_exp.h>
#include <chrono>
#include <cstddef>

extern "C" double original(double a, double b) {
    return a + xsf::log1mexp(b - a);
}

extern "C" double current(double a, double b) {
    return xsf::logdiffexp(a, b);
}

extern "C" double benchmark(const double* a, const double* b,
                            std::size_t count, std::size_t rounds,
                            int corrected) {
    volatile double sink = 0;
    const auto start = std::chrono::steady_clock::now();
    if (corrected) {
        for (std::size_t r = 0; r < rounds; ++r) {
            for (std::size_t i = 0; i < count; ++i) {
                sink = xsf::logdiffexp(a[i], b[i]);
            }
        }
    } else {
        for (std::size_t r = 0; r < rounds; ++r) {
            for (std::size_t i = 0; i < count; ++i) {
                sink = a[i] + xsf::log1mexp(b[i] - a[i]);
            }
        }
    }
    const auto elapsed = std::chrono::steady_clock::now() - start;
    return std::chrono::duration<double, std::nano>(elapsed).count()
           / (count * rounds);
}
'''


def sample_inputs():
    # generate the zero-result curve and its neighboring doubles
    for b in np.linspace(-100.0, 30.0, 201):
        b = float(b)
        center = float(mp.log1p(mp.exp(mp.mpf(b))))
        for a in (math.nextafter(center, -math.inf), center,
                  math.nextafter(center, math.inf)):
            yield "near zero", b, a, b

    # sample offsets and separations without claiming full domain coverage
    for a in (-1000.0, -10.0, -1.0, 0.0, 1.0, 10.0, 1000.0):
        for gap in np.geomspace(1e-12, 700.0, 90):
            b = a - float(gap)
            if a > b:
                yield "general", a - b, a, b

    # keep subnormal results separate from ordinary relative-error summaries
    for b in (-745.0, -740.0, -730.0, -710.0, -701.0, -700.0, -699.0):
        center = float(mp.log1p(mp.exp(mp.mpf(b))))
        for a in (math.nextafter(center, -math.inf), center,
                  math.nextafter(center, math.inf)):
            yield "extreme tail", b, a, b


def measure_accuracy(loaded):
    functions = {}
    for name in ("original", "current"):
        function = getattr(loaded, name)
        function.argtypes = [ctypes.c_double, ctypes.c_double]
        function.restype = ctypes.c_double
        functions[name] = function

    rows = []
    for region, coordinate, a, b in sample_inputs():
        # pass exact binary inputs to mpmath and avoid reference cancellation
        aa, bb = mp.mpf(a), mp.mpf(b)
        difference = bb - aa
        correction = (mp.log1p(-mp.exp(difference)) if difference < -1
                      else mp.log(-mp.expm1(difference)))
        reference = aa + correction
        spacing = mp.mpf(math.ulp(float(reference)))
        for name, function in functions.items():
            observed = function(a, b)
            absolute = abs(mp.mpf(observed) - reference)
            relative = absolute / abs(reference) if reference != 0 else None
            rows.append({
                "region": region, "method": name, "coordinate": coordinate,
                "a_hex": a.hex(), "b_hex": b.hex(),
                "reference": mp.nstr(reference, 60),
                "observed": repr(observed),
                "absolute_error": mp.nstr(absolute, 25),
                "relative_error": mp.nstr(relative, 25) if relative is not None else "",
                "spacing_error": mp.nstr(absolute / spacing, 25),
                "rounded_reference_match": observed == float(reference),
            })
    return rows


def measure_timings(loaded, rounds, repeats):
    # time compiled loops so Python call overhead is outside each measurement
    groups = {
        "ordinary": [(1000.0, 999.0), (1000.0, 999.999999999),
                     (0.0, -1e-12), (0.0, -50.0)],
        "near zero": [(a, 0.0) for a in
                      (0.6931471805599452, 0.6931471805599453, 0.6931471805599454)],
        "extreme tail": [(float.fromhex(a), -701.0) for a in
                         ("0x1.9788af0d610e5p-1012", "0x1.9788af0d610e6p-1012",
                          "0x1.9788af0d610e7p-1012")],
    }
    bench = loaded.benchmark
    pointer = ctypes.POINTER(ctypes.c_double)
    bench.argtypes = [pointer, pointer, ctypes.c_size_t, ctypes.c_size_t, ctypes.c_int]
    bench.restype = ctypes.c_double
    results = {}
    print("Median nanoseconds per scalar evaluation")
    print(f"{'region':<16} {'original':>12} {'current':>12} {'ratio':>10}")
    for region, pairs in groups.items():
        array = ctypes.c_double * len(pairs)
        a, b = array(*(p[0] for p in pairs)), array(*(p[1] for p in pairs))
        timings = {}
        for method, name in enumerate(("original", "current")):
            bench(a, b, len(pairs), 1000, method)
            samples = [bench(a, b, len(pairs), rounds, method) for _ in range(repeats)]
            timings[name] = {"samples_ns": samples, "median_ns": statistics.median(samples)}
        original = timings["original"]["median_ns"]
        current = timings["current"]["median_ns"]
        timings["ratio"] = current / original
        results[region] = timings
        print(f"{region:<16} {original:12.1f} {current:12.1f} {current / original:9.2f}x")
    return results


def make_plot(rows, destination):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for row_index, region in enumerate(("near zero", "general")):
        for column, metric in enumerate(("relative_error", "spacing_error")):
            ax = axes[row_index, column]
            for method, color in (("original", "tab:orange"), ("current", "tab:blue")):
                selected = [r for r in rows if r["region"] == region
                            and r["method"] == method and r[metric] != ""]
                ax.scatter([r["coordinate"] for r in selected],
                           [float(r[metric]) for r in selected],
                           s=10, alpha=0.55, label=method, color=color)
            ax.set_yscale("symlog", linthresh=1e-18 if column == 0 else 0.1)
            ax.yaxis.get_major_locator().set_params(numticks=9)
            if region == "general":
                ax.set_xscale("log")
                ax.set_xlabel("Argument separation a − b")
            else:
                ax.set_xlabel("b, with a near log(1 + exp(b))")
            ax.set_ylabel("Relative error" if column == 0
                          else "Absolute error / output spacing")
            ax.set_title(f"{region.capitalize()} — sampled inputs")
            ax.grid(True, which="both", alpha=0.2)
            ax.legend()
    fig.suptitle(f"Compiled logdiffexp versus {mp.mp.dps}-digit mpmath references")
    fig.tight_layout(rect=[0, 0.06, 1, 0.96])
    fig.text(0.5, 0.015,
             "Output spacing is measured at the rounded reference\n"
             "Error axes are linear below 10⁻¹⁸ for relative error "
             "and below 0.1 for spacing error",
             ha="center", fontsize=9)
    fig.savefig(destination, dpi=200)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=pathlib.Path, default=pathlib.Path.cwd())
    parser.add_argument("--output", type=pathlib.Path,
                        default=pathlib.Path("../logdiffexp-analysis/reproducible"))
    parser.add_argument("--cxx", default=shutil.which("clang++") or "c++")
    parser.add_argument("--rounds", type=int, default=100000)
    parser.add_argument("--repeats", type=int, default=5)
    args = parser.parse_args()
    if args.rounds <= 0 or args.repeats <= 0:
        parser.error("rounds and repeats must be positive")
    if sys.platform not in ("darwin", "linux"):
        parser.error("the shared-library build currently supports macOS and Linux")
    repo, output = args.repo.resolve(), args.output.resolve()
    if not (repo / "include/xsf/log_exp.h").is_file():
        parser.error("--repo must identify an XSF checkout")
    output.mkdir(parents=True, exist_ok=True)
    mp.mp.dps = 500
    compiler = shlex.split(args.cxx)
    flags = ["-std=c++17", "-O3"]
    flags += ["-dynamiclib"] if sys.platform == "darwin" else ["-shared", "-fPIC"]
    with tempfile.TemporaryDirectory(prefix="logdiffexp-analysis-") as directory:
        directory = pathlib.Path(directory)
        source, library = directory / "probe.cpp", directory / "probe.so"
        source.write_text(SOURCE)
        subprocess.run(compiler + flags + ["-I", str(repo / "include"),
                       str(source), "-o", str(library)], check=True)
        loaded = ctypes.CDLL(str(library))
        timings = measure_timings(loaded, args.rounds, args.repeats)
        rows = measure_accuracy(loaded)
    with (output / "accuracy.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    make_plot(rows, output / "accuracy.png")
    git = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"],
                         capture_output=True, text=True, check=True)
    status = subprocess.run(["git", "-C", str(repo), "status", "--short"],
                            capture_output=True, text=True, check=True)
    version = subprocess.run(compiler + ["--version"], capture_output=True,
                             text=True, check=True)
    metadata = {
        "commit": git.stdout.strip(), "working_tree_status": status.stdout,
        "header_sha256": hashlib.sha256((repo / "include/xsf/log_exp.h").read_bytes()).hexdigest(),
        "script_sha256": hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),
        "platform": platform.platform(), "machine": platform.machine(),
        "python": sys.version, "compiler": version.stdout.splitlines()[0],
        "compiler_flags": flags, "mpmath": mp.__version__, "dps": mp.mp.dps,
        "numpy": np.__version__, "matplotlib": matplotlib.__version__,
        "rounds": args.rounds, "repeats": args.repeats, "timings": timings,
        "input_pairs": len(rows) // 2,
        "limitations": "selected double inputs only; timings include loop and volatile-store overhead",
    }
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"\nCompared {len(rows) // 2} input pairs")
    for method in ("original", "current"):
        for region in ("near zero", "general", "extreme tail"):
            selected = [r for r in rows if r["region"] == region and r["method"] == method]
            worst = max(mp.mpf(r["relative_error"]) for r in selected if r["relative_error"])
            matches = sum(r["rounded_reference_match"] for r in selected)
            print(f"{method}, {region}: max relative error {mp.nstr(worst, 6)}, "
                  f"{matches}/{len(selected)} rounded matches")
    print(f"\nOutputs: {output}")


if __name__ == "__main__":
    main()
