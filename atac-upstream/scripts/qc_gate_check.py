#!/usr/bin/env python3
import argparse
import json
import os
import sys


COLORS = {
    "PASS": "\033[0;32m",
    "WARN": "\033[1;33m",
    "FAIL": "\033[0;31m",
    "RESET": "\033[0m",
}


DEFAULT_THRESHOLDS = {
    "tss_enrichment": {"direction": "higher", "pass": 10.0, "warn": 6.0, "label": "TSS enrichment"},
    "frip": {"direction": "higher", "pass": 0.30, "warn": 0.15, "label": "FRiP"},
    "nsc": {"direction": "higher", "pass": 1.05, "warn": 1.00, "label": "NSC"},
    "rsc": {"direction": "higher", "pass": 0.80, "warn": 0.50, "label": "RSC"},
    "nrf": {"direction": "higher", "pass": 0.70, "warn": 0.50, "label": "NRF"},
    "pbc": {"direction": "higher", "pass": 0.50, "warn": 0.30, "label": "PBC"},
    "mitochondrial_pct": {"direction": "lower", "pass": 20.0, "warn": 50.0, "label": "Mitochondrial %"},
    "idr_peak_count": {"direction": "higher", "pass": 70000.0, "warn": 50000.0, "label": "IDR peak count"},
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Evaluate ATAC-seq QC metrics against ENCODE-style thresholds."
    )
    parser.add_argument("--metrics-file", help="JSON file containing QC metrics")
    parser.add_argument("--thresholds-file", help="JSON file with custom thresholds")
    parser.add_argument("--json", action="store_true", help="Emit JSON output")
    parser.add_argument("--test", action="store_true", help="Run built-in test fixtures and exit")

    parser.add_argument("--tss-enrichment", type=float)
    parser.add_argument("--frip", type=float)
    parser.add_argument("--nsc", type=float)
    parser.add_argument("--rsc", type=float)
    parser.add_argument("--nrf", type=float)
    parser.add_argument("--pbc", type=float)
    parser.add_argument("--mitochondrial-pct", type=float)
    parser.add_argument("--idr-peak-count", type=float)
    return parser.parse_args()


def load_json_file(path):
    if not os.path.isfile(path):
        raise ValueError(f"File not found: {path}")
    with open(path, "r", encoding="utf-8") as handle:
        try:
            return json.load(handle)
        except json.JSONDecodeError as exc:
            raise ValueError(f"Invalid JSON in {path}: {exc}") from exc


def normalize_metrics(raw):
    if not isinstance(raw, dict):
        raise ValueError("Metrics input must be a JSON object")
    if "metrics" in raw and isinstance(raw["metrics"], dict):
        raw = raw["metrics"]

    metrics = {}
    for key in DEFAULT_THRESHOLDS:
        if key in raw and raw[key] is not None:
            try:
                metrics[key] = float(raw[key])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Metric '{key}' must be numeric") from exc
    return metrics


def load_thresholds(path):
    if not path:
        return DEFAULT_THRESHOLDS
    user = load_json_file(path)
    if not isinstance(user, dict):
        raise ValueError("Thresholds file must contain a JSON object")

    merged = {}
    for key, defaults in DEFAULT_THRESHOLDS.items():
        entry = user.get(key, {})
        if entry and not isinstance(entry, dict):
            raise ValueError(f"Threshold entry for '{key}' must be an object")
        merged[key] = {
            "direction": entry.get("direction", defaults["direction"]),
            "pass": float(entry.get("pass", defaults["pass"])),
            "warn": float(entry.get("warn", defaults["warn"])),
            "label": defaults["label"],
        }
        if merged[key]["direction"] not in {"higher", "lower"}:
            raise ValueError(f"Threshold direction for '{key}' must be 'higher' or 'lower'")
    return merged


def evaluate_metric(value, threshold):
    if threshold["direction"] == "higher":
        if value >= threshold["pass"]:
            return "PASS"
        if value >= threshold["warn"]:
            return "WARN"
        return "FAIL"

    if value <= threshold["pass"]:
        return "PASS"
    if value <= threshold["warn"]:
        return "WARN"
    return "FAIL"


def collect_metrics_from_args(args):
    return {
        "tss_enrichment": args.tss_enrichment,
        "frip": args.frip,
        "nsc": args.nsc,
        "rsc": args.rsc,
        "nrf": args.nrf,
        "pbc": args.pbc,
        "mitochondrial_pct": args.mitochondrial_pct,
        "idr_peak_count": args.idr_peak_count,
    }


def print_table(results):
    metric_width = max(len(item["label"]) for item in results)
    value_width = max(len(f"{item['value']:.4f}") for item in results)

    header = f"{'Metric'.ljust(metric_width)}  {'Value'.rjust(value_width)}  Status"
    print(header)
    print("-" * len(header))

    for item in results:
        status = item["status"]
        color = COLORS[status]
        reset = COLORS["RESET"]
        print(
            f"{item['label'].ljust(metric_width)}  "
            f"{item['value']:.4f}".rjust(value_width + 2)
            + f"  {color}{status}{reset}"
        )


def summary_to_exit(summary):
    if summary == "FAIL":
        return 2
    if summary == "WARN":
        return 1
    return 0


def summarize_status(results):
    statuses = {item["status"] for item in results}
    if "FAIL" in statuses:
        return "FAIL"
    if "WARN" in statuses:
        return "WARN"
    return "PASS"


def run_once(metrics, thresholds, output_json=False):
    missing = [key for key in DEFAULT_THRESHOLDS if key not in metrics]
    if missing:
        raise ValueError(
            "Missing required metrics: " + ", ".join(missing)
        )

    results = []
    for key, cfg in thresholds.items():
        value = float(metrics[key])
        status = evaluate_metric(value, cfg)
        results.append({
            "metric": key,
            "label": cfg["label"],
            "value": value,
            "status": status,
            "threshold": {
                "direction": cfg["direction"],
                "pass": cfg["pass"],
                "warn": cfg["warn"],
            },
        })

    summary = summarize_status(results)
    exit_code = summary_to_exit(summary)

    if output_json:
        payload = {
            "summary": summary,
            "exit_code": exit_code,
            "results": results,
        }
        print(json.dumps(payload, indent=2, sort_keys=False))
    else:
        print_table(results)
        color = COLORS[summary]
        reset = COLORS["RESET"]
        print(f"\nOverall: {color}{summary}{reset} (exit code {exit_code})")
    return exit_code


def run_tests(output_json=False):
    thresholds = DEFAULT_THRESHOLDS
    fixtures = [
        (
            "all_pass",
            {
                "tss_enrichment": 12,
                "frip": 0.35,
                "nsc": 1.1,
                "rsc": 0.9,
                "nrf": 0.8,
                "pbc": 0.6,
                "mitochondrial_pct": 15,
                "idr_peak_count": 75000,
            },
            0,
        ),
        (
            "warn_case",
            {
                "tss_enrichment": 8,
                "frip": 0.16,
                "nsc": 1.02,
                "rsc": 0.7,
                "nrf": 0.6,
                "pbc": 0.4,
                "mitochondrial_pct": 35,
                "idr_peak_count": 60000,
            },
            1,
        ),
        (
            "fail_case",
            {
                "tss_enrichment": 4,
                "frip": 0.10,
                "nsc": 0.95,
                "rsc": 0.4,
                "nrf": 0.4,
                "pbc": 0.2,
                "mitochondrial_pct": 60,
                "idr_peak_count": 40000,
            },
            2,
        ),
    ]

    failures = 0
    for name, metrics, expected in fixtures:
        print(f"[TEST] {name}")
        code = run_once(metrics, thresholds, output_json=output_json)
        if code != expected:
            failures += 1
            print(f"Expected exit code {expected}, got {code}", file=sys.stderr)
        print()

    if failures:
        print(f"{failures} test fixture(s) failed.", file=sys.stderr)
        return 1
    print("All test fixtures passed.")
    return 0


def main():
    args = parse_args()

    if args.test:
        return run_tests(output_json=args.json)

    try:
        thresholds = load_thresholds(args.thresholds_file)

        metrics = {}
        if args.metrics_file:
            metrics.update(normalize_metrics(load_json_file(args.metrics_file)))

        for key, value in collect_metrics_from_args(args).items():
            if value is not None:
                metrics[key] = float(value)

        if not metrics:
            raise ValueError(
                "No metrics provided. Use --metrics-file or individual --key=value arguments."
            )

        return run_once(metrics, thresholds, output_json=args.json)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
