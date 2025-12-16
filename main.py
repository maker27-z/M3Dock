import argparse
from pathlib import Path

from lshade_nsga2.config import load_config
from lshade_nsga2.pipeline import run_docking


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Python reimplementation of LSHADE + NSGA2 based docking pipeline."
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to Vina-style configuration file (e.g. example/1bzc/1bzc_config.txt).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config_path = Path(args.config)
    if not config_path.is_file():
        raise SystemExit(f"Config file not found: {config_path}")

    cfg = load_config(config_path)
    run_docking(cfg)


if __name__ == "__main__":
    main()



