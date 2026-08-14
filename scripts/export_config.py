#!/usr/bin/env python3
"""Port of the upstream genome_tracks 'config_export' rule run block.

Upstream json.dump()s the loaded Snakemake config dict into
{result_path}/configs/{project_name}_config.yaml. The oxo-flow equivalent of
"the loaded config" is the [config] table of the workflow file.

Usage:
  export_config.py WORKFLOW.oxoflow OUT.yaml
"""
import json
import sys
import tomllib


def main():
    if len(sys.argv) != 3:
        print("usage: export_config.py WORKFLOW.oxoflow OUT.yaml", file=sys.stderr)
        return 1
    with open(sys.argv[1], "rb") as f:
        workflow = tomllib.load(f)
    config = workflow.get("config", {})
    with open(sys.argv[2], "w") as outfile:
        json.dump(config, outfile)
    return 0


if __name__ == "__main__":
    sys.exit(main())
