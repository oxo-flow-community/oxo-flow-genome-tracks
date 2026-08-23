#!/usr/bin/env python3
"""Port of the upstream genome_tracks 'ucsc_hub' rule run block.

Creates a UCSC genome browser track hub for all bigWig files:
  * relative symlink per group inside bigWigs/{genome}/
  * genomes.txt, hub.txt, trackDb.txt with per-group colors (hex -> RGB).

Usage:
  ucsc_hub.py --bw-dir DIR --genome LABEL --project-name NAME --email MAIL \
      --groups g1,g2 --colors "g1=#hex,g2=#hex"
"""
import argparse
import os
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bw-dir", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--project-name", required=True)
    ap.add_argument("--email", required=True)
    ap.add_argument("--groups", required=True)
    ap.add_argument("--colors", required=True)
    args = ap.parse_args()

    groups = [g for g in args.groups.split(",") if g]

    # config["track_colors"][group] if group in config["track_colors"] else "#000000"
    color_dict = {}
    for entry in [e for e in args.colors.split(",") if e]:
        g, _, c = entry.partition("=")
        color_dict[g] = c

    genome_dir = os.path.join(args.bw_dir, args.genome)
    os.makedirs(genome_dir, exist_ok=True)

    # create bigWig symlinks (as upstream: os.symlink('../' + basename, ...));
    # upstream's plain symlink is not idempotent — re-runs die with
    # FileExistsError (live on resume) — unlink first.
    for group in groups:
        dst = os.path.join(genome_dir, "{}.bw".format(group))
        if os.path.lexists(dst):
            os.unlink(dst)
        os.symlink(os.path.join("..", "{}.bw".format(group)), dst)

    # create genomes.txt
    with open(os.path.join(args.bw_dir, "genomes.txt"), "w") as gf:
        gf.write(
            "genome {}\ntrackDb {}/trackDb.txt\n".format(args.genome, args.genome)
        )

    # create hub file
    with open(os.path.join(args.bw_dir, "hub.txt"), "w") as hf:
        hub_text = [
            "hub {}".format(args.project_name),
            "shortLabel {}".format(args.project_name),
            "longLabel {}".format(args.project_name),
            "genomesFile genomes.txt",
            "email {}\n".format(args.email),
        ]
        hf.write("\n".join(hub_text))

    # create trackdb file
    with open(os.path.join(genome_dir, "trackDb.txt"), "w") as tf:
        track_db = [
            "track {}".format(args.project_name),
            "type bigWig",
            "compositeTrack on",
            "autoScale on",
            "maxHeightPixels 32:32:8",
            "shortLabel {}".format(args.project_name[:8]),
            "longLabel {}".format(args.project_name),
            "visibility full",
            "",
            "",
        ]
        for group in groups:
            hex_color = color_dict.get(group, "#000000")
            # convert to RGB: tuple(int(hex[i:i+2], 16) for i in (1, 3, 5))
            track_color = ",".join(
                str(int(hex_color[i : i + 2], 16)) for i in (1, 3, 5)
            )
            track = [
                "track {}".format(group),
                "shortLabel {}".format(group),
                "longLabel {}".format(group),
                "bigDataUrl {}.bw".format(group),
                "parent {} on".format(args.project_name),
                "type bigWig",
                "windowingFunction mean",
                "color {}".format(track_color),
                "",
                "",
            ]
            track_db += track
        tf.write("\n".join(track_db))


if __name__ == "__main__":
    sys.exit(main())
