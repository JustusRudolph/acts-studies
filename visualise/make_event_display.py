#!/usr/bin/env python3
"""Build an interactive HTML event display from an ACTSO2 output directory.

Reads (all optional except hits.root and particles_simulation.root):
    hits.root                                  simulated hits (sensitive layers)
    particles_simulation.root                  Geant4 particles: primaries + kept secondaries
    trackstates_ambi-tracks-merged.root        reconstructed tracks (measurements + smoothed states)
    tracksummary_ambi-tracks-merged.root       majority particle and classification of each track
    performance_merged_ambi_tracks.root        per-particle reconstruction matching (matchingdetails)
    iter_N/*performance_seeding.root           per-particle seed matching of each iteration

and writes a single self-contained HTML file with x-y, r-z and 3D views, layer toggles,
and an event > primaries/secondaries > particle browser with filters on the particles'
number of hits, eta, pT and seeding/reconstruction status.

    alienv setenv ACTSO2/latest -c python3 make_event_display.py \
        /home/justus/projects/alice/ACTSO2/output/test_loopers_10000ev \
        --geo-config /home/justus/projects/alice/ACTSO2/geometries/geometry_2026_10_02_segmented_realistic/geom/gen3_geometry_config.json

The output goes to figures/<sample>_event_display.html (gitignored) unless -o is given.
"""

import argparse
import json
import math
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)

BARCODE_FIELDS = ("vertex_primary", "vertex_secondary", "particle", "generation", "sub_particle")
INVALID_POSITION = 1e5  # trackstates use large sentinels for missing positions
MM_PER_NS = 299.792458  # the particle (vt) and hit (tt) writers store time as a length in mm (c = 1)
# events with more simulated particles than this (PbPb) only get the particles with hits in the decay tree
DECAY_TREE_MAX_PARTICLES = 50000
# no-hit secondaries kept per event for the production-vertex markers
MAX_VERTICES_PER_EVENT = 5000
# a track measurement belongs to a simulated hit of its particle if it is on the same layer and
# this close; measured positions are within ~2 mm of the hit, others come from other particles
HIT_MATCH_DISTANCE = 5.0  # mm


def progress(iterable, total, label):
    """Yield from iterable while drawing a progress bar on stderr (plain lines when not a terminal)."""
    is_tty = sys.stderr.isatty()
    width = 30
    start = time.monotonic()
    last_draw = 0.0
    last_percent = 0

    def render(count):
        fraction = count / total if total else 1.0
        elapsed = time.monotonic() - start
        eta = elapsed / fraction - elapsed if 0 < fraction < 1 else 0.0
        filled = int(width * fraction)
        bar = "█" * filled + "·" * (width - filled)
        return f"{label:<12} |{bar}| {100 * fraction:5.1f}%  {count}/{total}  {elapsed:5.1f}s  eta {eta:5.1f}s"

    count = 0
    for item in iterable:
        yield item
        count += 1
        now = time.monotonic()
        if is_tty and now - last_draw > 0.1:
            sys.stderr.write("\r" + render(count))
            sys.stderr.flush()
            last_draw = now
        elif not is_tty and total and count < total and (percent := 100 * count // total) >= last_percent + 25:
            sys.stderr.write(render(count) + "\n")
            last_percent = percent - percent % 25
    sys.stderr.write(("\r" if is_tty else "") + render(count) + "\n")
    sys.stderr.flush()


def round_list(values, ndigits=2):
    return [round(float(value), ndigits) for value in values]


def is_valid(*values):
    return all(math.isfinite(value) and abs(value) < INVALID_POSITION for value in values)


def barcode_key(barcode):
    return "-".join(str(part) for part in barcode)


def family_key(barcode):
    """Geant4 descendants keep (vertex_primary, vertex_secondary, particle) of their primary."""
    return "-".join(str(part) for part in barcode[:3])


def read_particles(path):
    """Per event: the particle columns as numpy arrays (PbPb events have millions of particles)."""
    columns = ["event_id", "particle_type", "q", "m", "p", "pt", "eta", "phi", "px", "py", "pz",
               "vx", "vy", "vz", "vt", "e_loss", "process", *BARCODE_FIELDS]
    data = ROOT.RDataFrame("particles", str(path)).AsNumpy(columns)
    particles_by_event = {}
    for i_event in progress(range(len(data["event_id"])), len(data["event_id"]), "particles"):
        arrays = {column: np.asarray(data[column][i_event]) for column in columns[1:]}
        arrays["vt"] = arrays["vt"] / MM_PER_NS
        arrays["tkin"] = np.sqrt(arrays["p"] ** 2 + arrays["m"] ** 2) - arrays["m"]
        arrays["barcode"] = list(zip(*(arrays[field].tolist() for field in BARCODE_FIELDS)))
        particles_by_event[int(data["event_id"][i_event])] = arrays
    return particles_by_event


def particle_dict(arrays, i_part):
    return {
        "barcode": arrays["barcode"][i_part],
        "pdg": int(arrays["particle_type"][i_part]),
        "q": float(arrays["q"][i_part]),
        "vx": float(arrays["vx"][i_part]), "vy": float(arrays["vy"][i_part]), "vz": float(arrays["vz"][i_part]),
        "vt": float(arrays["vt"][i_part]),
        "px": float(arrays["px"][i_part]), "py": float(arrays["py"][i_part]), "pz": float(arrays["pz"][i_part]),
        "pt": float(arrays["pt"][i_part]), "eta": float(arrays["eta"][i_part]), "phi": float(arrays["phi"][i_part]),
        "p": float(arrays["p"][i_part]),
        "tkin": float(arrays["tkin"][i_part]),
        "eloss": float(arrays["e_loss"][i_part]),
        "process": int(arrays["process"][i_part]),
    }


def read_hits(path):
    """Per event and barcode: hit arrays, plus per-layer position samples."""
    columns = ["event_id", "tx", "ty", "tz", "tt", "tpx", "tpy", "tpz", "volume_id", "layer_id",
               *[f"barcode_{field}" for field in BARCODE_FIELDS]]
    data = ROOT.RDataFrame("hits", str(path)).AsNumpy(columns)
    data = {column: values.tolist() for column, values in data.items()}
    hits = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    layer_samples = defaultdict(lambda: {"r": [], "z": []})
    for i_hit in progress(range(len(data["event_id"])), len(data["event_id"]), "hits"):
        event_id = data["event_id"][i_hit]
        barcode = tuple(data[f"barcode_{field}"][i_hit] for field in BARCODE_FIELDS)
        x, y, z = data["tx"][i_hit], data["ty"][i_hit], data["tz"][i_hit]
        layer_key = (data["volume_id"][i_hit], data["layer_id"][i_hit])
        p_hit = math.sqrt(data["tpx"][i_hit] ** 2 + data["tpy"][i_hit] ** 2 + data["tpz"][i_hit] ** 2)
        particle_hits = hits[event_id][barcode]
        particle_hits["x"].append(x)
        particle_hits["y"].append(y)
        particle_hits["z"].append(z)
        particle_hits["t"].append(data["tt"][i_hit] / MM_PER_NS)
        particle_hits["p"].append(p_hit)
        particle_hits["layer"].append(f"{layer_key[0]}/{layer_key[1]}")
        samples = layer_samples[layer_key]
        if len(samples["r"]) < 5000:
            samples["r"].append(math.hypot(x, y))
            samples["z"].append(z)
    # hits.root is ordered by detector id; number each particle's hits in time order instead
    for event_hits in hits.values():
        for particle_hits in event_hits.values():
            time_order = sorted(range(len(particle_hits["t"])), key=lambda i_hit: particle_hits["t"][i_hit])
            for column, values in particle_hits.items():
                particle_hits[column] = [values[i_hit] for i_hit in time_order]
    return hits, layer_samples


def build_sensitive_layers(layer_samples):
    """Cylinders (barrel) or discs (endcap), estimated from where the hits are."""
    layers = []
    for (volume_id, layer_id), samples in sorted(layer_samples.items()):
        radii = sorted(samples["r"])
        zs = sorted(samples["z"])
        r_spread = radii[-1] - radii[0]
        z_spread = zs[-1] - zs[0]
        median_r = radii[len(radii) // 2]
        median_z = zs[len(zs) // 2]
        is_disc = r_spread > z_spread
        layers.append({
            "id": f"{volume_id}/{layer_id}",
            "name": f"vol {volume_id} / layer {layer_id} "
                    + (f"(disc z≈{median_z:.0f} mm)" if is_disc else f"(r≈{median_r:.1f} mm)"),
            "kind": "disc" if is_disc else "cylinder",
            "sensitive": True,
            "r": round(median_r, 2), "rMin": round(radii[0], 2), "rMax": round(radii[-1], 2),
            "z": round(median_z, 2), "zMin": round(zs[0], 2), "zMax": round(zs[-1], 2),
        })
    return layers


def build_passive_layers(geo_config_path):
    with open(geo_config_path) as file:
        config = json.load(file)
    layers = []
    for cylinder in config.get("passiveCylinders", []):
        z_centre = cylinder.get("zCentre", 0.0)
        layers.append({"id": cylinder["name"], "name": f"{cylinder['name']} (r={cylinder['r']} mm)",
                       "kind": "cylinder", "sensitive": False, "r": cylinder["r"],
                       "zMin": z_centre - cylinder["halfZ"], "zMax": z_centre + cylinder["halfZ"]})
    for cylinder in config.get("forwardCylinders", []):
        layers.append({"id": cylinder["name"], "name": f"{cylinder['name']} (r={cylinder['r']} mm)",
                       "kind": "cylinder", "sensitive": False, "r": cylinder["r"],
                       "zMin": cylinder["zMin"], "zMax": cylinder["zMax"]})
    for disc in config.get("passiveDiscs", []):
        layers.append({"id": disc["name"], "name": f"{disc['name']} (z={disc['z']} mm)",
                       "kind": "disc", "sensitive": False, "z": disc["z"],
                       "rMin": disc["rMin"], "rMax": disc["rMax"]})
    return layers


def read_tracks(trackstates_path, tracksummary_path):
    """Per event: list of tracks with majority particle, classification, measurements and smoothed states."""
    summary = {}
    if tracksummary_path.exists():
        file = ROOT.TFile.Open(str(tracksummary_path))
        summary_tree = file.Get("tracksummary")
        summary_tree.SetBranchStatus("*", 0)
        for branch in ["event_nr", "track_nr", "nMeasurements", "t_pT", "trackClassification",
                       *[f"majorityParticleId_{field}" for field in BARCODE_FIELDS]]:
            summary_tree.SetBranchStatus(branch, 1)
        for entry in progress(summary_tree, summary_tree.GetEntries(), "tracksummary"):
            for i_track in range(len(entry.track_nr)):
                barcode = tuple(int(getattr(entry, f"majorityParticleId_{field}")[i_track])
                                for field in BARCODE_FIELDS)
                summary[(int(entry.event_nr), int(entry.track_nr[i_track]))] = {
                    "barcode": barcode,
                    "nMeas": int(entry.nMeasurements[i_track]),
                    "pt": float(entry.t_pT[i_track]),
                    # ActsExamples::TrackMatchClassification: 0 unknown, 1 matched, 2 duplicate, 3 fake
                    "cls": int(entry.trackClassification[i_track]),
                }
        file.Close()

    tracks = defaultdict(list)
    file = ROOT.TFile.Open(str(trackstates_path))
    states_tree = file.Get("trackstates")
    states_tree.SetBranchStatus("*", 0)
    for branch in ["event_nr", "track_nr", "nMeasurements", "volume_id", "layer_id", "g_x_hit", "g_y_hit", "g_z_hit",
                   "g_x_smt", "g_y_smt", "g_z_smt"]:
        states_tree.SetBranchStatus(branch, 1)
    for entry in progress(states_tree, states_tree.GetEntries(), "trackstates"):
        event_id, track_nr = int(entry.event_nr), int(entry.track_nr)
        measurements = {"x": [], "y": [], "z": [], "layer": []}
        smoothed = {"x": [], "y": [], "z": []}
        for i_state in range(len(entry.g_x_hit)):
            hit = (float(entry.g_x_hit[i_state]), float(entry.g_y_hit[i_state]), float(entry.g_z_hit[i_state]))
            if is_valid(*hit):
                for axis, value in zip("xyz", hit):
                    measurements[axis].append(value)
                measurements["layer"].append(f"{int(entry.volume_id[i_state])}/{int(entry.layer_id[i_state])}")
            if i_state < len(entry.g_x_smt):
                state = (float(entry.g_x_smt[i_state]), float(entry.g_y_smt[i_state]), float(entry.g_z_smt[i_state]))
                if is_valid(*state):
                    for axis, value in zip("xyz", state):
                        smoothed[axis].append(value)
        info = summary.get((event_id, track_nr))
        tracks[event_id].append({
            "nr": track_nr,
            "key": barcode_key(info["barcode"]) if info else None,
            "cls": info["cls"] if info else 0,
            "nMeas": info["nMeas"] if info else int(entry.nMeasurements),
            "pt": round(info["pt"], 4) if info else None,
            "meas": {axis: round_list(measurements[axis]) for axis in "xyz"},
            "measLayer": measurements["layer"],  # used for the hit matching, not written out
            "smt": {axis: round_list(values) for axis, values in smoothed.items()},
        })
    file.Close()
    return tracks


def read_matching(path):
    """Per (event, barcode): the matched flag of a matchingdetails tree, or None without the tree."""
    if not path.exists():
        return None
    file = ROOT.TFile.Open(str(path))
    details = file.Get("matchingdetails")
    if not details:
        file.Close()
        return None
    matched = {}
    for entry in details:
        key = (int(entry.event_nr),
               tuple(int(getattr(entry, f"particle_id_{field}")) for field in BARCODE_FIELDS))
        matched[key] = bool(entry.matched)
    file.Close()
    return matched


def read_seed_matching(run_dir):
    """Per (event, barcode): the seeding iterations whose seeds matched the particle.

    Only the per-particle matching booleans are written (matchingdetails), not the seeds' space
    points, so the display can say whether a particle was seeded but cannot draw the seed.
    """
    iteration_files = sorted(run_dir.glob("iter_*/*performance_seeding.root"),
                             key=lambda path: int(path.parent.name.split("_")[1]))
    seeded_in = None
    for path in iteration_files:
        iteration = int(path.parent.name.split("_")[1])
        matched = read_matching(path)
        if matched is None:
            continue
        seeded_in = seeded_in if seeded_in is not None else {}
        for key, is_matched in matched.items():
            seeded_in.setdefault(key, [])
            if is_matched:
                seeded_in[key].append(iteration)
    return seeded_in


def match_hits_to_tracks(particle_hits, particle_tracks):
    """Which simulated hits of a particle its tracks picked up, and which track measurements are its own.

    The trackstates do not say which simulated hit a measurement came from, so each measurement
    is paired with the nearest unpaired hit on the same layer within HIT_MATCH_DISTANCE, closest
    pairs first. Returns per hit 1/0 (found by a track or not) and per track a 1/0 list over its
    measurements (own hit or not, i.e. from another particle or noise).
    """
    n_hits = len(particle_hits.get("x", []))
    pairs = []
    for i_track, track in enumerate(particle_tracks):
        for i_meas, layer in enumerate(track["measLayer"]):
            point = (track["meas"]["x"][i_meas], track["meas"]["y"][i_meas], track["meas"]["z"][i_meas])
            for i_hit in range(n_hits):
                if particle_hits["layer"][i_hit] != layer:
                    continue
                hit = (particle_hits["x"][i_hit], particle_hits["y"][i_hit], particle_hits["z"][i_hit])
                distance = math.dist(point, hit)
                if distance < HIT_MATCH_DISTANCE:
                    pairs.append((distance, i_track, i_meas, i_hit))
    found = [0] * n_hits
    own = [[0] * len(track["measLayer"]) for track in particle_tracks]
    # duplicate tracks may share a hit, so a hit is only used once per track
    used_hits = set()
    for _, i_track, i_meas, i_hit in sorted(pairs):
        if own[i_track][i_meas] or (i_track, i_hit) in used_hits:
            continue
        own[i_track][i_meas] = 1
        used_hits.add((i_track, i_hit))
        found[i_hit] = 1
    return found, own


def distance_to_ray(point, origin, direction):
    """Distance from point to the half-line origin + t * direction, t >= 0."""
    offset = [point[i_axis] - origin[i_axis] for i_axis in range(3)]
    norm = math.sqrt(sum(component * component for component in direction)) or 1.0
    unit = [component / norm for component in direction]
    along = max(0.0, sum(offset[i_axis] * unit[i_axis] for i_axis in range(3)))
    return math.dist(point, [origin[i_axis] + along * unit[i_axis] for i_axis in range(3)])


def assign_parents(particles, particle_hits):
    """Parent key per particle, and whether it had to be inferred.

    Geant4 gives a secondary the parent's barcode with generation + 1 and a fresh
    sub-particle index, so the parent's own sub-particle index is lost. Generation 1
    always has a unique parent (the primary); from generation 2 on, several particles
    of the previous generation can qualify. Then the parent is the candidate produced
    before the child whose trajectory passes closest to the child's production vertex:
    its hits and vertex if it is charged, its straight line if it is neutral.
    """
    candidates = defaultdict(list)
    for part in particles:
        barcode = part["barcode"]
        candidates[(barcode[0], barcode[1], barcode[2], barcode[3])].append(part)

    parents = {}
    for part in particles:
        vertex_primary, vertex_secondary, particle_nr, generation, _ = part["barcode"]
        if generation == 0:
            continue
        options = candidates.get((vertex_primary, vertex_secondary, particle_nr, generation - 1), [])
        if len(options) == 1:
            parents[part["barcode"]] = (options[0]["barcode"], False)
            continue
        if not options:
            parents[part["barcode"]] = (None, False)  # parent not stored (e.g. it left no hits)
            continue
        child_vertex = (part["vx"], part["vy"], part["vz"])

        def score(option):
            produced_later = option["vt"] > part["vt"] + 1e-6
            origin = (option["vx"], option["vy"], option["vz"])
            if option["q"] != 0:
                option_hits = particle_hits.get(option["barcode"], {})
                points = [origin] + list(zip(option_hits.get("x", []), option_hits.get("y", []),
                                             option_hits.get("z", [])))
                distance = min(math.dist(child_vertex, point) for point in points)
            else:
                distance = distance_to_ray(child_vertex, origin, (option["px"], option["py"], option["pz"]))
            return (produced_later, distance)

        parents[part["barcode"]] = (min(options, key=score)["barcode"], True)
    return parents


def decay_tree_particles(particles, event_hits):
    """Particles of one event with their (possibly inferred) parent, for the decay tree."""
    parents = assign_parents(particles, event_hits)
    particle_list = []
    for part in particles:
        parent_barcode, is_inferred = parents.get(part["barcode"], (None, False))
        particle_list.append({
            "key": barcode_key(part["barcode"]),
            "parent": barcode_key(parent_barcode) if parent_barcode else None,
            "inferred": is_inferred,
            "generation": part["barcode"][3],
            "pdg": part["pdg"], "q": part["q"], "process": part["process"],
            "p": round(part["p"], 6), "tkin": round(part["tkin"], 6), "eloss": round(part["eloss"], 6),
            "vx": round(part["vx"], 2), "vy": round(part["vy"], 2), "vz": round(part["vz"], 2),
            "vt": round(part["vt"], 4),
            "nHits": len(event_hits.get(part["barcode"], {}).get("x", [])),
        })
    particle_list.sort(key=lambda entry: (entry["generation"], entry["vt"]))
    return particle_list


def build_events(particles_by_event, hits, tracks, seed_matching, reco_matching, with_decay_tree=True):
    """One entry per event: the particles with hits (or in the efficiency selection) with their
    hits and status, the tracks, the production vertices of no-hit secondaries and the decay tree."""
    events = []
    for event_id in progress(sorted(particles_by_event), len(particles_by_event), "events"):
        arrays = particles_by_event[event_id]
        event_hits = hits.get(event_id, {})
        barcodes = arrays["barcode"]
        is_shown = [barcode in event_hits
                    or (reco_matching is not None and (event_id, barcode) in reco_matching)
                    or (seed_matching is not None and (event_id, barcode) in seed_matching)
                    for barcode in barcodes]
        if not any(is_shown):
            continue
        event_tracks = tracks.get(event_id, [])
        tracks_of = defaultdict(list)
        for track in event_tracks:
            if track["key"] is not None and track["cls"] != 3:  # fakes do not count as finding a hit
                tracks_of[track["key"]].append(track)

        parts = []
        for i_part in np.flatnonzero(is_shown):
            part = particle_dict(arrays, i_part)
            barcode = part["barcode"]
            particle_hits = event_hits.get(barcode, {})
            match_key = (event_id, barcode)
            in_selection = ((reco_matching is not None and match_key in reco_matching)
                            or (seed_matching is not None and match_key in seed_matching))
            found, own = match_hits_to_tracks(particle_hits, tracks_of.get(barcode_key(barcode), []))
            for track, own_flags in zip(tracks_of.get(barcode_key(barcode), []), own):
                track["own"] = own_flags
            parts.append({
                "k": barcode_key(barcode),
                "fam": family_key(barcode),
                "g": barcode[3],
                "pdg": part["pdg"], "q": part["q"], "proc": part["process"],
                "pt": round(part["pt"], 5), "eta": round(part["eta"], 4), "phi": round(part["phi"], 4),
                "px": round(part["px"], 5), "py": round(part["py"], 5), "pz": round(part["pz"], 5),
                "p": round(part["p"], 5), "tkin": round(part["tkin"], 5), "eloss": round(part["eloss"], 5),
                "vx": round(part["vx"], 2), "vy": round(part["vy"], 2), "vz": round(part["vz"], 2),
                "vt": round(part["vt"], 4),
                # Geant4 sets the momentum to zero when a track ends through a physics process
                # (stopping, decay, inelastic interaction); a kill by time or volume leaves it
                "atRest": abs(part["eloss"] - part["tkin"]) < 1e-3 * part["tkin"],
                "sel": in_selection,  # in the efficiency denominator (matchingdetails)
                # None: unknown (no matching details, or not selected for the efficiency)
                "seed": seed_matching.get(match_key) if seed_matching is not None else None,
                "reco": reco_matching.get(match_key) if reco_matching is not None else None,
                "nHits": len(particle_hits.get("x", [])),
                "hits": {
                    "x": round_list(particle_hits.get("x", [])), "y": round_list(particle_hits.get("y", [])),
                    "z": round_list(particle_hits.get("z", [])), "t": round_list(particle_hits.get("t", []), 4),
                    "p": round_list(particle_hits.get("p", []), 5), "layer": particle_hits.get("layer", []),
                    "found": found,  # 1: a (non-fake) track of this particle has a measurement on it
                },
            })

        # production vertices of secondaries that left no hits (drawn together with their family)
        no_hit_secondaries = [i_part for i_part in np.flatnonzero(arrays["generation"] > 0) if not is_shown[i_part]]
        no_hit_secondaries.sort(key=lambda i_part: arrays["vt"][i_part])
        vertices = [{"fam": family_key(barcodes[i_part]), "pdg": int(arrays["particle_type"][i_part]),
                     "vx": round(float(arrays["vx"][i_part]), 2), "vy": round(float(arrays["vy"][i_part]), 2),
                     "vz": round(float(arrays["vz"][i_part]), 2), "p": round(float(arrays["p"][i_part]), 5),
                     "proc": int(arrays["process"][i_part])}
                    for i_part in no_hit_secondaries[:MAX_VERTICES_PER_EVENT]]

        # the full particle list (decay tree) dominates the file size; None leaves the tree out
        tree = None
        is_pruned = with_decay_tree and len(barcodes) > DECAY_TREE_MAX_PARTICLES
        if with_decay_tree:
            tree_indices = np.flatnonzero(is_shown) if is_pruned else range(len(barcodes))
            tree = decay_tree_particles([particle_dict(arrays, i_part) for i_part in tree_indices], event_hits)

        events.append({
            "id": event_id,
            "nSim": len(barcodes),
            "parts": parts,
            "tracks": [{key: value for key, value in track.items() if key != "measLayer"} for track in event_tracks],
            "vertices": vertices,
            "nVerticesDropped": max(0, len(no_hit_secondaries) - MAX_VERTICES_PER_EVENT),
            "tree": tree,
            "treePruned": is_pruned,
        })
    return events


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("output_dir", type=Path, help="ACTSO2 output directory of one run")
    parser.add_argument("-o", "--output", type=Path, default=None, help="HTML file to write")
    parser.add_argument("--geo-config", type=Path, default=None,
                        help="gen3_geometry_config.json, to add the passive cylinders and discs as layers")
    parser.add_argument("--bfield", type=float, default=2.0, help="B_z in T used for the ideal helix")
    parser.add_argument("--no-decay-tree", action="store_true",
                        help="leave out the per-event particle list behind the decay tree (most of the file size)")
    parser.add_argument("--max-events", type=int, default=0, help="keep only the first N events (0 = all)")
    args = parser.parse_args()

    run_dir = args.output_dir
    particles_by_event = read_particles(run_dir / "particles_simulation.root")
    if args.max_events:
        kept = sorted(particles_by_event)[:args.max_events]
        particles_by_event = {event_id: particles_by_event[event_id] for event_id in kept}
    hits, layer_samples = read_hits(run_dir / "hits.root")
    trackstates_path = run_dir / "trackstates_ambi-tracks-merged.root"
    tracks = read_tracks(trackstates_path, run_dir / "tracksummary_ambi-tracks-merged.root") \
        if trackstates_path.exists() else {}

    seed_matching = read_seed_matching(run_dir)
    if seed_matching is None:
        print("no iter_*/*performance_seeding.root with matchingdetails: seeding status not shown", file=sys.stderr)
    reco_matching = read_matching(run_dir / "performance_merged_ambi_tracks.root")
    if reco_matching is None:
        print("no matchingdetails in performance_merged_ambi_tracks.root: reconstruction status not shown",
              file=sys.stderr)
    events = build_events(particles_by_event, hits, tracks, seed_matching, reco_matching,
                          with_decay_tree=not args.no_decay_tree)
    layers = build_sensitive_layers(layer_samples)
    if args.geo_config:
        layers += build_passive_layers(args.geo_config)

    payload = {"sample": run_dir.name, "bField": args.bfield, "layers": layers, "events": events,
               "hasSeedStatus": seed_matching is not None, "hasRecoStatus": reco_matching is not None}
    template = (Path(__file__).parent / "event_display_template.html").read_text()
    html = template.replace("/*__EVENT_DATA__*/null", json.dumps(payload, separators=(",", ":")))

    output = args.output or Path(__file__).parent / "figures" / f"{run_dir.name}_event_display.html"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html)
    n_particles = sum(len(event["parts"]) for event in events)
    print(f"{len(events)} events, {n_particles} particles, {len(layers)} layers -> {output} "
          f"({output.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
