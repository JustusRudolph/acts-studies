#!/usr/bin/env python3
"""Build an interactive HTML event display from an ACTSO2 output directory.

Reads (all optional except hits.root and particles_simulation.root):
    hits.root                                  simulated hits (sensitive layers)
    particles_simulation.root                  Geant4 particles: primaries + kept secondaries
    trackstates_ambi-tracks-merged.root        reconstructed tracks (measurements + smoothed states)
    tracksummary_ambi-tracks-merged.root       majority particle and classification of each track
    performance_merged_ambi_tracks.root        per-particle reconstruction matching (matchingdetails)
    iter_N/*performance_seeding.root           per-particle seed matching of each iteration
    iter_N/seeds.root, *estimatedparameters.root   the seeds (RootSeedWriter, writingMode "full") and their parameters
    iter_N/seed_segments_belowTOF.root         the seeds whose helix stays below the TOF (TofSeedExtensionAlgorithm),
                                               which are not in seeds.root
    measurements_with_noise.root               the particles behind each measurement, to assign the seeds
    iter_N/trackstates*ambi_tracks_writer.root the measurements iteration N removes (cfg.tracking.writeAmbiTrackStates),
                                               for the "after iteration" filter
    iter_N/propagation_steps_*ckf_tracks_*_steps.root   the CKF propagation steps (first/second pass, succeeded
                                               and failed, failed extrapolations; writePropagationSteps), under
                                               the seed they started from

and writes a single self-contained HTML file with x-y, r-z and 3D views, layer toggles,
and an event > primaries/secondaries > particle browser with filters on the particles'
number of hits, eta, pT and seeding/reconstruction status.

    # in the ACTSO2/latest environment with ACTSO2's setup.sh sourced
    python3 /home/justus/projects/alice/acts-studies/visualise/make_event_display.py \
        /home/justus/projects/alice/ACTSO2/output/test_loopers_10000ev \
        --geometry /home/justus/projects/alice/ACTSO2/geometries/geometry_2026_10_02_segmented_realistic/geom

With --geometry the sensitive layers come from the ACTS tracking geometry built from
o2sim_geometry.root (as alice3_full_chain.py does) and the passive ones from
gen3_geometry_config.json; without it the sensitive layers are estimated from the hits.

With --sim-only only hits.root and particles_simulation.root are read: all simulated particles
with hits are shown, without tracks, seeds or seeding/reconstruction status.

With --reco-dir the reconstruction files (tracks, matching, seeds, measurements) are read from
that directory instead, while hits and particles still come from output_dir. This compares
several reconstructions of one simulation, e.g. reco_only runs on the same sim run.

The output goes to figures/<sample>_event_display.html (gitignored) unless -o is given; with
--reco-dir the name is <sim dir>__<reco dir>.
"""

import argparse
import json
import math
import os
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path


# Needs ACTSO2's setup.sh sourced (alice3, acts and AliceActsPythonBindings on PYTHONPATH).
# Restarts with a static TLS reserve if the aliBuild Geant4 is built with the initial-exec TLS model.
from alice3.static_tls import ensure_static_tls_reserve

ensure_static_tls_reserve()

# The acts bindings (with Geant4) are loaded before ROOT, and acts.examples before AliceActsPythonBindings,
# since otherwise the static TLS reserve is used up first. They are only needed for --geometry, so a failed
# import is kept and raised there.
try:
    import acts
    import acts.examples  # noqa: F401
    from AliceActsPythonBindings import buildALICE3GeometryGen3
    acts_import_error = None
except ImportError as error:
    acts_import_error = error

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
# RootTrackStatesWriter::StateType: 0 measurement, 1 outlier, 2 hole, 3 material
STATE_MEASUREMENT = 0
# particles read from particles_simulation.root at a time (PbPb with all secondaries: millions per event)
PARTICLE_CHUNK = 5_000_000
# output above this size gets a warning: browsers struggle to open it
LARGE_OUTPUT_MB = 300


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
        return f"{label:<22} |{bar}| {100 * fraction:5.1f}%  {count}/{total}  {elapsed:5.1f}s  eta {eta:5.1f}s"

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


_n_event_sets = 0


def filter_events(frame, column, event_ids):
    """The RDataFrame restricted to the given event ids (None: all events)."""
    global _n_event_sets
    if event_ids is None:
        return frame
    _n_event_sets += 1
    name = f"keptEventIds{_n_event_sets}"
    ROOT.gInterpreter.Declare(f"#include <unordered_set>\nstd::unordered_set<long long> {name};")
    kept = getattr(ROOT, name)
    for event_id in event_ids:
        kept.insert(int(event_id))
    return frame.Filter(f"{name}.count({column}) > 0")


_tree_checks = {}


def has_tree(path, name):
    """Whether path is a readable ROOT file with a TTree called name. A missing file is skipped
    silently; an existing one without the tree (e.g. created but never written, as when a run
    stops in that iteration) once with a note, so every optional input that cannot be read is
    simply left out of the display."""
    key = (str(path), name)
    if key not in _tree_checks:
        is_readable = False
        if Path(path).exists():
            level, ROOT.gErrorIgnoreLevel = ROOT.gErrorIgnoreLevel, ROOT.kError  # no "has no keys" warnings
            file = ROOT.TFile.Open(str(path))
            ROOT.gErrorIgnoreLevel = level
            if file and not file.IsZombie():
                tree = file.Get(name)
                is_readable = bool(tree) and tree.InheritsFrom("TTree")
            if file:
                file.Close()
            if not is_readable:
                print(f"{path}: no '{name}' tree (empty or unfinished file), skipped", file=sys.stderr)
        _tree_checks[key] = is_readable
    return _tree_checks[key]


def tree_entries(tree, event_branch, event_ids, label):
    """Yield the tree loaded with each entry of the kept events; other entries only read the event number."""
    branch, leaf = tree.GetBranch(event_branch), tree.GetLeaf(event_branch)
    for i_entry in progress(range(tree.GetEntries()), tree.GetEntries(), label):
        if event_ids is not None:
            branch.GetEntry(i_entry)
            if int(leaf.GetValue()) not in event_ids:
                continue
        tree.GetEntry(i_entry)
        yield tree


def read_event_ids(path):
    """The event ids in particles_simulation.root, in file order."""
    return [int(event_id) for event_id in ROOT.RDataFrame("particles", str(path)).AsNumpy(["event_id"])["event_id"]]


def reduce_particles(arrays, needed, with_decay_tree):
    """The particles of one event that the display uses: the needed ones (with hits or in the matching),
    the earliest MAX_VERTICES_PER_EVENT other secondaries (production-vertex markers), and all of
    them for an unpruned decay tree. Keeps the original order; nSim and nNoHitSecondaries count all."""
    n_sim = len(arrays["particle"])
    is_needed = np.zeros(n_sim, dtype=bool)
    if needed:
        # cheap pre-selection on the particle number, then the full barcode
        numbers = np.fromiter({barcode[2] for barcode in needed}, dtype=np.int64)
        for i_part in np.flatnonzero(np.isin(arrays["particle"].astype(np.int64), numbers)):
            is_needed[i_part] = tuple(int(arrays[field][i_part]) for field in BARCODE_FIELDS) in needed
    no_hit_secondaries = np.flatnonzero((arrays["generation"] > 0) & ~is_needed)
    if with_decay_tree and n_sim <= DECAY_TREE_MAX_PARTICLES:
        keep = np.ones(n_sim, dtype=bool)
    else:
        # earliest first, ties in file order (as build_events sorts them)
        earliest = no_hit_secondaries[np.lexsort((no_hit_secondaries, arrays["vt"][no_hit_secondaries]))]
        keep = is_needed.copy()
        keep[earliest[:MAX_VERTICES_PER_EVENT]] = True
    reduced = {column: np.array(values[keep]) for column, values in arrays.items()}
    reduced["barcode"] = list(zip(*(reduced[field].tolist() for field in BARCODE_FIELDS)))
    reduced["nSim"] = n_sim
    reduced["nNoHitSecondaries"] = len(no_hit_secondaries)
    return reduced


def read_particles(path, event_ids, needed_of, with_decay_tree):
    """Per event: the columns of the particles the display uses (reduce_particles) as numpy arrays.

    PbPb events with all Geant4 secondaries have millions of particles, so the file is read in
    chunks of about PARTICLE_CHUNK particles and each event is reduced before the next is read.
    """
    columns = ["event_id", "particle_type", "q", "m", "p", "pt", "eta", "phi", "px", "py", "pz",
               "vx", "vy", "vz", "vt", "e_loss", "process", *BARCODE_FIELDS]
    sizes = filter_events(ROOT.RDataFrame("particles", str(path)), "event_id", event_ids) \
        .Define("i_entry", "static_cast<long long>(rdfentry_)") \
        .Define("n_particles", "static_cast<long long>(particle_type.size())") \
        .AsNumpy(["i_entry", "n_particles"])
    chunks, chunk, n_in_chunk = [], [], 0
    for i_entry, n_particles in zip(sizes["i_entry"].tolist(), sizes["n_particles"].tolist()):
        if chunk and n_in_chunk + n_particles > PARTICLE_CHUNK:
            chunks.append(chunk)
            chunk, n_in_chunk = [], 0
        chunk.append(i_entry)
        n_in_chunk += n_particles
    if chunk:
        chunks.append(chunk)

    particles_by_event = {}
    bar = progress(range(len(sizes["i_entry"])), len(sizes["i_entry"]), "particles (events)")
    for chunk in chunks:
        frame = ROOT.RDataFrame("particles", str(path)).Range(chunk[0], chunk[-1] + 1)
        data = filter_events(frame, "event_id", event_ids).AsNumpy(columns)
        for i_event in range(len(data["event_id"])):
            event_id = int(data["event_id"][i_event])
            arrays = {column: np.asarray(data[column][i_event]) for column in columns[1:]}
            arrays["vt"] = arrays["vt"] / MM_PER_NS
            arrays["tkin"] = np.sqrt(arrays["p"] ** 2 + arrays["m"] ** 2) - arrays["m"]
            particles_by_event[event_id] = reduce_particles(arrays, needed_of.get(event_id, set()), with_decay_tree)
            next(bar, None)
        del data
    for _ in bar:
        pass
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


def read_hits(path, event_ids=None):
    """Per event and barcode: hit arrays, plus the r and z of all hits per layer (of the kept events)."""
    columns = ["event_id", "tx", "ty", "tz", "tt", "tpx", "tpy", "tpz", "deltapx", "deltapy", "deltapz",
               "volume_id", "layer_id",
               *[f"barcode_{field}" for field in BARCODE_FIELDS]]
    data = filter_events(ROOT.RDataFrame("hits", str(path)), "event_id", event_ids).AsNumpy(columns)
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
        # momentum entering the sensor (momentum4Before: pre-step point of the first Geant4 step in
        # it) and its change up to leaving it (momentum4After - momentum4Before), in GeV
        for axis in "xyz":
            particle_hits[f"p{axis}"].append(float(f"{data[f'tp{axis}'][i_hit]:.5g}"))
            particle_hits[f"dp{axis}"].append(float(f"{data[f'deltap{axis}'][i_hit]:.4g}"))
        particle_hits["layer"].append(f"{layer_key[0]}/{layer_key[1]}")
        # all hits: hits.root is sorted by module, so a subset would only cover part of a layer
        samples = layer_samples[layer_key]
        samples["r"].append(math.hypot(x, y))
        samples["z"].append(z)
    # hits.root is ordered by detector id; number each particle's hits in time order instead
    for event_hits in hits.values():
        for particle_hits in event_hits.values():
            time_order = sorted(range(len(particle_hits["t"])), key=lambda i_hit: particle_hits["t"][i_hit])
            for column, values in particle_hits.items():
                particle_hits[column] = [values[i_hit] for i_hit in time_order]
    return hits, layer_samples


def build_geometry_layers(geom_dir):
    """Sensitive layers from the ACTS tracking geometry, keyed by the same volume/layer ids as the hits.

    The geometry is built like alice3_full_chain.py does for gen3, with the material layers
    (withMaterial=True), which take part in the layer numbering. The vertex-detector layers are
    single cylinder surfaces; the other layers are made of planar modules, whose corners and edge
    midpoints give the layer's radial and z extent.
    """
    if acts_import_error is not None:
        raise acts_import_error

    geometry = buildALICE3GeometryGen3(str(geom_dir / "o2sim_geometry.root"), acts.logging.WARNING, "", True,
                                       str(geom_dir / "gen3_geometry_config.json"))
    gctx = acts.GeometryContext.dangerouslyDefaultConstruct()
    surfaces = []
    geometry.visitSurfaces(lambda surface: surfaces.append(surface))

    extents = defaultdict(lambda: {"r": [], "z": [], "centreR": [], "cylinder": None})
    for surface in progress(surfaces, len(surfaces), "surfaces"):
        if not surface.isSensitive:
            continue
        geometry_id = surface.geometryId
        extent = extents[(geometry_id.volume, geometry_id.layer)]
        bounds = list(surface.bounds.values())
        if str(surface.type) == "SurfaceType.Cylinder":
            centre_z = surface.center(gctx)[2]
            extent["cylinder"] = (bounds[0], centre_z - bounds[1], centre_z + bounds[1])
            continue
        # rectangle bounds: (min x, min y, max x, max y) in the local frame
        x_lo, y_lo, x_hi, y_hi = bounds[:4]
        x_mid, y_mid = (x_lo + x_hi) / 2, (y_lo + y_hi) / 2
        for local in [(x_lo, y_lo), (x_lo, y_hi), (x_hi, y_lo), (x_hi, y_hi),
                      (x_mid, y_lo), (x_mid, y_hi), (x_lo, y_mid), (x_hi, y_mid)]:
            point = surface.localToGlobal(gctx, acts.Vector2(*local), acts.Vector3(0, 0, 1))
            extent["r"].append(math.hypot(point[0], point[1]))
            extent["z"].append(point[2])
        centre = surface.center(gctx)
        extent["centreR"].append(math.hypot(centre[0], centre[1]))

    layers = []
    for (volume_id, layer_id), extent in sorted(extents.items()):
        layer = {"id": f"{volume_id}/{layer_id}", "sensitive": True}
        if extent["cylinder"]:
            radius, z_min, z_max = extent["cylinder"]
            layer.update(kind="cylinder", r=round(radius, 2), rMin=round(radius, 2), rMax=round(radius, 2),
                         z=round((z_min + z_max) / 2, 2), zMin=round(z_min, 2), zMax=round(z_max, 2))
        else:
            radii, zs = np.asarray(extent["r"]), np.asarray(extent["z"])
            is_disc = radii.max() - radii.min() > zs.max() - zs.min()
            layer.update(kind="disc" if is_disc else "cylinder",
                         # barrel: modules sit at slightly different radii; draw the median centre radius
                         r=round(float(np.median(extent["centreR"])), 2),
                         rMin=round(float(radii.min()), 2), rMax=round(float(radii.max()), 2),
                         z=round(float((zs.min() + zs.max()) / 2), 2),
                         zMin=round(float(zs.min()), 2), zMax=round(float(zs.max()), 2))
        layer["name"] = f"vol {volume_id} / layer {layer_id} " + (
            f"(disc z≈{layer['z']:.0f} mm, r {layer['rMin']:.0f}–{layer['rMax']:.0f} mm)" if layer["kind"] == "disc"
            else f"(r≈{layer['r']:.1f} mm, |z| ≤ {max(abs(layer['zMin']), abs(layer['zMax'])):.0f} mm)")
        layers.append(layer)
    return layers


def build_sensitive_layers(layer_samples):
    """Cylinders (barrel) or discs (endcap), estimated from where the hits are (no geometry given)."""
    layers = []
    for (volume_id, layer_id), samples in sorted(layer_samples.items()):
        radii = np.asarray(samples["r"])
        zs = np.asarray(samples["z"])
        r_spread = radii.max() - radii.min()
        z_spread = zs.max() - zs.min()
        median_r = float(np.median(radii))
        median_z = float(np.median(zs))
        is_disc = r_spread > z_spread
        # the detector is symmetric in z, a sample (a forward gun) may not be: barrel layers
        # span +-max|z| of their hits
        z_half = float(np.abs(zs).max())
        layers.append({
            "id": f"{volume_id}/{layer_id}",
            "name": f"vol {volume_id} / layer {layer_id} "
                    + (f"(disc z≈{median_z:.0f} mm)" if is_disc else f"(r≈{median_r:.1f} mm)"),
            "kind": "disc" if is_disc else "cylinder",
            "sensitive": True,
            "r": round(median_r, 2), "rMin": round(float(radii.min()), 2), "rMax": round(float(radii.max()), 2),
            "z": round(median_z, 2),
            "zMin": round(float(zs.min()) if is_disc else -z_half, 2),
            "zMax": round(float(zs.max()) if is_disc else z_half, 2),
        })
    # discs seen on one side only get a mirrored copy on the other (no hits there in this sample)
    discs = [layer for layer in layers if layer["kind"] == "disc"]
    for disc in discs:
        # the two sensor faces of a disc are a few mm apart, so compare centres generously
        centre = (disc["zMin"] + disc["zMax"]) / 2
        has_mirror = any(abs((other["zMin"] + other["zMax"]) / 2 + centre) < 20
                         and abs(other["rMax"] - disc["rMax"]) < 0.2 * disc["rMax"] for other in discs)
        if not has_mirror:
            layers.append({**disc, "id": f"{disc['id']} mirrored", "z": -disc["z"],
                           "zMin": -disc["zMax"], "zMax": -disc["zMin"],
                           "name": f"{disc['name']} mirrored to z≈{-disc['z']:.0f} mm (no hits there)"})
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


def read_tracks(trackstates_path, tracksummary_path, event_ids=None):
    """Per event: list of tracks with majority particle, classification, measurements and smoothed states."""
    summary = {}
    if has_tree(tracksummary_path, "tracksummary"):
        file = ROOT.TFile.Open(str(tracksummary_path))
        summary_tree = file.Get("tracksummary")
        summary_tree.SetBranchStatus("*", 0)
        for branch in ["event_nr", "track_nr", "nMeasurements", "nMajorityHits", "t_pT", "trackClassification",
                       *[f"majorityParticleId_{field}" for field in BARCODE_FIELDS]]:
            summary_tree.SetBranchStatus(branch, 1)
        for entry in tree_entries(summary_tree, "event_nr", event_ids, "tracksummary (events)"):
            for i_track in range(len(entry.track_nr)):
                barcode = tuple(int(getattr(entry, f"majorityParticleId_{field}")[i_track])
                                for field in BARCODE_FIELDS)
                summary[(int(entry.event_nr), int(entry.track_nr[i_track]))] = {
                    "barcode": barcode,
                    "nMeas": int(entry.nMeasurements[i_track]),
                    "nMaj": int(entry.nMajorityHits[i_track]),
                    "pt": float(entry.t_pT[i_track]),
                    # ActsExamples::TrackMatchClassification: 0 unknown, 1 matched, 2 duplicate, 3 fake
                    "cls": int(entry.trackClassification[i_track]),
                }
        file.Close()

    tracks = defaultdict(list)
    file = ROOT.TFile.Open(str(trackstates_path))
    states_tree = file.Get("trackstates")
    states_tree.SetBranchStatus("*", 0)
    id_branches = [f"particle_ids_{field}" for field in BARCODE_FIELDS]
    has_particle_ids = all(states_tree.GetBranch(branch) for branch in id_branches + ["stateType"])
    for branch in ["event_nr", "track_nr", "nMeasurements", "volume_id", "layer_id", "g_x_hit", "g_y_hit", "g_z_hit",
                   "g_x_smt", "g_y_smt", "g_z_smt", *(id_branches + ["stateType"] if has_particle_ids else [])]:
        states_tree.SetBranchStatus(branch, 1)
    for entry in tree_entries(states_tree, "event_nr", event_ids, "trackstates (tracks)"):
        event_id, track_nr = int(entry.event_nr), int(entry.track_nr)
        measurements = {"x": [], "y": [], "z": [], "layer": []}
        smoothed = {"x": [], "y": [], "z": []}
        # per particle: its hits behind the track's measurements, as the ACTS truth matcher counts
        # them (TrackClassification.cpp: measurement states only, outliers skipped)
        hit_counts = Counter()
        id_columns = [getattr(entry, branch) for branch in id_branches] if has_particle_ids else None
        for i_state in range(len(entry.g_x_hit)):
            hit = (float(entry.g_x_hit[i_state]), float(entry.g_y_hit[i_state]), float(entry.g_z_hit[i_state]))
            if is_valid(*hit):
                for axis, value in zip("xyz", hit):
                    measurements[axis].append(value)
                measurements["layer"].append(f"{int(entry.volume_id[i_state])}/{int(entry.layer_id[i_state])}")
                if id_columns is not None and entry.stateType[i_state] == STATE_MEASUREMENT:
                    hit_counts.update(zip(*(map(int, column[i_state]) for column in id_columns)))
            if i_state < len(entry.g_x_smt):
                state = (float(entry.g_x_smt[i_state]), float(entry.g_y_smt[i_state]), float(entry.g_z_smt[i_state]))
                if is_valid(*state):
                    for axis, value in zip("xyz", state):
                        smoothed[axis].append(value)
        info = summary.get((event_id, track_nr))
        # the majority particle: from the summary, which leaves it empty (all-zero) for fakes;
        # then from the particles behind the track's measurements in the trackstates
        majority, n_majority = None, None
        if info and any(info["barcode"]):
            majority, n_majority = info["barcode"], info["nMaj"]
        elif hit_counts:
            # most measurements first, ties to the first in the track
            majority, n_majority = hit_counts.most_common(1)[0]
        tracks[event_id].append({
            "nr": track_nr,
            "key": barcode_key(majority) if majority else None,
            "nMaj": n_majority,  # measurements from the majority particle
            "cls": info["cls"] if info else 0,
            "nMeas": info["nMeas"] if info else int(entry.nMeasurements),
            "pt": round(info["pt"], 4) if info else None,
            "meas": {axis: round_list(measurements[axis]) for axis in "xyz"},
            "measLayer": measurements["layer"],  # used for the hit matching, not written out
            "smt": {axis: round_list(values) for axis, values in smoothed.items()},
        })
    file.Close()
    return tracks


def iteration_number(directory):
    """3 for iter_3; None for other names."""
    suffix = directory.name.removeprefix("iter_")
    return int(suffix) if suffix.isdigit() else None


def read_removed_measurements(run_dir, event_ids=None):
    """The measurements each iteration removes before the next: its ambiguity-resolved tracks'
    measurement states (MeasurementFilterAlgorithm skips outliers), in
    iter_<n>/trackstates_*ambi_tracks_writer.root (cfg.tracking.writeAmbiTrackStates).

    Returns the iterations covered (0..n-1 without gaps) and per (event, barcode) the list of
    (iteration, layer, x, y, z) of the measurements of that particle; None without the files.
    """
    iteration_dirs = sorted((path for path in run_dir.glob("iter_*") if iteration_number(path) is not None),
                            key=iteration_number)
    files = {iteration_number(path.parent): path
             for path in run_dir.glob("iter_*/trackstates*ambi_tracks_writer.root") if has_tree(path, "trackstates")}
    iterations = []
    for directory in iteration_dirs:
        if iteration_number(directory) not in files:
            break
        iterations.append(iteration_number(directory))
    if not iterations:
        return None, None
    if len(iterations) < len(iteration_dirs):
        print(f"no ambi trackstates for iteration {len(iterations)}: 'after iteration' only up to "
              f"{iterations[-1]}", file=sys.stderr)
    removed = defaultdict(list)
    id_branches = [f"particle_ids_{field}" for field in BARCODE_FIELDS]
    for iteration in iterations:
        file = ROOT.TFile.Open(str(files[iteration]))
        tree = file.Get("trackstates")
        tree.SetBranchStatus("*", 0)
        for branch in ["event_nr", "volume_id", "layer_id", "stateType", "g_x_hit", "g_y_hit", "g_z_hit", *id_branches]:
            tree.SetBranchStatus(branch, 1)
        for entry in tree_entries(tree, "event_nr", event_ids, f"iteration {iteration} trackstates"):
            event_id = int(entry.event_nr)
            id_columns = [getattr(entry, branch) for branch in id_branches]
            for i_state in range(len(entry.stateType)):
                if entry.stateType[i_state] != STATE_MEASUREMENT:
                    continue
                hit = (float(entry.g_x_hit[i_state]), float(entry.g_y_hit[i_state]), float(entry.g_z_hit[i_state]))
                if not is_valid(*hit):
                    continue
                layer = f"{int(entry.volume_id[i_state])}/{int(entry.layer_id[i_state])}"
                # a measurement from several particles (merged cluster) is removed for each of them
                for barcode in zip(*(map(int, column[i_state]) for column in id_columns)):
                    removed[(event_id, barcode)].append((iteration, layer, *hit))
        file.Close()
    return iterations, removed


def removing_iterations(particle_hits, measurements):
    """Per hit of a particle: the iteration whose tracks removed its measurement, or -1.

    Each measurement is paired with the nearest unpaired hit on its layer within
    HIT_MATCH_DISTANCE, closest pairs first (as match_hits_to_tracks).
    """
    n_hits = len(particle_hits.get("x", []))
    removed_in = [-1] * n_hits
    if not measurements:
        return removed_in
    # tracks of one iteration may share a measurement (maximumSharedHits): pair it once
    measurements = list(dict.fromkeys(measurements))
    pairs = []
    for i_meas, (_, layer, *point) in enumerate(measurements):
        for i_hit in range(n_hits):
            if particle_hits["layer"][i_hit] != layer:
                continue
            hit = (particle_hits["x"][i_hit], particle_hits["y"][i_hit], particle_hits["z"][i_hit])
            distance = math.dist(point, hit)
            if distance < HIT_MATCH_DISTANCE:
                pairs.append((distance, i_meas, i_hit))
    used_measurements = set()
    for _, i_meas, i_hit in sorted(pairs):
        if i_meas in used_measurements or removed_in[i_hit] >= 0:
            continue
        used_measurements.add(i_meas)
        removed_in[i_hit] = measurements[i_meas][0]
    return removed_in


def read_matching(path, event_ids=None):
    """Per (event, barcode): the matched flag of a matchingdetails tree, or None without the tree."""
    if not has_tree(path, "matchingdetails"):
        return None
    file = ROOT.TFile.Open(str(path))
    details = file.Get("matchingdetails")
    if not details:
        file.Close()
        return None
    matched = {}
    for entry in details:
        if event_ids is not None and int(entry.event_nr) not in event_ids:
            continue
        key = (int(entry.event_nr),
               tuple(int(getattr(entry, f"particle_id_{field}")) for field in BARCODE_FIELDS))
        matched[key] = bool(entry.matched)
    file.Close()
    return matched


def read_seed_matching(run_dir, event_ids=None):
    """Per (event, barcode): the seeding iterations whose seeds matched the particle (matchingdetails)."""
    iteration_files = sorted(run_dir.glob("iter_*/*performance_seeding.root"),
                             key=lambda path: int(path.parent.name.split("_")[1]))
    seeded_in = None
    for path in iteration_files:
        iteration = int(path.parent.name.split("_")[1])
        matched = read_matching(path, event_ids)
        if matched is None:
            continue
        seeded_in = seeded_in if seeded_in is not None else {}
        for key, is_matched in matched.items():
            seeded_in.setdefault(key, [])
            if is_matched:
                seeded_in[key].append(iteration)
    return seeded_in


SEGMENTS_FILE = "seed_segments_belowTOF.root"

# iter_N/propagation_steps_[iter_N_]ckf_tracks_[failed_]<kind>_steps.root (RootPropagationStepsWriter)
STEPS_FILE = re.compile(r"^propagation_steps_(?:iter_\d+_)?ckf_tracks_(failed_)?(firstPass|secondPass|extrapolation)_steps\.root$")
STEP_KINDS = {"firstPass": "FP", "secondPass": "SP", "extrapolation": "EX"}
# a second pass starts on a measurement of its first pass, an extrapolation on a state of the
# track: the propagation is linked to the seed whose recorded steps pass closest to its start
LINK_DISTANCE = 5.0  # mm
LINK_TIE = 0.01  # mm: seeds this much further away than the closest are as good a candidate
STEP_COLUMNS = ["g_x", "g_y", "g_z", "type", "step_acc", "volume_id", "layer_id", "sensitive_id", "nStepTrials"]


def read_step_file(path, event_ids):
    """One row per recorded propagation (segment) of a RootPropagationStepsWriter file."""
    frame = filter_events(ROOT.RDataFrame("propagation_steps", str(path)), "event_nr", event_ids)
    # seed_index: written since PropagationSummary::seedIndex (-1: not from a seed)
    has_seed_index = "seed_index" in {str(column) for column in frame.GetColumnNames()}
    data = frame.AsNumpy(["event_nr", "d0", "z0", "phi", "theta", "qOverP", *STEP_COLUMNS,
                          *(["seed_index"] if has_seed_index else [])])
    for column in STEP_COLUMNS:
        data[column] = [np.asarray(values) for values in data[column]]
    return data


def step_payload(data, i_row, max_steps):
    """The steps of one propagation for the display; above max_steps only the first and last
    max_steps / 2 (a propagation can run into maxSteps), "gap" says where and how many are left out."""
    n_steps = len(data["g_x"][i_row])
    keep = slice(None)
    gap = None
    if max_steps and n_steps > max_steps:
        half = max_steps // 2
        keep = np.r_[0:half, n_steps - half:n_steps]
        gap = [half, n_steps - 2 * half]

    def column(name, ndigits=None):
        values = data[name][i_row][keep]
        return [round(float(value), ndigits) for value in values] if ndigits is not None else values.astype(int).tolist()

    return {"n": n_steps, **({"gap": gap} if gap else {}),
            "x": column("g_x", 2), "y": column("g_y", 2), "z": column("g_z", 2),
            "h": column("type"), "acc": [float(f"{value:.3g}") for value in data["step_acc"][i_row][keep]],
            "v": column("volume_id"), "l": column("layer_id"), "s": column("sensitive_id"), "tr": column("nStepTrials")}


def read_propagations(run_dir, event_ids=None, max_steps=20000):
    """Per event: {seed ref "iteration:index in event": [propagation, ...]}; None without step files.

    With a seed_index branch (PropagationSummary::seedIndex) every propagation names its seed.
    Older files are linked as follows. A first pass starts from the seed's estimated parameters (estimatedparams, one row per seed in
    seeds.root order), so it is linked to its seed exactly. Second passes and extrapolations start
    on the first-pass (or second-pass) trajectory: they go to the seed with recorded steps closest
    to their start, the earliest seed not before the previous one of the same file when several are
    as close (the algorithm processes the seeds in order); "amb" lists the other candidates.
    """
    props_of = defaultdict(lambda: defaultdict(list))
    n_unlinked = Counter()
    found_any = False
    for directory in sorted(run_dir.glob("iter_*"), key=lambda path: iteration_number(path) or 0):
        iteration = iteration_number(directory)
        if iteration is None:
            continue
        step_files = {}
        for path in directory.glob("propagation_steps_*.root"):
            match = STEPS_FILE.match(path.name)
            if match and has_tree(path, "propagation_steps"):
                step_files[(STEP_KINDS[match.group(2)], 0 if match.group(1) else 1)] = path
        if not step_files:
            continue
        found_any = True
        param_files = [path for path in directory.glob("*estimatedparameters.root") if has_tree(path, "estimatedparams")]
        if not param_files:
            print(f"{directory}: no estimatedparameters.root, the CKF steps cannot be linked to seeds", file=sys.stderr)
            continue
        params = filter_events(ROOT.RDataFrame("estimatedparams", str(param_files[0])), "event_nr", event_ids).AsNumpy(
            ["event_nr", "loc0", "loc1", "phi", "theta", "qop"])
        seed_params = defaultdict(list)  # event -> rows (index in event = seed index)
        for i_row, event_id in enumerate(params["event_nr"].tolist()):
            seed_params[event_id].append([float(params[name][i_row]) for name in ("loc0", "loc1", "phi", "theta", "qop")])
        # recorded step positions per (event, seed index), to link the later propagations
        trajectory_points = defaultdict(lambda: defaultdict(list))

        def add_propagations(kind, is_ok, data, links):
            for i_row, seed_index, segment, n_segments, ambiguous in links:
                event_id = int(data["event_nr"][i_row])
                props_of[event_id][f"{iteration}:{seed_index}"].append({
                    "kind": kind, "ok": is_ok, "br": segment, "nBr": n_segments,
                    **({"amb": ambiguous} if ambiguous else {}), **step_payload(data, i_row, max_steps)})
                trajectory_points[event_id][seed_index].append(
                    np.stack([data["g_x"][i_row], data["g_y"][i_row], data["g_z"][i_row]], axis=1))

        def segment_groups(data):
            """Rows of one propagation (the CKF branches, same start parameters, written one after the other)."""
            groups = []
            for i_row in range(len(data["event_nr"])):
                start = (int(data["event_nr"][i_row]), *(float(data[name][i_row]) for name in ("d0", "z0", "phi", "theta", "qOverP")))
                if groups and groups[-1][0] == start:
                    groups[-1][1].append(i_row)
                else:
                    groups.append((start, [i_row]))
            return groups

        # first passes: the seed with exactly these start parameters
        for is_ok in (1, 0):
            path = step_files.get(("FP", is_ok))
            if path is None:
                continue
            data = read_step_file(path, event_ids)
            links = []
            for (event_id, *start), rows in segment_groups(data):
                if "seed_index" in data:
                    linked = [(i_row, int(data["seed_index"][i_row])) for i_row in rows]
                    links += [(i_row, seed_index, i_segment, len(rows), None)
                              for i_segment, (i_row, seed_index) in enumerate(linked) if seed_index >= 0]
                    n_unlinked["FP"] += sum(seed_index < 0 for _, seed_index in linked)
                    continue
                candidates = seed_params.get(event_id, [])
                deviations = [max(abs(a - b) / max(1e-6, abs(b)) for a, b in zip(start, row)) for row in candidates]
                if not deviations or min(deviations) > 1e-4:
                    n_unlinked["FP"] += len(rows)
                    continue
                seed_index = int(np.argmin(deviations))
                links += [(i_row, seed_index, i_segment, len(rows), None) for i_segment, i_row in enumerate(rows)]
            add_propagations("FP", is_ok, data, links)

        # second passes, then extrapolations: the seed whose trajectory passes through the start
        for kind, is_ok in (("SP", 1), ("SP", 0), ("EX", 0)):
            path = step_files.get((kind, is_ok))
            if path is None:
                continue
            data = read_step_file(path, event_ids)
            links = []
            previous = {}
            for (event_id, *_), rows in segment_groups(data):
                if "seed_index" in data:
                    linked = [(i_row, int(data["seed_index"][i_row])) for i_row in rows]
                    links += [(i_row, seed_index, i_segment, len(rows), None)
                              for i_segment, (i_row, seed_index) in enumerate(linked) if seed_index >= 0]
                    n_unlinked[kind] += sum(seed_index < 0 for _, seed_index in linked)
                    continue
                first = rows[0]
                start = np.array([data[axis][first][0] for axis in ("g_x", "g_y", "g_z")]) if len(data["g_x"][first]) else None
                distances = {seed_index: min(float(np.min(np.linalg.norm(points - start, axis=1))) for points in point_sets)
                             for seed_index, point_sets in trajectory_points[event_id].items()} if start is not None else {}
                if not distances or min(distances.values()) > LINK_DISTANCE:
                    n_unlinked[kind] += len(rows)
                    continue
                closest = min(distances.values())
                candidates = sorted(seed_index for seed_index, distance in distances.items() if distance <= closest + LINK_TIE)
                later = [seed_index for seed_index in candidates if seed_index >= previous.get(event_id, -1)]
                seed_index = later[0] if later else candidates[0]
                previous[event_id] = seed_index
                ambiguous = [other for other in candidates if other != seed_index] or None
                links += [(i_row, seed_index, i_segment, len(rows), ambiguous) for i_segment, i_row in enumerate(rows)]
            add_propagations(kind, is_ok, data, links)
    if not found_any:
        return None
    if n_unlinked:
        print(f"CKF steps not linked to a seed: {dict(n_unlinked)}", file=sys.stderr)
    return {event_id: dict(props) for event_id, props in props_of.items()}


def read_segments(path, event_ids=None):
    """The below-TOF seeds of iter_N/seed_segments_belowTOF.root (RootSeedSegmentWriter), one row per seed."""
    columns = ["event_nr", "seed_index", "partner_index", "pair_same_particle", "pt", "cot_theta",
               "sp_x", "sp_y", "sp_z", "measurement_id"]
    data = filter_events(ROOT.RDataFrame("segments", str(path)), "event_nr", event_ids).AsNumpy(columns)
    for column in ("sp_x", "sp_y", "sp_z", "measurement_id"):
        data[column] = [list(values) for values in data[column]]
    return data


def seed_measurement_keys(run_dir, event_ids=None):
    """(event, measurement id) of the space points of all seeds (iter_N/seeds.root, writingMode "full",
    and the below-TOF seeds of iter_N/seed_segments_belowTOF.root)."""
    keys = set()
    for path in run_dir.glob(f"iter_*/{SEGMENTS_FILE}"):
        if not has_tree(path, "segments"):
            continue
        segments = read_segments(path, event_ids)
        for event_id, measurement_ids in zip(segments["event_nr"].tolist(), segments["measurement_id"]):
            keys.update((event_id, int(measurement_id)) for measurement_id in measurement_ids)
    for path in run_dir.glob("iter_*/seeds.root"):
        if not has_tree(path, "seeds"):
            continue
        frame = ROOT.RDataFrame("seeds", str(path))
        columns = ["event_id", *[f"measurement_id_{i_sp}" for i_sp in (1, 2, 3)]]
        if not set(columns) <= {str(column) for column in frame.GetColumnNames()}:
            continue
        data = filter_events(frame, "event_id", event_ids).AsNumpy(columns)
        for i_sp in (1, 2, 3):
            keys.update(zip(data["event_id"].tolist(), data[f"measurement_id_{i_sp}"].tolist()))
    return keys


_n_measurement_sets = 0


def read_measurement_particles(run_dir, event_ids=None, wanted_keys=None):
    """Per (event, measurement id): the barcodes of the particles behind the measurement (none for noise).

    With wanted_keys only those measurements are read: the file has tens of millions of
    measurements for PbPb, of which only the seeds' space points are needed.
    """
    path = next((path for name in ("measurements_with_noise.root", "measurements.root")
                 if has_tree(path := run_dir / name, "measurements")), None)
    if path is None:
        return None
    columns = ["event_nr", "measurement_id", *[f"particles_{field}" for field in BARCODE_FIELDS]]
    frame = filter_events(ROOT.RDataFrame("measurements", str(path)), "event_nr", event_ids)
    if wanted_keys is not None:
        global _n_measurement_sets
        _n_measurement_sets += 1
        name = f"wantedMeasurements{_n_measurement_sets}"
        ROOT.gInterpreter.Declare(f"#include <unordered_set>\nstd::unordered_set<unsigned long long> {name};")
        wanted = getattr(ROOT, name)
        wanted.reserve(len(wanted_keys))
        for event_id, measurement_id in wanted_keys:
            wanted.insert((int(event_id) << 32) | int(measurement_id))
        frame = frame.Filter(f"{name}.count((static_cast<unsigned long long>(event_nr) << 32) "
                             f"| static_cast<unsigned long long>(measurement_id)) > 0")
    data = frame.AsNumpy(columns)
    particles_of = {}
    for i_meas in progress(range(len(data["event_nr"])), len(data["event_nr"]), "measurements"):
        barcodes = zip(*(data[f"particles_{field}"][i_meas] for field in BARCODE_FIELDS))
        particles_of[(int(data["event_nr"][i_meas]), int(data["measurement_id"][i_meas]))] = \
            {tuple(int(part) for part in barcode) for barcode in barcodes}
    return particles_of


def read_seeds(run_dir, measurement_particles, event_ids=None, orphan_refs=None):
    """Per (event, barcode): the seeds with at least two of their three space points from the particle.

    The seeds come from iter_N/seeds.root (RootSeedWriter, writingMode "full" for the space
    points, quality and z_vertex), the estimated pT and eta from the estimatedparams tree,
    which has one row per seed in the same order. The particles behind each space point come
    from the measurements file; without it only seeds whose hits all belong to one particle
    (truthMatched in estimatedparams) are assigned.

    Every seed has a "ref" (iteration:index in its event). The seeds under no particle whose ref
    is in orphan_refs[event] (those with CKF steps) are kept under (event, None).
    """
    seeds_of = defaultdict(list)
    orphan_refs = orphan_refs or {}
    seed_files = sorted((path for path in run_dir.glob("iter_*/seeds.root") if has_tree(path, "seeds")),
                        key=lambda path: int(path.parent.name.split("_")[1]))
    for path in seed_files:
        iteration = int(path.parent.name.split("_")[1])
        frame = ROOT.RDataFrame("seeds", str(path))
        if "x_1" not in [str(column) for column in frame.GetColumnNames()]:
            print(f"{path} has no space points (writingMode \"small\"): skipped", file=sys.stderr)
            continue
        columns = ["event_id", "seed_quality", "z_vertex",
                   *[f"{column}_{i_sp}" for column in ("measurement_id", "x", "y", "z") for i_sp in (1, 2, 3)]]
        seeds = filter_events(frame, "event_id", event_ids).AsNumpy(columns)
        n_seeds = len(seeds["event_id"])
        params = None
        param_files = [param_path for param_path in path.parent.glob("*estimatedparameters.root")
                       if has_tree(param_path, "estimatedparams")]
        if param_files:
            # one row per seed in the same order, so the same event filter keeps them lined up
            params = filter_events(ROOT.RDataFrame("estimatedparams", str(param_files[0])), "event_nr", event_ids).AsNumpy(
                ["event_nr", "pt", "eta", "truthMatched", *[f"particleId_{field}" for field in BARCODE_FIELDS]])
            if not np.array_equal(params["event_nr"], seeds["event_id"]):
                print(f"{param_files[0]} does not line up with {path}: no estimated parameters", file=sys.stderr)
                params = None
        index_in_event = Counter()
        for i_seed in progress(range(n_seeds), n_seeds, f"seeds iter {iteration}"):
            event_id = int(seeds["event_id"][i_seed])
            ref = f"{iteration}:{index_in_event[event_id]}"
            index_in_event[event_id] += 1
            if measurement_particles is not None:
                sp_particles = [measurement_particles.get((event_id, int(seeds[f"measurement_id_{i_sp}"][i_seed])), set())
                                for i_sp in (1, 2, 3)]
            elif params is not None and params["truthMatched"][i_seed]:
                barcode = tuple(int(params[f"particleId_{field}"][i_seed]) for field in BARCODE_FIELDS)
                sp_particles = [{barcode}] * 3
            else:
                sp_particles = [set()] * 3
            seed = {
                "it": iteration,
                "ref": ref,
                "q": round(float(seeds["seed_quality"][i_seed]), 1),
                "zv": round(float(seeds["z_vertex"][i_seed]), 1),
                **({"pt": round(float(params["pt"][i_seed]), 5), "eta": round(float(params["eta"][i_seed]), 3)}
                   if params is not None else {}),
                **{axis: [round(float(seeds[f"{axis}_{i_sp}"][i_seed]), 2) for i_sp in (1, 2, 3)] for axis in "xyz"},
            }
            n_owners = 0
            for barcode in set().union(*sp_particles):
                own = [int(barcode in particles) for particles in sp_particles]
                if sum(own) >= 2:
                    seeds_of[(event_id, barcode)].append({**seed, "own": own})
                    n_owners += 1
            if not n_owners and ref in orphan_refs.get(event_id, ()):
                seeds_of[(event_id, None)].append({**seed, "own": [int(bool(particles)) for particles in sp_particles]})
        add_segments(path.parent / SEGMENTS_FILE, iteration, measurement_particles, event_ids, seeds_of)
    return seeds_of if seed_files else None


def add_segments(path, iteration, measurement_particles, event_ids, seeds_of):
    """The below-TOF seeds of one iteration, under each particle with ≥ 2 of their space points,
    as read_seeds does for seeds.root (the writer's own majority particle is not used: it names one
    even when the space points all come from different particles)."""
    if not has_tree(path, "segments") or measurement_particles is None:
        if has_tree(path, "segments"):
            print(f"{path}: no measurements file to assign the space points, skipped", file=sys.stderr)
        return
    segments = read_segments(path, event_ids)
    n_segments = len(segments["event_nr"])
    for i_seed in progress(range(n_segments), n_segments, f"below-TOF seeds iter {iteration}"):
        event_id = int(segments["event_nr"][i_seed])
        measurement_ids = segments["measurement_id"][i_seed]
        sp_particles = [measurement_particles.get((event_id, int(measurement_id)), set())
                        for measurement_id in measurement_ids]
        if not sp_particles:
            continue
        partner = int(segments["partner_index"][i_seed])
        seed = {
            "it": iteration,
            "belowTof": 1,
            "nr": int(segments["seed_index"][i_seed]),
            # index of the paired segment (the other half of a loop) in the same event, -1: unpaired
            "partner": partner,
            **({"partnerSame": int(segments["pair_same_particle"][i_seed])} if partner >= 0 else {}),
            "pt": round(float(segments["pt"][i_seed]), 5),
            "eta": round(math.asinh(float(segments["cot_theta"][i_seed])), 3),
            **{axis: round_list(segments[f"sp_{axis}"][i_seed]) for axis in "xyz"},
        }
        for barcode in set().union(*sp_particles):
            own = [int(barcode in particles) for particles in sp_particles]
            if sum(own) >= 2:
                seeds_of[(event_id, barcode)].append({**seed, "own": own})


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


def build_events(particles_by_event, hits, tracks, seed_matching, reco_matching, seeds=None, with_decay_tree=True,
                 removed=None, propagations=None):
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
        event_props = (propagations or {}).get(event_id, {})
        if not any(is_shown) and not event_props:
            continue
        event_tracks = tracks.get(event_id, [])
        # every track under its majority particle, whatever its classification
        tracks_of = defaultdict(list)
        for track in event_tracks:
            if track["key"] is not None:
                tracks_of[track["key"]].append(track)

        parts = []
        for i_part in np.flatnonzero(is_shown):
            part = particle_dict(arrays, i_part)
            barcode = part["barcode"]
            particle_hits = event_hits.get(barcode, {})
            match_key = (event_id, barcode)
            in_selection = ((reco_matching is not None and match_key in reco_matching)
                            or (seed_matching is not None and match_key in seed_matching))
            # the pairing is per track: "own" for all tracks, "found" only from the non-fake ones
            particle_tracks = tracks_of.get(barcode_key(barcode), [])
            _, own = match_hits_to_tracks(particle_hits, particle_tracks)
            for track, own_flags in zip(particle_tracks, own):
                track["own"] = own_flags
            found, _ = match_hits_to_tracks(particle_hits, [track for track in particle_tracks if track["cls"] != 3])
            hit_columns = {}
            if removed is not None:
                # -1: still there after the last iteration
                hit_columns["it"] = removing_iterations(particle_hits, removed.get(match_key))
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
                "seeds": seeds.get(match_key, []) if seeds is not None else [],
                "nHits": len(particle_hits.get("x", [])),
                "hits": {
                    "x": round_list(particle_hits.get("x", [])), "y": round_list(particle_hits.get("y", [])),
                    "z": round_list(particle_hits.get("z", [])), "t": round_list(particle_hits.get("t", []), 4),
                    "p": round_list(particle_hits.get("p", []), 5), "layer": particle_hits.get("layer", []),
                    # momentum entering the sensor and its change up to leaving it (GeV)
                    **{column: particle_hits.get(column, []) for column in ("px", "py", "pz", "dpx", "dpy", "dpz")},
                    "found": found,  # 1: a (non-fake) track of this particle has a measurement on it
                    **hit_columns,
                },
            })

        # production vertices of secondaries that left no hits (drawn together with their family)
        # reduce_particles kept the earliest of them; the count is of all
        no_hit_secondaries = [i_part for i_part in np.flatnonzero(arrays["generation"] > 0) if not is_shown[i_part]]
        no_hit_secondaries.sort(key=lambda i_part: arrays["vt"][i_part])
        vertices = [{"fam": family_key(barcodes[i_part]), "pdg": int(arrays["particle_type"][i_part]),
                     "vx": round(float(arrays["vx"][i_part]), 2), "vy": round(float(arrays["vy"][i_part]), 2),
                     "vz": round(float(arrays["vz"][i_part]), 2), "p": round(float(arrays["p"][i_part]), 5),
                     "proc": int(arrays["process"][i_part])}
                    for i_part in no_hit_secondaries[:MAX_VERTICES_PER_EVENT]]

        # the full particle list (decay tree) dominates the file size; None leaves the tree out
        tree = None
        is_pruned = with_decay_tree and arrays["nSim"] > DECAY_TREE_MAX_PARTICLES
        if with_decay_tree:
            tree_indices = np.flatnonzero(is_shown) if is_pruned else range(len(barcodes))
            tree = decay_tree_particles([particle_dict(arrays, i_part) for i_part in tree_indices], event_hits)

        events.append({
            "id": event_id,
            "nSim": arrays["nSim"],
            "parts": parts,
            "tracks": [{key: value for key, value in track.items() if key != "measLayer"} for track in event_tracks],
            "vertices": vertices,
            "nVerticesDropped": max(0, arrays["nNoHitSecondaries"] - MAX_VERTICES_PER_EVENT),
            "tree": tree,
            "treePruned": is_pruned,
            # CKF propagations per seed ref, and the seeds with some but without a majority particle
            "props": event_props,
            "nseeds": seeds.get((event_id, None), []) if seeds is not None else [],
        })
    return events


def read_reconstruction(run_dir, sim_dir=None, event_ids=None, max_steps=20000):
    """Tracks, seed matching, reconstruction matching and seeds of a run; each is empty/None when missing.

    The measurements (the particles behind the seeds' space points) are taken from run_dir,
    or from sim_dir when run_dir has none.
    """
    trackstates_path = run_dir / "trackstates_ambi-tracks-merged.root"
    tracks = read_tracks(trackstates_path, run_dir / "tracksummary_ambi-tracks-merged.root", event_ids) \
        if has_tree(trackstates_path, "trackstates") else {}
    seed_matching = read_seed_matching(run_dir, event_ids)
    if seed_matching is None:
        print("no iter_*/*performance_seeding.root with matchingdetails: seeding status not shown", file=sys.stderr)
    reco_matching = read_matching(run_dir / "performance_merged_ambi_tracks.root", event_ids)
    if reco_matching is None:
        print("no matchingdetails in performance_merged_ambi_tracks.root: reconstruction status not shown",
              file=sys.stderr)
    # the measurements only serve to assign the seeds' space points
    measurement_particles = None
    if any(has_tree(path, "seeds") for path in run_dir.glob("iter_*/seeds.root")):
        wanted_keys = seed_measurement_keys(run_dir, event_ids)
        measurement_particles = read_measurement_particles(run_dir, event_ids, wanted_keys)
        if measurement_particles is None and sim_dir is not None:
            measurement_particles = read_measurement_particles(sim_dir, event_ids, wanted_keys)
    propagations = read_propagations(run_dir, event_ids, max_steps)
    if propagations is None:
        print("no iter_*/propagation_steps_*.root (writePropagationSteps): no CKF steps", file=sys.stderr)
    seeds = read_seeds(run_dir, measurement_particles, event_ids,
                       orphan_refs={event_id: set(props) for event_id, props in (propagations or {}).items()})
    if seeds is None:
        print("no iter_*/seeds.root: seeds not drawn", file=sys.stderr)
    iterations, removed = read_removed_measurements(run_dir, event_ids)
    if iterations is None:
        print("no iter_*/trackstates*ambi_tracks_writer.root (cfg.tracking.writeAmbiTrackStates): "
              "no 'after iteration' filter", file=sys.stderr)
    return tracks, seed_matching, reco_matching, seeds, iterations, removed, propagations


def recorded_sim_dir(reco_dir):
    """The simulation a reco_only run read (its --simdir, else its -f), as logged in
    configs/script_command_line_args_used.log; None when the run did not log a reco_only call.
    The path is as given on the command line, relative to the directory the chain ran in."""
    log = reco_dir / "configs" / "script_command_line_args_used.log"
    if not log.exists():
        return None
    calls = [line.split() for line in log.read_text().splitlines() if "reco_only" in line]
    if not calls:
        return None
    call = calls[-1]
    for flag in ("--simdir", "-f"):
        if flag in call and call.index(flag) + 1 < len(call):
            return call[call.index(flag) + 1]
    return None


def is_same_run_dir(logged, run_dir):
    """Whether a logged (relative) run path names run_dir: run_dir's absolute path ends with it."""
    logged_parts = [part for part in Path(logged).parts if part not in (".", "")]
    return list(run_dir.resolve().parts[-len(logged_parts):]) == logged_parts


def parse_event_list(text, parser):
    """Event ids from "3,7,10-12"."""
    event_ids = set()
    for item in text.split(","):
        item = item.strip()
        try:
            if "-" in item:
                first, last = (int(bound) for bound in item.split("-"))
                event_ids.update(range(first, last + 1))
            elif item:
                event_ids.add(int(item))
        except ValueError:
            parser.error(f"--events: cannot read {item!r} (expected e.g. 3,7,10-12)")
    return event_ids


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("output_dir", type=Path, help="ACTSO2 output directory of one run")
    parser.add_argument("-o", "--output", type=Path, default=None, help="HTML file to write")
    parser.add_argument("--geometry", "--geo-config", type=Path, default=None,
                        help="geometry directory with o2sim_geometry.root and gen3_geometry_config.json "
                             "(or a file in it): sensitive layers from the ACTS geometry, passive ones from the config")
    parser.add_argument("--bfield", type=float, default=2.0, help="B_z in T (for the bending radius in the particle details)")
    parser.add_argument("--no-decay-tree", action="store_true",
                        help="leave out the per-event particle list behind the decay tree (most of the file size)")
    parser.add_argument("--max-events", type=int, default=0, help="keep only the first N events (0 = all)")
    parser.add_argument("--events", type=str, default=None,
                        help="only these event ids, e.g. 3,7,10-12 (with --max-events: the first N of them)")
    parser.add_argument("--sim-only", "--sim_only", action="store_true",
                        help="only the simulation (hits.root, particles_simulation.root): no tracks, seeds "
                             "or matching status")
    parser.add_argument("--reco-dir", "--reco_dir", "--reco-loc", "--reco_loc", type=Path, default=None,
                        help="directory with the reconstruction output (tracks, matching, seeds) of the "
                             "simulation in output_dir; default: output_dir itself")
    parser.add_argument("--max-propagation-steps", type=int, default=20000,
                        help="keep at most this many steps of a CKF propagation, the first and last half (0 = all)")
    args = parser.parse_args()
    if args.sim_only and args.reco_dir:
        parser.error("--sim-only and --reco-dir exclude each other")

    run_dir = args.output_dir
    if not run_dir.is_dir():
        parser.error(f"output directory {run_dir} does not exist")
    reco_dir = args.reco_dir or run_dir
    if not reco_dir.is_dir():
        parser.error(f"reconstruction directory {reco_dir} does not exist")

    # the events first: every reader then skips the others
    particles_path = run_dir / "particles_simulation.root"
    all_event_ids = sorted(read_event_ids(particles_path))
    kept_event_ids = all_event_ids
    if args.events:
        wanted = parse_event_list(args.events, parser)
        missing = sorted(wanted - set(all_event_ids))
        if missing:
            print(f"warning: events {missing} are not in {particles_path}", file=sys.stderr)
        kept_event_ids = [event_id for event_id in all_event_ids if event_id in wanted]
    if args.max_events:
        kept_event_ids = kept_event_ids[:args.max_events]
    if not kept_event_ids:
        parser.error("no event selected")
    event_ids = None if len(kept_event_ids) == len(all_event_ids) else set(kept_event_ids)
    print(f"{len(kept_event_ids)} of {len(all_event_ids)} events", file=sys.stderr)

    hits, layer_samples = read_hits(run_dir / "hits.root", event_ids)
    if args.sim_only:
        tracks, seed_matching, reco_matching, seeds, iterations, removed, propagations = {}, None, None, None, None, None, None
    else:
        if args.reco_dir:
            print(f"reconstruction from {reco_dir}", file=sys.stderr)
        tracks, seed_matching, reco_matching, seeds, iterations, removed, propagations = read_reconstruction(
            reco_dir, sim_dir=run_dir, event_ids=event_ids, max_steps=args.max_propagation_steps)
        logged_sim_dir = recorded_sim_dir(reco_dir) if args.reco_dir else None
        if logged_sim_dir and not is_same_run_dir(logged_sim_dir, run_dir):
            print(f"warning: {reco_dir} was reconstructed from {logged_sim_dir} (its configs/"
                  f"script_command_line_args_used.log), not {run_dir}", file=sys.stderr)
    # particles last, keeping only those with hits or in the matching (and the vertex markers)
    needed_of = defaultdict(set)
    for event_id, event_hits in hits.items():
        needed_of[event_id].update(event_hits)
    for matching in (reco_matching, seed_matching):
        for event_id, barcode in matching or {}:
            needed_of[event_id].add(barcode)
    particles_by_event = read_particles(particles_path, event_ids, needed_of, with_decay_tree=not args.no_decay_tree)
    events = build_events(particles_by_event, hits, tracks, seed_matching, reco_matching, seeds,
                          with_decay_tree=not args.no_decay_tree, removed=removed, propagations=propagations)
    if args.geometry:
        geom_dir = args.geometry if args.geometry.is_dir() else args.geometry.parent
        layers = build_geometry_layers(geom_dir) + build_passive_layers(geom_dir / "gen3_geometry_config.json")
        layer_source = "geometry"
    else:
        print("no --geometry: sensitive layers are estimated from the hits", file=sys.stderr)
        layers = build_sensitive_layers(layer_samples)
        layer_source = "hits"

    sample = f"{run_dir.resolve().name}__{reco_dir.resolve().name}" if args.reco_dir else run_dir.resolve().name
    payload = {"sample": sample, "bField": args.bfield, "layers": layers, "layerSource": layer_source,
               "events": events,
               "hasSeedStatus": seed_matching is not None, "hasRecoStatus": reco_matching is not None,
               "hasSeeds": seeds is not None, "simOnly": args.sim_only,
               "iterations": iterations or [], "hasSteps": propagations is not None}
    template = (Path(__file__).parent / "event_display_template.html").read_text()
    html = template.replace("/*__EVENT_DATA__*/null", json.dumps(payload, separators=(",", ":")))

    output = args.output or Path(__file__).parent / "figures" / f"{sample}_event_display.html"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html)
    n_particles = sum(len(event["parts"]) for event in events)
    size_mb = output.stat().st_size / 1e6
    print(f"{len(events)} events, {n_particles} particles, {len(layers)} layers -> {output} ({size_mb:.1f} MB)")
    if size_mb > LARGE_OUTPUT_MB:
        print(f"warning: {size_mb:.0f} MB is a lot for a browser; fewer events (--events, --max-events) "
              f"or --no-decay-tree make it smaller", file=sys.stderr)


if __name__ == "__main__":
    main()
