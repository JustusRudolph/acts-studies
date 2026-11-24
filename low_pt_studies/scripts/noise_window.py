"""Expected noise hits in the CKF search window of every measurement on a track.

For each measurement of the tracks in trackstates_ambi-tracks-merged.root (ACTSO2
RootTrackStatesWriter) the window is the region the CKF chi2 cut accepts around
the prediction, an ellipse of area pi * chi2_cut * sigma_loc0 * sigma_loc1 (the
loc0-loc1 correlation is not stored and taken as 0, an upper bound of the area).
Times the noise density of the surface (noise probability per pixel / pixel
area, from the digi config) this is N, the expected number of noise hits in the
window: the quantity NoiseWindowTrackStateCreator cuts on.

Every measurement is classified by truth (own = majority particle of the track,
other particle, noise = no particle) and by stage (the first three measurements
in time are the seed, the others were added by the CKF), and the script prints
per region the distribution of N, and for the CKF hits how much noise and how
many own hits a cut on N would remove.

Usage (needs ROOT, e.g. via alienv):
    python3 scripts/noise_window.py <run dir> <digi config json> [--chi2-cut 45] [--noise-prob 1e-7]
"""

import argparse
import collections
import json
import math

import ROOT

# volumes with injected per-pixel noise (cfg.digi.noiseLayersByVolume in ACTSO2, gen 3)
NOISE_VOLUMES = {20, 30, 40, 50, 60, 70}
ITOF = (40, 17)


def pixel_pitches(digi_file):
    """(volume, layer) -> (pitch loc0, pitch loc1) in mm from the digi config."""
    with open(digi_file) as digi_input:
        digi = json.load(digi_input)
    entries = digi["entries"] if isinstance(digi, dict) and "entries" in digi else digi
    pitches = {}
    for entry in entries:
        binning = entry.get("value", {}).get("geometric", {}).get("segmentation", {}).get("binningdata", [])
        entry_pitches = [(axis["max"] - axis["min"]) / axis["bins"] for axis in binning if axis.get("bins")]
        if len(entry_pitches) == 2:
            pitches[(entry.get("volume"), entry.get("layer"))] = entry_pitches
    return pitches


def region_name(volume, layer):
    if (volume, layer) == ITOF:
        return "iTOF"
    return {40: "ML", 50: "OT", 90: "oTOF"}.get(volume, f"vol{volume}")


def quantile(values, fraction):
    ordered = sorted(values)
    return ordered[int(fraction * (len(ordered) - 1))]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run_dir", help="ACTSO2 output directory with trackstates_ambi-tracks-merged.root")
    parser.add_argument("digi_file", help="digitisation config json of the run")
    parser.add_argument("--chi2-cut", type=float, default=45., help="CKF measurement chi2 cut (ckfChi2Measurement)")
    parser.add_argument("--noise-prob", type=float, default=1e-7, help="noise probability per pixel (noiseProbPerPixel)")
    args = parser.parse_args()

    pitches = pixel_pitches(args.digi_file)
    states_file = ROOT.TFile.Open(args.run_dir + "/trackstates_ambi-tracks-merged.root")
    states = states_file.Get("trackstates")

    rows = []  # (region, stage, truth class, expected noise N)
    for track in states:
        n_states = len(track.chi2)
        particle_of = {}
        for i_state in range(n_states):
            if len(track.particle_ids_sub_particle[i_state]):
                particle_of[i_state] = (track.particle_ids_vertex_primary[i_state][0],
                                        track.particle_ids_generation[i_state][0],
                                        track.particle_ids_sub_particle[i_state][0])
        if not particle_of:
            continue
        majority = collections.Counter(particle_of.values()).most_common(1)[0][0]
        measurements = [i_state for i_state in range(n_states)
                        if (track.chi2[i_state] > 0 or i_state in particle_of)
                        and math.isfinite(track.err_eLOC0_prt[i_state])]
        # states are stored from the tip to the stem: the last three are the seed
        seed_states = set(measurements[-3:])
        for i_state in measurements:
            surface = (track.volume_id[i_state], track.layer_id[i_state])
            if surface not in pitches:
                continue
            pitch_loc0, pitch_loc1 = pitches[surface]
            density = args.noise_prob / (pitch_loc0 * pitch_loc1) if surface[0] in NOISE_VOLUMES else 0.
            area = math.pi * args.chi2_cut * track.err_eLOC0_prt[i_state] * track.err_eLOC1_prt[i_state]
            if i_state not in particle_of:
                truth_class = "noise"
            elif particle_of[i_state] == majority:
                truth_class = "own"
            else:
                truth_class = "other"
            rows.append((region_name(*surface), "seed" if i_state in seed_states else "CKF", truth_class, density * area))

    groups = collections.defaultdict(list)
    for region, stage, truth_class, expected_noise in rows:
        groups[(region, stage, truth_class)].append(expected_noise)
    print("region stage class      n   median N      90% N   frac N>1  frac N>0.1")
    for key in sorted(groups):
        values = groups[key]
        print("%-6s %-5s %-6s %5d %10.3g %10.3g %10.2f %10.2f" % (
            key + (len(values), quantile(values, .5), quantile(values, .9),
                   sum(value > 1 for value in values) / len(values),
                   sum(value > .1 for value in values) / len(values))))

    own = [row[3] for row in rows if row[1] == "CKF" and row[2] == "own"]
    noise = [row[3] for row in rows if row[1] == "CKF" and row[2] == "noise"]
    if own and noise:
        for cut in (0.1, 0.3, 1, 3):
            n_noise = sum(value > cut for value in noise)
            n_own = sum(value > cut for value in own)
            print("CKF hits, cut N > %-4g: removes %4d/%4d noise (%3.0f%%), %4d/%4d own (%3.0f%%)" % (
                cut, n_noise, len(noise), 100 * n_noise / len(noise), n_own, len(own), 100 * n_own / len(own)))


if __name__ == "__main__":
    main()
