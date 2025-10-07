import uproot
import argparse

import numpy as np
import pandas as pd

import variables
import constants

def is_in_reduced_volume(x, y, z):
    return (constants.RminX < x < constants.RmaxX) and (constants.RminY < y < constants.RmaxY) and (constants.RminZ < z < constants.RmaxZ)

def compute_track_length(wx, wy, wz):
    points          = np.stack([wx, wy, wz], axis=1)
    diffs           = np.diff(points, axis=0)
    segment_lengths = np.linalg.norm(diffs, axis=1)
    return np.sum(segment_lengths)

def point_to_segment_dist(
    px, py, pz, # point to check
    ax, ay, az, # segment start
    bx, by, bz  # segment end
):
    # get AB 
    ABx = bx - ax
    ABy = by - ay
    ABz = bz - az

    # get AP
    APx = px - ax
    APy = py - ay
    APz = pz - az

    # magnitude of AB
    ab2 = ABx * ABx + ABy * ABy + ABz * ABz

    # find t by projecting AP into AB
    t = 0
    if (ab2 > 0):
        t = (APx * ABx + APy * ABy + APz * ABz) / ab2
        if (t < 0):
            t = 0
        elif (t > 1):
            t = 1
        
    # compute Q = A + t AB
    qx = ax + t * ABx
    qy = ay + t * ABy
    qz = az + t * ABz

    dx = px - qx
    dy = py - qy
    dz = pz - qz

    return dx * dx + dy * dy + dz * dz

def is_point_inside_cylinder(
    primary_x, primary_y, primary_z,
    x, y, z,
    radius
):
    num_points = len(primary_x)
    r2 = radius * radius
    best = float("inf")
    if (num_points == 1):
        dx = x - primary_x[0]
        dy = y - primary_y[0]
        dz = z - primary_z[0]
        best = dx * dx + dy * dy + dz * dz
    else:
        for i in range(num_points - 1):
            d2 = point_to_segment_dist(
                x, y, z,
                primary_x[i], primary_y[i], primary_z[i],
                primary_x[i + 1], primary_y[i + 1], primary_z[i + 1]
            )
            if (d2 < best):
                best = d2
                if (best <= r2):
                    return True

    return (best <= r2)

def get_avg_direction(x, y, z):
    n = len(x)

    direction = []
    for i in range(n - 1):
        dx = x[i + 1] - x[i]
        dy = y[i + 1] - y[i]
        dz = z[i + 1] - z[i]
        direction.append((dx, dy, dz))
    
    return np.average(direction, axis=0)


def clean_up(df):
    # Get rid of entries with no matched WC2TPC track
    df = df[df["WC2TPCtrkID"] != -99999].copy()

    # Keep track of events tagged as beamline electrons with current selection
    is_beamline_electron_idx = []

    for i, row in df.iterrows():
        # Get WC2TPC positions
        wcX = row["WC2TPCLocationsX"]
        wcY = row["WC2TPCLocationsY"]
        wcZ = row["WC2TPCLocationsZ"]

        # Remove repeated entries at the end of wcX, wcY, wcZ
        repeat_count = 0
        num_points   = len(wcX)
        if num_points > 1:
            last_x, last_y, last_z = wcX[-1], wcY[-1], wcZ[-1]
            for j in range(num_points - 2, -1, -1):
                if wcX[j] == last_x and wcY[j] == last_y and wcZ[j] == last_z:
                    repeat_count += 1
                else:
                    break
        if repeat_count > 0:
            wcX = wcX[:num_points - repeat_count] + [wcX[-1]]
            wcY = wcY[:num_points - repeat_count] + [wcY[-1]]
            wcZ = wcZ[:num_points - repeat_count] + [wcZ[-1]]

        # Get information for WC2TPC track, before extending
        df.at[i, "WC2TPCTrackLength"] = compute_track_length(wcX, wcY, wcZ)
        df.at[i, "WC2TPCBeginX"] = wcX[0] if len(wcX) > 0 else -9999
        df.at[i, "WC2TPCBeginY"] = wcY[0] if len(wcY) > 0 else -9999
        df.at[i, "WC2TPCBeginZ"] = wcZ[0] if len(wcZ) > 0 else -9999
        df.at[i, "WC2TPCEndX"]   = wcX[-1] if len(wcX) > 0 else -9999
        df.at[i, "WC2TPCEndY"]   = wcY[-1] if len(wcY) > 0 else -9999
        df.at[i, "WC2TPCEndZ"]   = wcZ[-1] if len(wcZ) > 0 else -9999

        # Get direction towards end of track
        num_points = len(wcX)
        n_tail = min(10, num_points - 1)
        if n_tail > 0:
            avg_direction = get_avg_direction(wcX[-(n_tail+1):], wcY[-(n_tail+1):], wcZ[-(n_tail+1):])
        else:
            avg_direction = np.array([np.nan, np.nan, np.nan])
        
        # Extrapolate WC2TPC track to end
        if num_points > 0:
            last_point = np.array([wcX[-1], wcY[-1], wcZ[-1]])
            extrapolated = last_point + avg_direction * (constants.maxZ - last_point[2]) / avg_direction[2]
            wcX = np.append(wcX, extrapolated[0])
            wcY = np.append(wcY, extrapolated[1])
            wcZ = np.append(wcZ, extrapolated[2])

        trks_in_cylinder = 0
        small_trks_in_cylinder = 0
        num_reco_trks = len(row["recoTrkID"])

        trk_idx_len = []
        for i_trk in range(num_reco_trks):
            if (row["recoTrkID"][i_trk] == row["WC2TPCtrkID"]):
                df.at[i, "WC2TPCdEdx"] = np.mean(np.array(row["recoDEDX"][i_trk]))
                continue

            bx = row["recoBeginX"][i_trk]
            by = row["recoBeginY"][i_trk]
            bz = row["recoBeginZ"][i_trk]

            ex = row["recoEndX"][i_trk]
            ey = row["recoEndY"][i_trk]
            ez = row["recoEndZ"][i_trk]

            start_in_cylinder = is_point_inside_cylinder(
                wcX, wcY, wcZ,
                bx, by, bz,
                constants.CYLINDER_RADIUS
            )

            end_in_cylinder = is_point_inside_cylinder(
                wcX, wcY, wcZ,
                ex, ey, ez,
                constants.CYLINDER_RADIUS
            )

            trk_len = np.sqrt((bx - ex) * (bx - ex) + (by - ey) * (by - ey) + (bz - ez) * (bz - ez))
            if start_in_cylinder and end_in_cylinder:
                trks_in_cylinder += 1
                trk_idx_len.append((i_trk, trk_len))
                if trk_len < constants.CYLINDER_SMALL_TRACK:
                    small_trks_in_cylinder += 1

        # Check if we would tag it as beamline electron
        if (small_trks_in_cylinder > constants.ALLOWED_CYLINDER_SMALL_TRACKS):
            is_beamline_electron_idx.append(i)
        
        trk_idx_len.sort(key=lambda x: x[1], reverse=True)
        for j in range(constants.RECO_TRKS_SAVE):
            if j < len(trk_idx_len):
                bx = row["recoBeginX"][trk_idx_len[j][0]]
                by = row["recoBeginY"][trk_idx_len[j][0]]
                bz = row["recoBeginZ"][trk_idx_len[j][0]]

                ex = row["recoEndX"][trk_idx_len[j][0]]
                ey = row["recoEndY"][trk_idx_len[j][0]]
                ez = row["recoEndZ"][trk_idx_len[j][0]]

                df.at[i, f"recoTrkBeginX_{j}"] = bx
                df.at[i, f"recoTrkBeginY_{j}"] = by
                df.at[i, f"recoTrkBeginZ_{j}"] = bz

                df.at[i, f"recoTrkEndX_{j}"] = ex
                df.at[i, f"recoTrkEndY_{j}"] = ey
                df.at[i, f"recoTrkEndZ_{j}"] = ez

                df.at[i, f"recoTrkLen_{j}"] = np.sqrt((bx - ex) * (bx - ex) + (by - ey) * (by - ey) + (bz - ez) * (bz - ez))
                if len(row["recoDEDX"][trk_idx_len[j][0]]) == 0:
                    df.at[i, f"recoTrkdEdx_{j}"] = -9999
                else:
                    df.at[i, f"recoTrkdEdx_{j}"] = np.mean(np.array(row["recoDEDX"][trk_idx_len[j][0]]))

            else:
                df.at[i, f"recoTrkBeginX_{j}"] = -9999
                df.at[i, f"recoTrkBeginY_{j}"] = -9999
                df.at[i, f"recoTrkBeginZ_{j}"] = -9999

                df.at[i, f"recoTrkEndX_{j}"] = -9999
                df.at[i, f"recoTrkEndY_{j}"] = -9999
                df.at[i, f"recoTrkEndZ_{j}"] = -9999

                df.at[i, f"recoTrkLen_{j}"]  = -9999
                df.at[i, f"recoTrkdEdx_{j}"] = -9999
        df.at[i, "numRecoTrksInCylinder"] = trks_in_cylinder
    
    # Drop events tagged as beamline electrons
    # df = df.drop(index=is_beamline_electron_idx)

    # Drop variables not needed for training
    df = df.drop(columns=[var for var in variables.bdt_import_variables if var in df.columns])

    return df


if (__name__ == "__main__"):
    filepath = "/exp/lariat/app/users/epelaez/files"
    filename = "RecoNNAllEval_histo.root:RecoNNAllEval"
    treename = "RecoNNAllEvalTree"

    with uproot.open(f"{filepath}/{filename}") as file:
        ttree         = file[treename]
        total_entries = ttree.num_entries
        
        print(f"Opened tree {treename} successfully.")
        print(f"Total entries: {total_entries}")
        print(f"Entries to read: {min(constants.EVENTS_TO_READ, total_entries)}")

        # Create dataframe with selected variables
        keep_vars = variables.bdt_import_variables + variables.bdt_keep_variables
        arrays    = ttree.arrays(keep_vars, library="np", entry_stop=min(constants.EVENTS_TO_READ, total_entries))
        df        = pd.DataFrame(arrays)
        
        # Clean up dataframe to get rid of events we reject
        # and put data in format suitable for BDT training
        df = clean_up(df)

        print(f"Entries after cleanup: {df.shape[0]}")

        # Save the DataFrame to a pickle file
        df.to_pickle("files/train_initial_data.pkl")