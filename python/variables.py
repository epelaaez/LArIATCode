# list of variables to import
bdt_import_variables = [
    "WC2TPCtrkID",      # WC2TPC track ID
    "WC2TPCLocationsX", # x position of the WC2TPC track
    "WC2TPCLocationsY", # y position of the WC2TPC track
    "WC2TPCLocationsZ", # z position of the WC2TPC track
    "recoTrkID",        # reconstructed track ID
    "recoBeginX",       # x position of reconstructed tracks start
    "recoBeginY",       # y position of reconstructed tracks start
    "recoBeginZ",       # z position of reconstructed tracks start
    "recoEndX",         # x position of reconstructed tracks end
    "recoEndY",         # y position of reconstructed tracks end
    "recoEndZ",         # z position of reconstructed tracks end
    "recoDEDX",         # dE/dx of reconstructed tracks
    "recoResR"          # residual range of reconstructed tracks
]

# other variables we want to keep for bookkeping and truth-tagging results
bdt_keep_variables = [
    "backgroundType"
]