#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TEllipse.h"
#include "TLine.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TPRegexp.h"
#include "TPaletteAxis.h"

#include <vector>
#include <map>
#include <fstream>
#include <algorithm>

#include "Helpers.cc"

// Hit-level event displays (wire vs drift coordinate) for events saved in the SelectedEvents
// tree by RecoClassify3Cat.cc / DataClassify.cc. Uses the same coordinates as the hit cluster cut:
// W = channel * 0.4 cm (induction), (channel - 240) * 0.4 cm (collection); X = driftT * 0.128 us * 0.14732 cm/us.
//
// Hit views:
//   Reco column:  raw-display-style charge map: each hit is drawn as a Gaussian pulse (width hit_rms) on its wire,
//                 colored by amplitude relative to the median primary-track hit in that plane (blue = no charge -> red =
//                 CHARGE_MAX x median). Reco tracks are lines connecting their hits: black = WC-matched primary up to the
//                 vertex, purple = primary past the vertex when the primary was stitched (that segment is PID'd as a
//                 secondary), white = other tracks. Tracks with an endpoint within VERTEX_RADIUS of the vertex are labeled
//                 with their PID tag.
//   Truth column: (MC only) colored by the true particle that deposited the hit.
// The dashed circle is the VERTEX_RADIUS sphere around the vertex (a sphere projects to a circle of the same radius).
// Right column: collection-plane dE/dx vs residual range (primary + tracks near the vertex, with templates),
//               and dE/dx along the primary vs z with the vertex marked.
//
// Category: 0 = abs 0p, 1 = abs Np, 2 = scatter. Background (MC only): true category, -99 = any.

namespace {
    const int kMaxHitsED  = 20000;
    const int kMaxG4ED    = 5000;

    const double SAMPLING  = 0.128;
    const double DRIFT_VEL = 0.14732;

    const double WIRE_PITCH = 0.4;
    const double DRIFT_BIN  = 0.05; // cm, charge-map binning along the drift direction

    const int    kNChargeColors = 60;
    const double CHARGE_MAX     = 3.0; // in units of the median primary hit charge

    // Blue (no charge) -> cyan -> green -> yellow -> red (high charge); returns the first palette color index
    int chargePalette() {
        static int first = -1;
        if (first < 0) {
            double stops[] = {0.00, 0.25, 0.50, 0.75, 1.00};
            double red[]   = {0.00, 0.00, 0.20, 1.00, 0.85};
            double green[] = {0.10, 0.75, 0.80, 0.85, 0.00};
            double blue[]  = {0.85, 0.95, 0.20, 0.00, 0.00};
            first = TColor::CreateGradientColorTable(5, stops, red, green, blue, kNChargeColors);
        }
        return first;
    }

    const std::map<int, TString> kBackgroundNames = {
        {0, "Abs 0p"}, {1, "Abs Np"}, {2, "Primary muon"}, {3, "Primary electron"}, {4, "Other primary"},
        {5, "Outside reduced volume"}, {6, "Inelastic scattering"}, {7, "Charge exchange"}, {8, "Double charge exchange"},
        {9, "Capture at rest"}, {10, "Decay"}, {11, "Other"}, {12, "Elastic scattering"}, {13, "Scattering 0p"}, {14, "Scattering Np"}
    };
    const char* kCategoryNames[]   = {"abs0p", "absNp", "scatter"};
    const char* kCategoryLabels[]  = {"abs 0p", "abs Np", "scatter"};
    const char* kStitchedLabels[]  = {"not stitched", "pi(chi2)", "p(chi2)", "pi(dEdx)", "p(dEdx)"};

    double hitW(int channel) { return (channel < 240) ? channel * 0.4 : (channel - 240) * 0.4; }
    double hitXcm(float driftT) { return driftT * SAMPLING * DRIFT_VEL; }

    int trackColor(int idx) {
        static const int cols[] = {
            TColor::GetColor(0, 114, 178), TColor::GetColor(213, 94, 0), TColor::GetColor(0, 158, 115),
            TColor::GetColor(230, 159, 0), TColor::GetColor(86, 180, 233), TColor::GetColor(160, 40, 40),
            TColor::GetColor(204, 121, 167), TColor::GetColor(120, 120, 0)
        };
        return cols[idx % 8];
    }
    const int kStitchColor = TColor::GetColor(140, 60, 200);

    int truthColor(int pdgCode, bool isPrimary) {
        if (isPrimary)                                return kBlack;
        if (pdgCode == 2212)                          return TColor::GetColor(213, 94, 0);
        if (std::abs(pdgCode) == 211)                 return TColor::GetColor(0, 114, 178);
        if (std::abs(pdgCode) == 13)                  return TColor::GetColor(0, 158, 115);
        if (std::abs(pdgCode) == 11 || pdgCode == 22) return TColor::GetColor(200, 185, 0);
        return TColor::GetColor(204, 121, 167);
    }

    double median(std::vector<double> v) {
        if (v.empty()) return 0;
        std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
        return v[v.size() / 2];
    }
}

void EventDisplay(
    TString Sample   = "MC",
    int Category     = 1,
    int NEvents      = 50,
    int Background   = -99,
    int Seed         = 1,
    TString OutDir   = "/exp/lariat/app/users/epelaez/analysis_matt/figs/SecondaryStudies/Displays/",
    double HalfW     = 25.0
) {
    gStyle->SetOptStat(0);
    chargePalette();
    TH1::AddDirectory(false);
    bool isMC = (Sample == "MC");

    TGraph* gProton = new TGraph(); initializeProtonPoints(gProton);
    TGraph* gPion   = new TGraph(); initializePionPoints(gPion);
    gProton->SetLineColor(kRed); gProton->SetLineWidth(2);
    gPion->SetLineColor(kGray + 2); gPion->SetLineWidth(2); gPion->SetLineStyle(2);

    TString selPath = isMC ? "/exp/lariat/app/users/epelaez/analysis_matt/files/SecondaryStudies_MC.root"
                           : "/exp/lariat/app/users/epelaez/analysis_matt/files/SecondaryStudies_Data.root";
    std::unique_ptr<TFile> fSel(TFile::Open(selPath, "READ"));
    TTree* tSel = (TTree*) fSel->Get("SelectedEvents");

    int category, run, subrun, event, background, stitched = 0; Long64_t entry;
    double vertexX, vertexY, vertexZ, energyAtVertex;
    std::string* file = nullptr; std::string* tracks = nullptr;
    tSel->SetBranchAddress("category", &category);
    tSel->SetBranchAddress("file", &file);
    tSel->SetBranchAddress("entry", &entry);
    tSel->SetBranchAddress("run", &run);
    tSel->SetBranchAddress("subrun", &subrun);
    tSel->SetBranchAddress("event", &event);
    tSel->SetBranchAddress("background", &background);
    tSel->SetBranchAddress("vertexX", &vertexX);
    tSel->SetBranchAddress("vertexY", &vertexY);
    tSel->SetBranchAddress("vertexZ", &vertexZ);
    tSel->SetBranchAddress("energyAtVertex", &energyAtVertex);
    tSel->SetBranchAddress("tracks", &tracks);
    if (tSel->GetBranch("stitched")) tSel->SetBranchAddress("stitched", &stitched);

    // Candidate events, then a reproducible random subset
    std::vector<Long64_t> candidates;
    for (Long64_t i = 0; i < tSel->GetEntries(); ++i) {
        tSel->GetEntry(i);
        if (category != Category) continue;
        if (isMC && Background != -99 && background != Background) continue;
        candidates.push_back(i);
    }
    TRandom3 rng(Seed);
    for (int i = (int) candidates.size() - 1; i > 0; --i) std::swap(candidates[i], candidates[rng.Integer(i + 1)]);
    if ((int) candidates.size() > NEvents) candidates.resize(NEvents);
    std::cout << "Selected " << candidates.size() << " events to display" << std::endl;

    TString tag = Form("%s_%s", Sample.Data(), kCategoryNames[Category]);
    if (isMC && Background != -99) tag += Form("_bg%d", Background);
    TString dir = OutDir + tag + "/";
    gSystem->mkdir(dir, true);
    std::ofstream log((dir + "events.txt").Data());

    // Anatree branches
    static int   nhits, hit_plane[kMaxHitsED], hit_channel[kMaxHitsED], hit_trkid[kMaxHitsED], hit_g4id[kMaxHitsED];
    static float hit_driftT[kMaxHitsED], hit_charge[kMaxHitsED], hit_rms[kMaxHitsED], hit_x[kMaxHitsED], hit_y[kMaxHitsED], hit_z[kMaxHitsED];
    static int   ntracks_reco, trkWCtoTPCMatch[kMaxTrack];
    static float trkvtxx[kMaxTrack], trkvtxy[kMaxTrack], trkvtxz[kMaxTrack], trkendx[kMaxTrack], trkendy[kMaxTrack], trkendz[kMaxTrack];
    static int   ntrkcalopts[kMaxTrack][2];
    static float trkdedx[kMaxTrack][2][kMaxTrackHits], trkrr[kMaxTrack][2][kMaxTrackHits], trkxyz[kMaxTrack][2][kMaxTrackHits][3];
    static int   geant_list_size, pdg[kMaxG4ED], TrackId[kMaxG4ED], Mother[kMaxG4ED], process_primary[kMaxG4ED];
    static float Eng[kMaxG4ED], Mass[kMaxG4ED], Px[kMaxG4ED], Py[kMaxG4ED], Pz[kMaxG4ED];

    TPRegexp trackRegex("trk (\\d+): .*tag=(\\S+)");

    for (size_t iCand = 0; iCand < candidates.size(); ++iCand) {
        tSel->GetEntry(candidates[iCand]);

        std::unique_ptr<TFile> fIn(TFile::Open(file->c_str(), "READ"));
        TTree* t = (TTree*) fIn->Get("anatree/anatree");
        t->SetBranchStatus("*", 0);
        auto bind = [&](const char* name, void* addr) { t->SetBranchStatus(name, 1); t->SetBranchAddress(name, addr); };
        bind("nhits", &nhits);             bind("hit_plane", hit_plane);   bind("hit_channel", hit_channel);
        bind("hit_trkid", hit_trkid);      bind("hit_driftT", hit_driftT); bind("hit_charge", hit_charge); bind("hit_rms", hit_rms);
        bind("hit_x", hit_x);              bind("hit_y", hit_y);           bind("hit_z", hit_z);
        bind("ntracks_reco", &ntracks_reco); bind("trkWCtoTPCMatch", trkWCtoTPCMatch);
        bind("trkvtxx", trkvtxx); bind("trkvtxy", trkvtxy); bind("trkvtxz", trkvtxz);
        bind("trkendx", trkendx); bind("trkendy", trkendy); bind("trkendz", trkendz);
        bind("ntrkcalopts", ntrkcalopts); bind("trkdedx", trkdedx); bind("trkrr", trkrr); bind("trkxyz", trkxyz);
        if (isMC) {
            bind("hit_g4id", hit_g4id); bind("geant_list_size", &geant_list_size); bind("pdg", pdg);
            bind("TrackId", TrackId);   bind("Mother", Mother);   bind("process_primary", process_primary);
            bind("Eng", Eng); bind("Mass", Mass); bind("Px", Px); bind("Py", Py); bind("Pz", Pz);
        }
        t->GetEntry(entry);
        int nHits   = std::min(nhits, kMaxHitsED);
        int nTracks = std::min(ntracks_reco, kMaxTrack);

        int primaryTrk = -1;
        for (int i = 0; i < nTracks; ++i) if (trkWCtoTPCMatch[i]) { primaryTrk = i; break; }

        // Tracks with an endpoint within VERTEX_RADIUS of the vertex (same criterion as the selection), with PID tags
        std::map<int, TString> nearTags;
        TString trackText = tracks->c_str();
        TObjArray* trackLines = trackText.Tokenize("\n");
        for (int i = 0; i < trackLines->GetEntries(); ++i) {
            TString line = ((TObjString*) trackLines->At(i))->GetString();
            TObjArray* m = trackRegex.MatchS(line);
            if (m->GetEntries() == 3) nearTags[((TObjString*) m->At(1))->GetString().Atoi()] = ((TObjString*) m->At(2))->GetString();
            delete m;
        }
        // 3D distance of the closest track endpoint to the vertex
        auto vertexDistance = [&](int trk) {
            double d1 = std::hypot(std::hypot(trkvtxx[trk] - vertexX, trkvtxy[trk] - vertexY), trkvtxz[trk] - vertexZ);
            double d2 = std::hypot(std::hypot(trkendx[trk] - vertexX, trkendy[trk] - vertexY), trkendz[trk] - vertexZ);
            return std::min(d1, d2);
        };
        auto isNear = [&](int trk) { return vertexDistance(trk) < VERTEX_RADIUS; };

        // Primary direction (pointing downstream) to split primary hits before / after the vertex
        double pdx = 0, pdy = 0, pdz = 1;
        if (primaryTrk >= 0) {
            pdx = trkendx[primaryTrk] - trkvtxx[primaryTrk]; pdy = trkendy[primaryTrk] - trkvtxy[primaryTrk]; pdz = trkendz[primaryTrk] - trkvtxz[primaryTrk];
            if (pdz < 0) { pdx = -pdx; pdy = -pdy; pdz = -pdz; }
        }
        auto pastVertex = [&](double x, double y, double z) { return ((x - vertexX) * pdx + (y - vertexY) * pdy + (z - vertexZ) * pdz) > 0.3; };

        // Vertex projection with the wire geometry (+/-60 deg planes), fitted from track hits (residual rms 0.05 cm):
        //   W = -/+ 0.5 y + 0.866 z + c_plane ;  X = x + offset (MC drift-time offset of 0.96 cm)
        const double wOffset[2] = {8.420, 9.264};
        const double xOffset    = isMC ? 0.96 : 0.0;
        double vtxW[2], vtxX[2], best[2] = {1e9, 1e9};
        for (int pl = 0; pl < 2; ++pl) {
            vtxW[pl] = (pl == 0 ? -0.5 : 0.5) * vertexY + 0.8660 * vertexZ + wOffset[pl];
            vtxX[pl] = vertexX + xOffset;
        }
        // Hit pulse amplitude (Gaussian peak with integral = charge), and its median over primary-track hits per plane
        auto hitAmp = [&](int i) { return hit_charge[i] / (std::sqrt(2 * M_PI) * std::max(1.0f, hit_rms[i])); };
        std::vector<double> primaryAmps[2];
        for (int i = 0; i < nHits; ++i) {
            int pl = hit_plane[i];
            if (pl < 0 || pl > 1 || hit_trkid[i] != primaryTrk) continue;
            primaryAmps[pl].push_back(hitAmp(i));
            double d = std::hypot(std::hypot(hit_x[i] - vertexX, hit_y[i] - vertexY), hit_z[i] - vertexZ);
            if (d < best[pl]) best[pl] = d;
        }
        double refAmp[2];
        for (int pl = 0; pl < 2; ++pl) refAmp[pl] = primaryAmps[pl].empty() ? 1.0 : median(primaryAmps[pl]);

        // MC truth
        std::map<int, int> g4Index;
        int primaryG4 = -1;
        if (isMC) {
            for (int i = 0; i < std::min(geant_list_size, kMaxG4ED); ++i) {
                g4Index[TrackId[i]] = i;
                if (primaryG4 < 0 && process_primary[i] == 1 && Mother[i] == 0) primaryG4 = TrackId[i];
            }
        }

        // Canvas layout: hit views (1 or 2 columns) + calorimetry column; text at the bottom
        int nViewCols = isMC ? 2 : 1;
        double viewFrac = isMC ? 0.70 : 0.62;
        int canvasW = isMC ? 2600 : 2000, canvasH = 2000;
        TCanvas* c = new TCanvas("c", "", canvasW, canvasH);
        const double textH = 0.22, rowH = (1 - textH) / 2;
        std::vector<TObject*> keep;

        c->cd();
        TPad* text = new TPad("text", "", 0.0, 0.0, 1.0, textH); text->Draw(); keep.push_back(text);

        for (int pl = 0; pl < 2; ++pl) {
            for (int col = 0; col < nViewCols; ++col) {
                c->cd();
                double x1 = viewFrac * col / nViewCols, x2 = viewFrac * (col + 1) / nViewCols;
                double y1 = textH + (1 - pl) * rowH, y2 = y1 + rowH;
                TPad* p = new TPad(Form("p%d%d", pl, col), "", x1, y1, x2, y2);
                double lm = 0.09, rm = (col == 0) ? 0.11 : 0.02, tm = 0.09, bm = 0.12;
                p->SetLeftMargin(lm); p->SetRightMargin(rm); p->SetTopMargin(tm); p->SetBottomMargin(bm);
                p->Draw(); p->cd(); keep.push_back(p);

                // Equal scale in both directions so the vertex-radius circle is round
                double padWpx = (x2 - x1) * canvasW * (1 - lm - rm), padHpx = (y2 - y1) * canvasH * (1 - tm - bm);
                double halfX  = HalfW * padHpx / padWpx;
                double cw = (vtxW[pl] > -900) ? vtxW[pl] : 45.0, cx = (vtxX[pl] > -900) ? vtxX[pl] : 23.5;

                TH2F* frame = new TH2F(Form("frame%d%d", pl, col),
                    Form("%s plane, %s;Wire coordinate [cm];Drift coordinate [cm]",
                         pl == 0 ? "Induction" : "Collection", col == 0 ? "color = charge (1 = median primary-track hit)" : "true particles"),
                    10, cw - HalfW, cw + HalfW, 10, cx - halfX, cx + halfX);
                frame->GetXaxis()->SetTitleSize(0.05); frame->GetYaxis()->SetTitleSize(0.05);
                frame->GetXaxis()->SetLabelSize(0.045); frame->GetYaxis()->SetLabelSize(0.045);
                frame->GetYaxis()->SetTitleOffset(0.8);
                if (col == 0) {
                    // "COL SAME" takes its z range from the first histogram in the pad, so the frame carries the charge scale
                    p->SetFrameFillColor(chargePalette());
                    frame->SetMinimum(0); frame->SetMaximum(CHARGE_MAX); frame->SetContour(kNChargeColors);
                }
                frame->Draw();
                keep.push_back(frame);

                // Vertex radius
                if (vtxW[pl] > -900) {
                    TEllipse* circ = new TEllipse(vtxW[pl], vtxX[pl], VERTEX_RADIUS, VERTEX_RADIUS);
                    circ->SetFillStyle(0); circ->SetLineStyle(2); circ->SetLineColor(kRed - 7); circ->SetLineWidth(2);
                    circ->Draw();
                    keep.push_back(circ);
                }

                if (col == 0) {
                    // Charge map, raw-display style: each hit is a Gaussian pulse (width hit_rms) on its wire, with
                    // amplitude charge / (sqrt(2 pi) rms); no deposit = palette minimum (frame filled in that color)
                    int    w0 = (int) std::floor((cw - HalfW) / WIRE_PITCH) - 1, nW = (int) std::ceil(2 * HalfW / WIRE_PITCH) + 3;
                    int    nX = (int) std::ceil(2 * halfX / DRIFT_BIN);
                    TH2F* qmap = new TH2F(Form("qmap%d", pl), "", nW, (w0 - 0.5) * WIRE_PITCH, (w0 + nW - 0.5) * WIRE_PITCH,
                                          nX, cx - halfX, cx + halfX);
                    const double tickCm = SAMPLING * DRIFT_VEL;
                    for (int i = 0; i < nHits; ++i) {
                        if (hit_plane[i] != pl) continue;
                        double w = hitW(hit_channel[i]), x = hitXcm(hit_driftT[i]);
                        if (std::abs(w - cw) > HalfW + 1 || std::abs(x - cx) > halfX + 1) continue;
                        double sigma = std::max(1.0, (double) hit_rms[i]) * tickCm, a = hitAmp(i) / refAmp[pl];
                        int bw = qmap->GetXaxis()->FindBin(w);
                        for (int bx = qmap->GetYaxis()->FindBin(x - 3 * sigma); bx <= qmap->GetYaxis()->FindBin(x + 3 * sigma); ++bx) {
                            double d = (qmap->GetYaxis()->GetBinCenter(bx) - x) / sigma;
                            qmap->SetBinContent(bw, bx, qmap->GetBinContent(bw, bx) + a * std::exp(-0.5 * d * d));
                        }
                    }
                    for (int bw = 1; bw <= nW; ++bw) for (int bx = 1; bx <= nX; ++bx)
                        if (qmap->GetBinContent(bw, bx) > 0.999 * CHARGE_MAX) qmap->SetBinContent(bw, bx, 0.999 * CHARGE_MAX);
                    qmap->SetMinimum(0); qmap->SetMaximum(CHARGE_MAX); qmap->SetContour(kNChargeColors);
                    qmap->GetZaxis()->SetLabelSize(0.04);
                    qmap->Draw("COL SAME"); keep.push_back(qmap);
                    double perNDC = 2 * HalfW / (1 - lm - rm);
                    TPaletteAxis* pal = new TPaletteAxis(cw + HalfW + 0.012 * perNDC, cx - halfX, cw + HalfW + 0.04 * perNDC, cx + halfX, qmap);
                    pal->Draw(); keep.push_back(pal);
                    gPad->RedrawAxis();

                    // Reco tracks: a line connecting each track's hits, ordered along the track
                    std::map<int, std::vector<std::pair<double, int>>> trackHits; // trk -> (position along track, hit)
                    for (int i = 0; i < nHits; ++i) {
                        int trk = hit_trkid[i];
                        if (hit_plane[i] != pl || trk < 0 || trk >= nTracks) continue;
                        double s = (hit_x[i] - trkvtxx[trk]) * (trkendx[trk] - trkvtxx[trk]) + (hit_y[i] - trkvtxy[trk]) * (trkendy[trk] - trkvtxy[trk])
                                 + (hit_z[i] - trkvtxz[trk]) * (trkendz[trk] - trkvtxz[trk]);
                        trackHits[trk].push_back({s, i});
                    }
                    for (auto& [trk, hits] : trackHits) {
                        std::sort(hits.begin(), hits.end());
                        // Split into segments where the before/after-vertex state changes (stitched primary)
                        std::vector<TGraph*> segs; int segState = -1;
                        for (size_t k = 0; k < hits.size(); ++k) {
                            int i = hits[k].second;
                            int state = (trk == primaryTrk && stitched > 0 && pastVertex(hit_x[i], hit_y[i], hit_z[i]));
                            if (state != segState) {
                                TGraph* seg = new TGraph();
                                if (!segs.empty()) { int j = hits[k - 1].second; seg->SetPoint(0, hitW(hit_channel[j]), hitXcm(hit_driftT[j])); }
                                segState = state;
                                seg->SetLineColor(trk == primaryTrk ? (state ? kStitchColor : kBlack) : kWhite);
                                seg->SetLineWidth(3);
                                segs.push_back(seg);
                            }
                            segs.back()->SetPoint(segs.back()->GetN(), hitW(hit_channel[i]), hitXcm(hit_driftT[i]));
                        }
                        for (TGraph* seg : segs) {
                            if (seg->GetN() >= 2) seg->Draw("L SAME");
                            keep.push_back(seg);
                        }
                    }
                } else {
                    // Truth: hits colored by the true particle, one graph per color
                    std::map<int, TGraph*> graphs;
                    for (int i = 0; i < nHits; ++i) {
                        if (hit_plane[i] != pl) continue;
                        double w = hitW(hit_channel[i]), x = hitXcm(hit_driftT[i]);
                        if (std::abs(w - cw) > HalfW || std::abs(x - cx) > halfX) continue;
                        auto it = g4Index.find(hit_g4id[i]);
                        int color = (it == g4Index.end()) ? kGray + 1 : truthColor(pdg[it->second], hit_g4id[i] == primaryG4);
                        if (!graphs.count(color)) graphs[color] = new TGraph();
                        graphs[color]->SetPoint(graphs[color]->GetN(), w, x);
                    }
                    for (auto& [color, g] : graphs) {
                        g->SetMarkerStyle(20); g->SetMarkerColor(color); g->SetMarkerSize(1.0);
                        g->Draw("P SAME");
                        keep.push_back(g);
                    }
                }

                // Vertex marker (explained in the legend)
                TMarker* vtxMarker = new TMarker(vtxW[pl], vtxX[pl], 34);
                vtxMarker->SetMarkerColor(kRed); vtxMarker->SetMarkerSize(1.4);
                if (vtxW[pl] > -900) vtxMarker->Draw();
                keep.push_back(vtxMarker);
                TLine* radiusLine = new TLine(); radiusLine->SetLineStyle(2); radiusLine->SetLineColor(kRed - 7); radiusLine->SetLineWidth(2);
                keep.push_back(radiusLine);

                if (col == 0) {
                    // Label each reco track at its hit farthest from the vertex (within the window)
                    std::map<int, std::pair<double, double>> labelPos;
                    for (int i = 0; i < nHits; ++i) {
                        if (hit_plane[i] != pl || hit_trkid[i] < 0 || hit_trkid[i] == primaryTrk) continue;
                        double w = hitW(hit_channel[i]), x = hitXcm(hit_driftT[i]);
                        if (std::abs(w - cw) > 0.95 * HalfW || std::abs(x - cx) > 0.92 * halfX) continue;
                        double dist = std::hypot(w - cw, x - cx);
                        if (!labelPos.count(hit_trkid[i]) || dist > std::hypot(labelPos[hit_trkid[i]].first - cw, labelPos[hit_trkid[i]].second - cx))
                            labelPos[hit_trkid[i]] = {w, x};
                    }
                    // Compact labels: "<trk>:<PID>" for tracks within the vertex radius, just "<trk>" otherwise
                    // (full PID info and vertex distances are in the text panel)
                    for (auto& [trk, pos] : labelPos) {
                        bool near = isNear(trk);
                        TString pid = nearTags.count(trk) ? nearTags[trk] : TString("?");
                        pid = pid.BeginsWith("pi") ? "#pi" : (pid.BeginsWith("p") ? "p" : pid);
                        TString lab = near ? Form("%d:%s", trk, pid.Data()) : Form("%d", trk);
                        TLatex* l = new TLatex(pos.first, pos.second, lab);
                        l->SetTextAlign(pos.first > cw + 0.2 * HalfW ? 32 : 12);
                        l->SetTextSize(near ? 0.04 : 0.03); l->SetTextColor(near ? kWhite : kGray); l->SetTextFont(near ? 62 : 42);
                        l->Draw(); keep.push_back(l);
                    }

                    // Transparent legend, white text on the charge map
                    TLegend* leg = new TLegend(0.10, stitched > 0 ? 0.62 : 0.66, 0.36, 0.905);
                    leg->SetTextSize(0.03); leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextColor(kWhite);
                    TLine* lp = new TLine(); lp->SetLineColor(kBlack); lp->SetLineWidth(3);
                    leg->AddEntry(lp, "Primary (to vertex)", "l"); keep.push_back(lp);
                    if (stitched > 0) {
                        TLine* ls = new TLine(); ls->SetLineColor(kStitchColor); ls->SetLineWidth(3);
                        leg->AddEntry(ls, "Primary past vertex (stitched)", "l"); keep.push_back(ls);
                    }
                    TLine* lo = new TLine(); lo->SetLineColor(kWhite); lo->SetLineWidth(3);
                    leg->AddEntry(lo, "Other reco tracks", "l"); keep.push_back(lo);
                    leg->AddEntry(vtxMarker, "Vertex", "p");
                    leg->AddEntry(radiusLine, Form("%.0f cm vertex radius", VERTEX_RADIUS), "l");
                    leg->AddEntry((TObject*) nullptr, "n:#pi / n:p = track n + PID", "");
                    leg->Draw(); keep.push_back(leg);
                } else {
                    TLegend* leg = new TLegend(0.86, 0.66, 0.975, 0.905);
                    leg->SetTextSize(0.032); leg->SetFillStyle(0); leg->SetBorderSize(0);
                    std::vector<std::pair<TString, int>> entries = {
                        {"Primary", kBlack}, {"p", truthColor(2212, false)}, {"#pi^{#pm}", truthColor(211, false)},
                        {"#mu", truthColor(13, false)}, {"e/#gamma", truthColor(11, false)}, {"Other", truthColor(0, false)}
                    };
                    for (auto& [lab, colr] : entries) {
                        TMarker* mk = new TMarker(0, 0, 20); mk->SetMarkerColor(colr); mk->SetMarkerSize(1.2);
                        leg->AddEntry(mk, lab, "p"); keep.push_back(mk);
                    }
                    leg->Draw(); keep.push_back(leg);
                }
            }
        }

        // Calorimetry column (collection plane)
        {
            // Top: dE/dx vs residual range for the primary and tracks near the vertex
            c->cd();
            TPad* p = new TPad("pcal1", "", viewFrac, textH + rowH, 1.0, 1.0);
            p->SetLeftMargin(0.14); p->SetRightMargin(0.03); p->SetTopMargin(0.09); p->SetBottomMargin(0.12);
            p->Draw(); p->cd(); keep.push_back(p);
            TH2F* frame = new TH2F("fcal1", "Collection dE/dx vs residual range;Residual range [cm];dE/dx [MeV/cm]", 10, 0, 30, 10, 0, 30);
            frame->GetXaxis()->SetTitleSize(0.05); frame->GetYaxis()->SetTitleSize(0.05);
            frame->GetXaxis()->SetLabelSize(0.045); frame->GetYaxis()->SetLabelSize(0.045);
            frame->Draw(); keep.push_back(frame);
            gProton->Draw("L SAME"); gPion->Draw("L SAME");
            TLegend* leg = new TLegend(0.45, 0.62, 0.97, 0.90);
            leg->SetTextSize(0.04); leg->SetBorderSize(0); leg->SetFillStyle(0);
            leg->AddEntry(gProton, "Proton template", "l");
            leg->AddEntry(gPion, "Pion template", "l");
            std::vector<int> drawTracks;
            if (primaryTrk >= 0) drawTracks.push_back(primaryTrk);
            for (int trk = 0; trk < nTracks; ++trk) if (trk != primaryTrk && isNear(trk)) drawTracks.push_back(trk);
            for (int trk : drawTracks) {
                int n = std::min(ntrkcalopts[trk][1], kMaxTrackHits);
                TGraph* g = new TGraph();
                for (int iPt = 0; iPt < n; ++iPt) if (trkrr[trk][1][iPt] < 30) g->SetPoint(g->GetN(), trkrr[trk][1][iPt], std::min(trkdedx[trk][1][iPt], 29.5f));
                if (g->GetN() == 0) { delete g; continue; }
                int colr = (trk == primaryTrk) ? kBlack : trackColor(trk);
                g->SetMarkerStyle(trk == primaryTrk ? 24 : 20); g->SetMarkerSize(1.1); g->SetMarkerColor(colr);
                g->Draw("P SAME"); keep.push_back(g);
                leg->AddEntry(g, trk == primaryTrk ? TString("Primary (full track)") : Form("trk %d", trk), "p");
            }
            leg->Draw(); keep.push_back(leg);

            // Bottom: dE/dx along the primary vs z, with the vertex marked
            c->cd();
            TPad* q = new TPad("pcal2", "", viewFrac, textH, 1.0, textH + rowH);
            q->SetLeftMargin(0.14); q->SetRightMargin(0.03); q->SetTopMargin(0.09); q->SetBottomMargin(0.12);
            q->Draw(); q->cd(); keep.push_back(q);
            TH2F* frame2 = new TH2F("fcal2", "Collection dE/dx along the primary;z [cm];dE/dx [MeV/cm]", 10, vertexZ - 30, vertexZ + 20, 10, 0, 30);
            frame2->GetXaxis()->SetTitleSize(0.05); frame2->GetYaxis()->SetTitleSize(0.05);
            frame2->GetXaxis()->SetLabelSize(0.045); frame2->GetYaxis()->SetLabelSize(0.045);
            frame2->Draw(); keep.push_back(frame2);
            if (primaryTrk >= 0) {
                int n = std::min(ntrkcalopts[primaryTrk][1], kMaxTrackHits);
                TGraph* gBefore = new TGraph(); TGraph* gAfter = new TGraph();
                for (int iPt = 0; iPt < n; ++iPt) {
                    float* pos = trkxyz[primaryTrk][1][iPt];
                    TGraph* g = pastVertex(pos[0], pos[1], pos[2]) ? gAfter : gBefore;
                    g->SetPoint(g->GetN(), pos[2], std::min(trkdedx[primaryTrk][1][iPt], 29.5f));
                }
                gBefore->SetMarkerStyle(20); gBefore->SetMarkerColor(kBlack); gBefore->SetMarkerSize(1.1);
                gAfter->SetMarkerStyle(20); gAfter->SetMarkerColor(kStitchColor); gAfter->SetMarkerSize(1.1);
                if (gBefore->GetN()) gBefore->Draw("P SAME");
                if (gAfter->GetN()) gAfter->Draw("P SAME");
                keep.push_back(gBefore); keep.push_back(gAfter);
            }
            TLine* vl = new TLine(vertexZ, 0, vertexZ, 30);
            vl->SetLineColor(kRed); vl->SetLineWidth(3); vl->SetLineStyle(2);
            vl->Draw(); keep.push_back(vl);
            TLatex* vt = new TLatex(vertexZ + 0.5, 27, "vertex");
            vt->SetTextColor(kRed); vt->SetTextSize(0.05); vt->SetTextFont(62);
            vt->Draw(); keep.push_back(vt);
        }

        // Text panel
        text->cd();
        std::vector<TString> lines;
        lines.push_back(Form("%s  run %d / subrun %d / event %d   |   selected as %s   |   KE at vertex %.0f MeV   |   vertex (%.1f, %.1f, %.1f) cm   |   stitched primary: %s",
            Sample.Data(), run, subrun, event, kCategoryLabels[Category], energyAtVertex, vertexX, vertexY, vertexZ, kStitchedLabels[stitched]));
        if (isMC) {
            TString truth = Form("True category: %s.  True daughters of primary (KE > 10 MeV):", kBackgroundNames.count(background) ? kBackgroundNames.at(background).Data() : "?");
            for (int i = 0; i < std::min(geant_list_size, kMaxG4ED); ++i) {
                if (Mother[i] != primaryG4) continue;
                double ke = 1000 * (Eng[i] - Mass[i]);
                if (ke < 10 || pdg[i] == 2112 || pdg[i] > 100000) continue;
                double pp = std::sqrt(Px[i] * Px[i] + Py[i] * Py[i] + Pz[i] * Pz[i]);
                truth += Form("  %d (%.0f MeV, cos#theta_{z}=%.2f)", pdg[i], ke, pp > 0 ? Pz[i] / pp : 1.0);
            }
            lines.push_back(truth);
        }
        TString dists = "Closest track endpoint to vertex (3D):";
        for (int trk = 0; trk < nTracks; ++trk) if (trk != primaryTrk) dists += Form("  trk %d %.1f cm%s", trk, vertexDistance(trk), isNear(trk) ? "" : " (out)");
        lines.push_back(dists);
        lines.push_back("The dashed circle is the projected vertex sphere: tracks offset in y look closer than they are.");
        lines.push_back(Form("Tracks with an endpoint within %.0f cm (3D) of the vertex, PID'd by the selection:", VERTEX_RADIUS));
        if (trackLines->GetEntries() == 0) lines.push_back("   none");
        for (int i = 0; i < trackLines->GetEntries(); ++i) lines.push_back("   " + ((TObjString*) trackLines->At(i))->GetString());

        double y = 0.90, step = 0.9 / std::max<size_t>(lines.size(), 8);
        for (const TString& line : lines) {
            TLatex* l = new TLatex(0.01, y, line);
            l->SetNDC(); l->SetTextSize(std::min(0.10, 0.8 * step)); l->SetTextFont(42);
            l->Draw(); keep.push_back(l);
            y -= step;
        }

        TString name = Form("%s%02zu_run%d_sub%d_evt%d.png", dir.Data(), iCand, run, subrun, event);
        c->SaveAs(name);

        log << "=== " << iCand << ": run " << run << " subrun " << subrun << " event " << event
            << " (" << *file << " entry " << entry << ")" << std::endl;
        for (const TString& line : lines) log << line << std::endl;
        log << std::endl;

        delete c;
        delete trackLines;
    }
    log.close();
}
