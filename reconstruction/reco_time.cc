#include <iostream>
#include <deque>
#include <cmath>

#include "TCanvas.h"
#include "TChain.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TVector3.h"

struct Hits
{
    double time;
    TVector3 pos;
    double energy;
};

// Function for reconstructing the Line of Response (LOR)
TVector3 Reconstruction(const Hits& hits1, const Hits& hits2)
{
    if (hits1.energy < 0.4 || hits1.energy > 0.6 ||
        hits2.energy < 0.4 || hits2.energy > 0.6)
    {
        return TVector3(sqrt(-1), sqrt(-1), sqrt(-1));
    }

    double dt = hits1.time - hits2.time;

    if (fabs(dt) > 0.2)
        return TVector3(sqrt(-1), sqrt(-1), sqrt(-1));

    double dist = 0.5 * (299.792458 * dt);

    TVector3 midPoint = 0.5 * (hits2.pos + hits1.pos);

    TVector3 LOR = (hits2.pos - hits1.pos).Unit();

    TVector3 annihilationPoint = midPoint + LOR * dist;

    return annihilationPoint;
}

int main(int argc, char **argv)
{
    gStyle->SetPalette(kRainBow);

    const int nRuns = 70;

    TH2F *hPos =
        new TH2F("hPos", "Reconstructed Position in XY Plane;X [mm]; Y[mm]", 100,
                 -30, 30, 100, -30, 30);

    TH2F *hPosTot =
        new TH2F("hPosTot", "Total Reconstructed Position in XY Plane;X [mm]; Y[mm]", 100,
                 -30, 30, 100, -30, 30);

    TH1F *hDist =
        new TH1F("hDist", "Distance to Center;Distance [mm];Counts", 100, 0, 50);

    hPos->SetMinimum(-1);
    hPosTot->SetMinimum(-1);

    TCanvas *cRun = new TCanvas();

    for (int iRun = 0; iRun < nRuns; iRun++)
    {
        TChain *fHits = new TChain("sorted_tree");

        fHits->AddFile(Form("../digi/sorted_output_run%d.root", iRun));

        Int_t entries = fHits->GetEntries();

        std::cout << "Number of hits: " << entries << std::endl;

        double time;
        double x, y, z;
        double energy;

        fHits->SetBranchAddress("fTime", &time);
        fHits->SetBranchAddress("fEnergy", &energy);
        fHits->SetBranchAddress("fX", &x);
        fHits->SetBranchAddress("fY", &y);
        fHits->SetBranchAddress("fZ", &z);

        std::deque<Hits> hits;

        for (int iEntry = 0; iEntry < entries; iEntry++)
        {
            fHits->GetEntry(iEntry);

            TVector3 pos(x, y, z);

            Hits hit;
            hit.time = time;
            hit.pos = pos;
            hit.energy = energy;

            hits.push_back(hit);

            // Remove old hits from the deque
            while (!hits.empty() && (time - hits.front().time > 1))
            {
                hits.pop_front();
            }

            // Reconstructing LOR for all neighboring hit combinations
            for (size_t iHit = 0; iHit < hits.size(); ++iHit)
            {
                for (size_t jHit = iHit + 1; jHit < hits.size(); ++jHit)
                {
                    TVector3 annihilationPoint = Reconstruction(hits[iHit], hits[jHit]);

                    // Only fill useful hits into histograms
                    if (annihilationPoint.X() == annihilationPoint.X() &&
                        annihilationPoint.Y() == annihilationPoint.Y() &&
                        annihilationPoint.Z() == annihilationPoint.Z())
                    {
                        hPos->Fill(annihilationPoint.X(), annihilationPoint.Y());
                        hPosTot->Fill(annihilationPoint.X(), annihilationPoint.Y());

                        double distToCenter = sqrt(pow(annihilationPoint.X(), 2) +
                                                   pow(annihilationPoint.Y(), 2));
                        hDist->Fill(distToCenter);
                    }
                }
            }
        }

        cRun->cd();
        hPos->Draw("colz");
        cRun->Print(Form("pos_run%d.pdf", iRun));
        cRun->Print(Form("pos_run%d.png", iRun));

        cRun->cd();
        hPosTot->Draw("colz");
        cRun->Print("pos_tot.pdf");
        cRun->Print("pos_tot.png");

        TCanvas *cDist = new TCanvas();
        hDist->Draw();
        cDist->Print(Form("dist_run%d.pdf", iRun));
        cDist->Print(Form("dist_run%d.png", iRun));

        // Reset the histogram for the next run
        hPos->Reset();
        hDist->Reset();
    }
}

