#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TString.h>
#include <iostream>
#include <set>
#include <tuple>
#include <fstream>
#include <sstream>

// Define type for (run, subrun, event)
typedef std::tuple<int,int,int> RSE;

void skim_by_meta(const char* infile, const char* outfile, const char* meta_list) {
    // ------------------------------------------------------------------
    // 1. Fill allowed (run,subrun,event) set
    //     You should generate a text file from your sqlite db
    //       like: run subrun event (space-separated per line)
    // ------------------------------------------------------------------
    std::set<RSE> allowed;
    {
        std::ifstream fin(meta_list);
        if (!fin.is_open()) {
            std::cerr << "Error: cannot open meta_list\n";
            return;
        }
        std::string line;
        // Skip header
        if (!std::getline(fin, line)) {
            std::cerr << "Error: meta_list is empty or missing header\n";
            return;
        }
        int r, sr, ev;
        while (std::getline(fin, line)) {
            std::istringstream iss(line);
            std::string token;
            // Read run
            if (!std::getline(iss, token, ',')) continue;
            r = std::stoi(token);
            // Read subrun
            if (!std::getline(iss, token, ',')) continue;
            sr = std::stoi(token);
            // Read event
            if (!std::getline(iss, token, ',')) continue;
            ev = std::stoi(token);
            allowed.emplace(r, sr, ev);
        }
        std::cout << "Loaded " << allowed.size() << " RSE entries from CSV list.\n";
    }

    // ------------------------------------------------------------------
    // 2. Open input file + TTree
    // ------------------------------------------------------------------
    // Read a list of input ROOT files from a text file (infile)
    // Each line in infile should be a path to a ROOT file
    std::vector<std::string> input_files;
    {
        std::ifstream fin(infile);
        if (!fin.is_open()) {
            std::cerr << "Error: cannot open input file list " << infile << "\n";
            return;
        }
        std::string line;
        while (std::getline(fin, line)) {
            if (!line.empty()) input_files.push_back(line);
        }
    }
    if (input_files.empty()) {
        std::cerr << "Error: no input files found in " << infile << "\n";
        return;
    }

    // Use TChain to chain all input files
    TChain tin("recTree"); // adjust tree name if needed
    for (const auto& fname : input_files) {
        int nadd = tin.Add(fname.c_str());
        if (nadd == 0) {
            std::cerr << "Warning: could not add file " << fname << " to TChain\n";
        }
    }
    if (tin.GetNtrees() == 0) {
        std::cerr << "Error: no valid trees found in input files.\n";
        return;
    }

    // Branches
    UInt_t run, subrun, event;
    tin.SetBranchAddress("rec.hdr.run",    &run);
    tin.SetBranchAddress("rec.hdr.subrun", &subrun);
    tin.SetBranchAddress("rec.hdr.evt",    &event);

    // ------------------------------------------------------------------
    // 3. Create output file + clone tree structure
    // ------------------------------------------------------------------
    TFile *fout = TFile::Open(outfile, "RECREATE");
    TTree *tout = tin.CloneTree(0, "fast"); // empty clone with same branches
    tout->SetDirectory(fout);

    // ------------------------------------------------------------------
    // 4. Loop over entries and skim
    // ------------------------------------------------------------------
    Long64_t nentries = tin.GetEntries();
    std::cout << "Scanning " << nentries << " entries across " << input_files.size() << " files...\n";

    // Implement a counter for files during skimming
    // We'll count how many unique files contributed at least one event to the output
    std::set<std::string> contributing_files;

    // Add a counter for skimming files to track progress
    std::string lastFileName = "";
    size_t fileCount = 0;
    size_t totalFiles = input_files.size();
    Long64_t lastEntry = -1;

    // For timing
    auto start_time = std::chrono::steady_clock::now();
    auto last_report_time = start_time;

    for (Long64_t i = 0; i < nentries; i++) {
        tin.GetEntry(i);

        // Get the file name for this entry
        TFile* currentFile = tin.GetCurrentFile();
        std::string currentFileName = currentFile ? currentFile->GetName() : "";

        // If we have moved to a new file, increment fileCount and print progress
        if (currentFileName != lastFileName) {
            fileCount++;
            lastFileName = currentFileName;

            // Print progress for file skimming
            auto now = std::chrono::steady_clock::now();
            double elapsed_sec = std::chrono::duration_cast<std::chrono::seconds>(now - start_time).count();
            double avg_per_file = (fileCount > 1) ? elapsed_sec / (fileCount - 1) : 0.0;
            double est_total = avg_per_file * totalFiles;
            double est_remaining = est_total - elapsed_sec;

            if ((i+1) % 500 == 0) {
                std::cout << "[Skim Progress] File " << fileCount << "/" << totalFiles
                        << ": " << currentFileName << std::endl;
                std::cout << "    Elapsed: " << elapsed_sec << "s, "
                          << "Avg/file: " << avg_per_file << "s, "
                          << "Est. remaining: " << (est_remaining > 0 ? est_remaining : 0) << "s"
                          << std::endl;
            }
        }

        RSE key(run, subrun, event);
        if (allowed.count(key)) {
            tout->Fill();
            if (!currentFileName.empty()) {
                contributing_files.insert(currentFileName);
            }
        }

        // Optionally, print entry progress every 100k entries
        // if ((i+1) % 10000 == 0) {{
        //     auto now = std::chrono::steady_clock::now();
        //     double elapsed_sec = std::chrono::duration_cast<std::chrono::seconds>(now - start_time).count();
        //     std::cout << "[Entry Progress] " << (i+1) << "/" << nentries
        //               << " entries (" << 100.0*(i+1)/nentries << "%), "
        //               << "Elapsed: " << elapsed_sec << "s" << std::endl;
        // }
    }

    // ------------------------------------------------------------------
    // 5. Save and close
    // ------------------------------------------------------------------
    fout->cd();
    tout->Write();
    
    // Get entry count before closing
    Long64_t nwritten = tout->GetEntries();
    
    // Close output file first
    fout->Close();
    delete fout;
    
    // Reset TChain to close all input files
    tin.Reset();
    
    std::cout << "Wrote skimmed file: " << outfile 
              << " with " << nwritten << " entries.\n";
    std::cout << "Number of files that contributed at least one event: " 
              << contributing_files.size() << std::endl;
    // if (!contributing_files.empty()) {
    //     std::cout << "Files that contributed:" << std::endl;
    //     for (const auto& fname : contributing_files) {
    //         std::cout << "  " << fname << std::endl;
    //     }
    // }
}
