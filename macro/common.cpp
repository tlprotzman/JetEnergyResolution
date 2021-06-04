#include <TROOT.h>
#include <TMath.h>

#include <string>
#include <list>
#include <fstream>

// Jet Regions
const int CENTRAL = 0;
const int FORWARD = 1;
const int BACKWARD = 2;

const uint8_t NUM_REGIONS = 3;

// TODO Doing this with iterators could clean things up

class jetData {
    public:
        bool loaded = false;
        std::string descriptiveName;
        float minEta;
        float maxEta;
        std::list<std::string> files;
        int color;
        int secondaryColor;
        int marker;
        float markerSize;
};

// Translate file list into list of file paths
int readFileList(std::string fileList, std::list<std::string> &list) {
    int numFiles = 0;
    std::ifstream files(fileList);
    std::string filePath;
    
    while (std::getline(files, filePath)) {
        list.push_back(filePath);
        numFiles++;
    }
    files.close();
    return numFiles;
}

// pos = [truthEta, truthPhi, recoEta, recoPhi]
// returns R2 = dEta * dEta + dPhi * dPhi
// Wraps phi 
float calculateDistance(float *pos) {
    for (uint8_t i = 0; i < 4; i++) {
        if (std::isnan(pos[i])) {
            return 9999;
        }
    }
    float dEta, dPhi;
    dEta = pos[0] - pos[2];
    dPhi = pos[1] - pos[3];
    if (dPhi > TMath::Pi()) {
        dPhi -= TMath::TwoPi();
        pos[3] += TMath::TwoPi();
    }
    if (dPhi < -1 * TMath::Pi()) {
        dPhi += TMath::TwoPi();
        pos[3] -= TMath::TwoPi();
    }
    return dEta *dEta + dPhi * dPhi;
}