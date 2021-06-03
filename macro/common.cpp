#include <TROOT.h>
#include <TMath.h>

#include <string>
#include <list>
#include <fstream>

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
    }
    if (dPhi < -1 * TMath::Pi()) {
        dPhi += TMath::TwoPi();
    }
    return dEta *dEta + dPhi * dPhi;
}