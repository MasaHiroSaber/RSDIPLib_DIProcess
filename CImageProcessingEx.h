//
// Created by MasaHiroSaber on 2023/11/29.
//

#ifndef DIP_CIMAGEPROCESSINGEX_H
#define DIP_CIMAGEPROCESSINGEX_H

#include <minwindef.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include <complex>
#include "RSDIPLib.h"


using namespace std;

class CImageProcessingEx :
        public CImageProcessing
{
public:
    CImageProcessingEx(void);

    ~CImageProcessingEx(void);

    static BOOL hisTeg(CImageDataset &imgIn, CImageDataset &imgOut);

    static BOOL hisMatch(CImageDataset &imgIn, CImageDataset &imgReference, CImageDataset &imgOut);

    static BOOL medFilter(CImageDataset &imgIn, CImageDataset &imgOut, int kerLen);

    static BOOL bilateralFilter(CImageDataset &imgIn, CImageDataset &imgOut, int kerLen, double sigma_S, double sigma_R);

    static BOOL laplaceSharpening(CImageDataset &imgIn, CImageDataset &imgOut);

    static BOOL TwoD_DFT(CImageDataset &imgIn, complex<double> *imgDFTdata);

    static BOOL TwoD_DFTshow(CImageDataset &imgIn, CImageDataset &imgOut, int doWhat, bool trans_option);

    static BOOL TwoD_IDFTshow(CImageDataset &imgIn, CImageDataset &imgOut);

    static BOOL TwoD_IDFT(complex<double> *imgDFTdata, CImageDataset &imgOut, int rows, int cols, int bands);

    static BOOL DFT2IDFT(CImageDataset &imgIn, CImageDataset &imgOut, CImageDataset &imgDif);

    static BOOL highpassFilter(CImageDataset &imgIn, CImageDataset &imgOut, int HPFtype, double distance);

    static BOOL stripProcessing(CImageDataset &imgIn, CImageDataset &imgOut, double inclination, double angle);

    static BOOL PCA(CImageDataset &imgIn, CImageDataset &imgOut, CImageDataset &imgPCA, CImageDataset &imgDif, int P);

    static BOOL RGB2IHS(CImageDataset &imgIn, CImageDataset &imgOut, CImageDataset &imgDif,double I_factor, double H_factor, double S_factor);

    static BOOL otsu(CImageDataset &imgIn, CImageDataset &imgOut, bool otsuSwitch, double thresholdValue);

    static BOOL normalThresholdprocessing(CImageDataset & imgIn, CImageDataset &imgOut);
};


#endif //DIP_CIMAGEPROCESSINGEX_H
