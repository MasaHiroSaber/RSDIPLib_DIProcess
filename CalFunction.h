//
// Created by MasaHiroSaber on 2023/12/5.
//

#ifndef DIPROCESS_CALFUNCTION_H
#define DIPROCESS_CALFUNCTION_H
#define PI 3.1415926535
#include <minwindef.h>


#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include <complex>
#include "RSDIPLib.h"

using namespace std;

class CalFunction
{
public:
    static void calHist(CImageDataset &imgIn, int band, vector<double> &hist);

    static void calCumHist(vector<double> &hist, vector<double> &cumHist, int area, int LEVEL);

    static void zeroFill(CImageDataset &imgIn, CImageDataset &Temp, int kerLen);

    static void OneD_DFT(vector<complex<double>> &complexDFT, bool inverse);

    static void
    visualDFT(complex<double> *imgDFTdata, double *imgOutput, int doWhat, bool trans_option, int rows, int cols,
              int bands);
    static void typeHPF(complex<double> *imgHPFdata, int rows, int cols, int bands, int type, double distance);
    static inline double amplitude(const complex<double> &x);//幅度
    static inline double angle(const complex<double> &x);//相角
    static inline double energy(const complex<double> &x);//功率
    static inline double idealHPF(int row, int col, double distance, int rows, int cols);
    static inline double butterworthHPF(int row, int col, double distance, int rows, int cols);
    static inline double gaussHPF(int row, int col, double distance, int rows, int cols);
    static inline double wedgeFilter(int row, int col, double inclination, double angle, int rows, int cols);
    static int insertSort(float *eigvalues, int *index, int length);
};



#endif //DIPROCESS_CALFUNCTION_H
