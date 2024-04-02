//
// Created by MasaHiroSaber on 2023/11/29.
//

#include "CImageProcessingEx.h"
#include "CalFunction.h"
#include "PCAMatrix.h"

using namespace std;
double separability;

CImageProcessingEx::CImageProcessingEx(void)
{}

CImageProcessingEx::~CImageProcessingEx(void)
{}

//直方图均衡化
BOOL CImageProcessingEx::hisTeg(CImageDataset &imgIn, CImageDataset &imgOut)
{
    const int LEVEL = 256;
    int k, row, col;
    double hist[LEVEL], sk[LEVEL];

    if (imgIn.empty())
    {
        return FALSE;
    }

    if (FALSE == imgIn.duplicate(imgOut))
    {
        return FALSE;
    }
    double *data = imgOut.m_data;

    for (k = 0; k < LEVEL; k++)
    {
        hist[k] = 0;
    }
    for (row = 0; row < imgIn.m_ysize; row++)
    {
        for (col = 0; col < imgIn.m_xsize; col++)
        {
            hist[UINT8(data[row * imgIn.m_xsize + col])]++;
        }
    }

    sk[0] = hist[0] / (imgIn.m_ysize * imgIn.m_xsize);
    for (k = 1; k < LEVEL - 1; k++)
    {
        sk[k] = sk[k - 1] + hist[k] / (imgIn.m_ysize * imgIn.m_xsize);
    }
    sk[255] = 1;

    for (row = 0; row < imgOut.m_ysize; row++)
    {
        for (col = 0; col < imgOut.m_xsize; col++)
        {
            for (k = 0; k < LEVEL; k++)
            {
                if (data[row * imgOut.m_xsize + col] == k)
                {
                    data[row * imgOut.m_xsize + col] = int((LEVEL - 1) * sk[k] + 0.5);//灰度变换
                    k = LEVEL;
                }
            }
        }
    }
    return TRUE;
}

//直方图匹配
BOOL CImageProcessingEx::hisMatch(CImageDataset &imgIn, CImageDataset &imgReference, CImageDataset &imgOut)
{
    if (imgIn.empty() || imgReference.empty()) return FALSE;

    const int LEVEL = 256;
    vector<double> imgInHist(LEVEL), imgRefHist(LEVEL);
    vector<double> imgInCumHist(LEVEL), imgRefCumHist(LEVEL);

    if (FALSE == imgOut.create(imgIn.m_xsize, imgIn.m_ysize, imgIn.m_rastercount))
        return FALSE;

    double *imgOutput = imgOut.m_data;
    const double *imgInput = imgIn.m_data;

    int area = imgIn.m_xsize * imgIn.m_ysize;

    for (int band = 0; band < imgIn.m_rastercount; band++)
    {
        CalFunction::calHist(imgIn, band, imgInHist);
        CalFunction::calHist(imgReference, band, imgRefHist);

        CalFunction::calCumHist(imgInHist, imgInCumHist, area, LEVEL);
        CalFunction::calCumHist(imgRefHist, imgRefCumHist, area, LEVEL);

        vector<double> temp(LEVEL);
        vector<double>::iterator minValue;
        int record[LEVEL];

        int i, j;
        for (i = 0; i < LEVEL; i++)
        {
            for (j = 0; j < LEVEL; j++)
            {
                temp[j] = abs(imgInCumHist[i] - imgRefCumHist[j]);
            }
            //#include <algorithm> 来使用 min_element
            minValue = min_element(temp.begin(), temp.end());
            record[i] = distance(temp.begin(), minValue);
        }


        for (int row = 0; row < imgIn.m_ysize; row++)
        {
            for (int col = 0; col < imgIn.m_xsize; col++)
            {
                int index = band * imgIn.m_ysize * imgIn.m_xsize + row * imgIn.m_xsize + col;
                imgOutput[index] = record[int(imgInput[index])];
            }
        }
    }
    return TRUE;
}


BOOL CImageProcessingEx::medFilter(CImageDataset &imgIn, CImageDataset &imgOut, int kerLen)
{
    if (imgIn.empty()) return FALSE;

    if (FALSE == imgOut.create(imgIn.m_xsize, imgIn.m_ysize, imgIn.m_rastercount))
        return FALSE;

    int count = 0;
    int kerArea = kerLen * kerLen;
    int border = kerLen / 2;

    double *imgOutput = imgOut.m_data;
    const double *imgInput = imgIn.m_data;

    CImageDataset Temp;
    CalFunction::zeroFill(imgIn, Temp, kerLen);
//    Temp.create(imgIn.m_xsize + 2 * border, imgIn.m_ysize + 2 * border, imgIn.m_rastercount);
    vector<double> kerValue(kerArea);

    //0填充
    for (int band = 0; band < imgIn.m_rastercount; band++)
    {
        for (int row = border; row < Temp.m_ysize - border; row++)
        {
            for (int col = border; col < Temp.m_xsize - border; col++)
            {
                int indexTemp = band * Temp.m_ysize * Temp.m_xsize + row * Temp.m_xsize + col;
                int indexIn = band * imgIn.m_ysize * imgIn.m_xsize + (row - border) * imgIn.m_xsize + (col - border);

                Temp.m_data[indexTemp] = imgInput[indexIn];
            }
        }
    }

    for (int band = 0; band < Temp.m_rastercount; band++)
    {
        for (int row = border; row < Temp.m_ysize - border; row++)
        {
            for (int col = border; col < Temp.m_xsize - border; col++)
            {
                for (int kerRow = row - border; kerRow <= row + border; kerRow++)
                {
                    for (int kerCol = col - border; kerCol <= col + border; kerCol++)
                    {
                        int index = band * Temp.m_ysize * Temp.m_xsize + kerRow * Temp.m_xsize + kerCol;
                        kerValue[count++] = Temp.m_data[index];
                    }
                }

                int indexOut = band * imgIn.m_ysize * imgIn.m_xsize + (row - border) * imgIn.m_xsize + (col - border);
                sort(kerValue.begin(), kerValue.end());
                imgOutput[indexOut] = kerValue[int(kerArea / 2 + 1)];

                count = 0;
            }
        }
    }
    return TRUE;
}

//双边滤波
BOOL CImageProcessingEx::bilateralFilter(CImageDataset &imgIn, CImageDataset &imgOut, int kerLen, double sigma_S, double sigma_R)
{
    if (imgIn.empty()) return FALSE;

    if (FALSE == imgOut.create(imgIn.m_xsize, imgIn.m_ysize, imgIn.m_rastercount)) return FALSE;

    int count = 0;
    int kerArea = kerLen * kerLen;
    int border = kerLen / 2;

    double Gsigma_S;
    double Gsigma_R;
    double weightP;

    double *imgOutput = imgOut.m_data;
//    const double *imgInput = imgIn.m_data;

    CImageDataset Temp;
    CalFunction::zeroFill(imgIn, Temp, kerLen);

    vector<double> Molecule(kerArea);
    vector<double> Denominator(kerArea);

    for (int band = 0; band < imgIn.m_rastercount; band++)
    {
        for (int row = border; row < Temp.m_ysize - border; row++)
        {
            for (int col = border; col < Temp.m_xsize - border; col++)
            {
                for (int kerRow = row - border; kerRow <= row + border; kerRow++)
                {
                    for (int kerCol = col - border; kerCol <= col + border; kerCol++)
                    {
                        int indexCenter = band * Temp.m_ysize * Temp.m_xsize + row * Temp.m_xsize + col;
                        int indexTemp = band * Temp.m_ysize * Temp.m_xsize + kerRow * Temp.m_xsize + kerCol;

                        Gsigma_S = exp(
                                -((pow(row - kerRow, 2) + pow(col - kerCol, 2)) / (2 * pow(sigma_S, 2))));
                        Gsigma_R = exp(
                                -(pow(Temp.m_data[indexTemp] - Temp.m_data[indexCenter], 2) / (2 * pow(sigma_R, 2))));
                        weightP = Gsigma_S * Gsigma_R;

                        Molecule[count] = weightP * Temp.m_data[indexTemp];
                        Denominator[count] = weightP;
                        count++;
                    }
                }

                int indexOut = band * imgIn.m_ysize * imgIn.m_xsize + (row - border) * imgIn.m_xsize + (col - border);
                double outValue = (accumulate(Molecule.begin(), Molecule.end(), 0.0) /
                                   accumulate(Denominator.begin(), Denominator.end(), 0.0));
                imgOutput[indexOut] = outValue;
                count = 0;
            }
        }
    }
    return TRUE;
}

//拉普拉斯锐化
BOOL CImageProcessingEx::laplaceSharpening(CImageDataset &imgIn, CImageDataset &imgOut)
{
    if (imgIn.empty()) return FALSE;

    if (FALSE == imgOut.create(imgIn.m_xsize, imgIn.m_ysize, imgIn.m_rastercount)) return FALSE;

    int count = 0;
    int border = 1;

    double *imgOutput = imgOut.m_data;
    const double *imgInput = imgIn.m_data;

    vector<double> kerValue(9);
    vector<int> Laplacian = {0, 1, 0,
                             1, -4, 1,
                             0, 1, 0};

    CImageDataset Temp;
    CalFunction::zeroFill(imgIn, Temp, 3);

    for (int band = 0; band < Temp.m_rastercount; band++)
    {
        for (int row = border; row < Temp.m_ysize - border; row++)
        {
            for (int col = border; col < Temp.m_xsize - border; col++)
            {
                for (int kerRow = row - border; kerRow <= row + border; kerRow++)
                {
                    for (int kerCol = col - border; kerCol <= col + border; kerCol++)
                    {
                        int indexTemp = band * Temp.m_ysize * Temp.m_xsize + kerRow * Temp.m_xsize + kerCol;
                        kerValue[count] = Temp.m_data[indexTemp] * Laplacian[count];
                        count++;
                    }
                }
                int indexOut = band * imgIn.m_ysize * imgIn.m_xsize + (row - border) * imgIn.m_xsize + (col - border);
                double outValue = imgInput[indexOut] - accumulate(kerValue.begin(), kerValue.end(), 0.0);
                imgOutput[indexOut] = outValue;

                count = 0;
            }
        }
    }
    return TRUE;
}


//傅里叶变换
BOOL CImageProcessingEx::TwoD_DFT(CImageDataset &imgIn, complex<double> *imgDFTdata)
{
    if (imgIn.empty()) return FALSE;

    const double *imgInput = imgIn.m_data;
    vector<complex<double>> dataDFT(imgIn.m_xsize);

    for (int band = 0; band < imgIn.m_rastercount; band++)
    {
        for (int row = 0; row < imgIn.m_ysize; row++)
        {
            for (int col = 0; col < imgIn.m_xsize; col++)
            {
                int index = band * imgIn.m_ysize * imgIn.m_xsize + row * imgIn.m_xsize + col;
                imgDFTdata[index] = imgInput[index] * pow(-1, (row + col));
            }
        }
    }

    for (int band = 0; band < imgIn.m_rastercount; band++)
    {
        for (int row = 0; row < imgIn.m_ysize; row++)
        {
            memcpy(dataDFT.data(), imgDFTdata + band * imgIn.m_ysize * imgIn.m_xsize + row * imgIn.m_xsize,
                   imgIn.m_xsize * sizeof(complex<double>));

            CalFunction::OneD_DFT(dataDFT, FALSE);

            memcpy(imgDFTdata + band * imgIn.m_ysize * imgIn.m_xsize + row * imgIn.m_xsize, dataDFT.data(),
                   imgIn.m_xsize * sizeof(complex<double>));
        }

        dataDFT.resize(imgIn.m_ysize);

        for (int col = 0; col < imgIn.m_xsize; col++)
        {
            for (int row = 0; row < imgIn.m_ysize; row++)
            {
                dataDFT[row] = imgDFTdata[band * imgIn.m_ysize * imgIn.m_xsize + row * imgIn.m_xsize + col];
            }

            CalFunction::OneD_DFT(dataDFT, FALSE);

            for (int row = 0; row < imgIn.m_ysize; row++)
            {
                imgDFTdata[band * imgIn.m_ysize * imgIn.m_xsize + row * imgIn.m_xsize + col] = dataDFT[row];
            }
        }
    }
    return TRUE;
}

BOOL CImageProcessingEx::TwoD_DFTshow(CImageDataset &imgIn, CImageDataset &imgOut, int doWhat, bool trans_option)
{
    if (imgIn.empty()) return FALSE;

    if (FALSE == imgOut.create(imgIn.m_xsize, imgIn.m_ysize, imgIn.m_rastercount)) return FALSE;

    complex<double> *imgDFTdata = new complex<double>[imgIn.m_rastercount * imgIn.m_xsize * imgIn.m_ysize];

    if (FALSE == TwoD_DFT(imgIn, imgDFTdata)) return FALSE;

    CalFunction::visualDFT(imgDFTdata, imgOut.m_data, doWhat, trans_option, imgIn.m_ysize, imgIn.m_xsize,
                           imgIn.m_rastercount);

    delete[]imgDFTdata;
    imgDFTdata = NULL;

    return TRUE;
}

//傅里叶反变换
BOOL CImageProcessingEx::TwoD_IDFT(complex<double> *imgDFTdata, CImageDataset &imgOut, int rows, int cols, int bands)
{
    if (imgDFTdata == NULL) return FALSE;

    double *imgOutput = imgOut.m_data;

    vector<complex<double>> dataIDFT(cols);

    for (int band = 0; band < bands; band++)
    {
        for (int row = 0; row < rows; row++)
        {
            memcpy(dataIDFT.data(), imgDFTdata + band * rows * cols + row * cols, cols * sizeof(complex<double>));
            CalFunction::OneD_DFT(dataIDFT, TRUE);
            memcpy(imgDFTdata + band * rows * cols + row * cols, dataIDFT.data(), cols * sizeof(complex<double>));
        }

        dataIDFT.resize(rows);

        for (int col = 0; col < cols; col++)
        {
            for (int row = 0; row < rows; row++)
            {
                dataIDFT[row] = imgDFTdata[band * rows * cols + row * cols + col];
            }

            CalFunction::OneD_DFT(dataIDFT, TRUE);

            for (int row = 0; row < rows; row++)
            {
                imgOutput[band * rows * cols + row * cols + col] = dataIDFT[row].real();
            }
        }
    }

    for (int band = 0; band < bands; band++)
    {
        for (int row = 0; row < rows; row++)
        {
            for (int col = 0; col < cols; col++)
            {
                int indexOut = band * rows * cols + row * cols + col;
                imgOutput[indexOut] = imgOutput[indexOut] * pow(-1, ((row + col) % 2));
            }
        }
    }
    return TRUE;
}

BOOL CImageProcessingEx::DFT2IDFT(CImageDataset &imgIn, CImageDataset &imgOut, CImageDataset &imgDif)
{
    if (imgIn.empty()) return FALSE;

    if (FALSE == imgOut.create(imgIn.m_xsize, imgIn.m_ysize, imgIn.m_rastercount)) return FALSE;

    if (FALSE == imgDif.create(imgIn.m_xsize, imgIn.m_ysize, imgIn.m_rastercount)) return FALSE;

    complex<double> *imgDFTdata = new complex<double>[imgIn.m_rastercount * imgIn.m_xsize * imgIn.m_ysize];

    if (FALSE == TwoD_DFT(imgIn, imgDFTdata)) return FALSE;

    if (FALSE == TwoD_IDFT(imgDFTdata, imgOut, imgIn.m_ysize, imgIn.m_xsize, imgIn.m_rastercount)) return FALSE;

    for (int count = 0; count < imgIn.m_xsize * imgIn.m_ysize * imgIn.m_rastercount; count++)
    {
        imgDif.m_data[count] =imgOut.m_data[count] - imgIn.m_data[count];
    }

    return TRUE;
}

BOOL CImageProcessingEx::highpassFilter(CImageDataset &imgIn, CImageDataset &imgOut, int HPFtype, double distance)
{
    if (imgIn.empty()) return FALSE;

    if (FALSE == imgOut.create(imgIn.m_xsize, imgIn.m_ysize, imgIn.m_rastercount)) return FALSE;

    complex<double> *imgHPFdata = new complex<double>[imgIn.m_rastercount * imgIn.m_xsize * imgIn.m_ysize];

    if (FALSE == TwoD_DFT(imgIn, imgHPFdata)) return FALSE;

    enum type
    {
        IDEAL = 0, BUTTERWORTH = 1, GUASSIAN = 2
    };
    double (*filter)(int row, int col, double cut, int rows, int cols);
    switch (HPFtype)
    {
        case IDEAL://理想高通滤波器
            filter = CalFunction::idealHPF;
            break;
        case BUTTERWORTH://butterworth
            filter = CalFunction::butterworthHPF;
            break;
        case GUASSIAN://高斯
            filter = CalFunction::gaussHPF;
            break;
    }

    //做点乘
    for (int band = 0; band < imgIn.m_rastercount; band++)
    {
        for (int row = 0; row < imgIn.m_ysize; row++)
        {
            for (int col = 0; col < imgIn.m_xsize; col++)
            {//遍历全图
                int index = band * imgIn.m_ysize * imgIn.m_xsize + row * imgIn.m_xsize + col;
                imgHPFdata[index] = imgHPFdata[index] * filter(row, col, distance, imgIn.m_ysize, imgIn.m_xsize);
            }
        }
    }

    if (FALSE == TwoD_IDFT(imgHPFdata, imgOut, imgIn.m_ysize, imgIn.m_xsize, imgIn.m_rastercount)) return FALSE;

    delete[]imgHPFdata;
    imgHPFdata = NULL;

    return TRUE;
}

BOOL CImageProcessingEx::stripProcessing(CImageDataset &imgIn, CImageDataset &imgOut, double inclination, double angle)
{
    if (imgIn.empty()) return FALSE;

    if (FALSE == imgOut.create(imgIn.m_xsize, imgIn.m_ysize, imgIn.m_rastercount)) return FALSE;

    complex<double> *stripPrc = new complex<double>[imgIn.m_rastercount * imgIn.m_xsize * imgIn.m_ysize];

    if (FALSE == TwoD_DFT(imgIn, stripPrc)) return FALSE;

    for (int band = 0; band < imgIn.m_rastercount; band++)
    {
        for (int row = 0; row < imgIn.m_ysize; row++)
        {
            for (int col = 0; col < imgIn.m_xsize; col++)
            {
                int index = band * imgIn.m_ysize * imgIn.m_xsize + row * imgIn.m_xsize + col;
                stripPrc[index] = stripPrc[index] *
                                  CalFunction::wedgeFilter(row, col, inclination, angle, imgIn.m_ysize, imgIn.m_xsize);
            }
        }
    }

    if (FALSE == TwoD_IDFT(stripPrc, imgOut, imgIn.m_ysize, imgIn.m_xsize, imgIn.m_rastercount)) return FALSE;

    delete[]stripPrc;
    stripPrc = NULL;

    return TRUE;
}

BOOL CImageProcessingEx::PCA(CImageDataset &imgIn, CImageDataset &imgOut, CImageDataset &imgPCA, CImageDataset &imgDif, int P)
{
    if (imgIn.empty()) return FALSE;

    if (FALSE == imgOut.create(imgIn.m_xsize, imgIn.m_ysize, imgIn.m_rastercount)) return FALSE;

    if (FALSE == imgPCA.create(imgIn.m_xsize, imgIn.m_ysize, imgIn.m_rastercount)) return FALSE;

    if (FALSE == imgDif.create(imgIn.m_xsize, imgIn.m_ysize, imgIn.m_rastercount)) return FALSE;


////////1
    PCAMatrix Matrix_X;
    int pixels = imgIn.m_xsize * imgIn.m_ysize * imgIn.m_rastercount;
    Matrix_X.create(imgIn.m_rastercount, imgIn.m_xsize * imgIn.m_ysize);
    for (int i = 0; i < pixels; i++)
    {
        Matrix_X.m_data[i] = imgIn.m_data[i];
    }
////////2
    PCAMatrix Matrix_S = Matrix_X.covariance();
////////3
    float *eigvalues = new float[Matrix_S.m_rows];
    float *eigvectors = new float[Matrix_S.m_rows * Matrix_S.m_cols];
    Matrix_S.eig(eigvalues, eigvectors);
////////4
    float average;
    float sigma = 0;
    for (int i = 0; i < Matrix_S.m_rows; i++)
    {
        average = 0;
        sigma = 0;

        for (int j = 0; j < Matrix_S.m_cols; j++)
        {
            average += eigvectors[i * Matrix_S.m_cols + j] / Matrix_S.m_cols;
        }

        for (int k = 0; k < Matrix_S.m_cols; k++)
        {
            float temp = average - eigvectors[i * Matrix_S.m_cols + k];
            sigma += temp * temp / Matrix_S.m_cols;
        }

        sigma = sqrt(sigma);

        for (int t = 0; t < Matrix_S.m_cols; t++)
        {
            eigvectors[i * Matrix_S.m_cols + t] = (eigvectors[i * Matrix_S.m_cols + t] - average) / sigma;
        }

    }

    int *index = new int[Matrix_S.m_rows];
    for (int i = 0; i < Matrix_S.m_rows; i++)
    {
        index[i] = -1;
    }

    int head = CalFunction::insertSort(eigvalues, index, Matrix_S.m_rows);

    PCAMatrix Matrix_U(Matrix_S.m_cols, Matrix_S.m_rows);

    int k = 0;
    for (int i = head; i != -1; i = index[i])
    {
        for (int j = 0; j < Matrix_S.m_cols; j++)
        {
            Matrix_U.m_data[j * Matrix_S.m_rows + k] = eigvectors[i * Matrix_S.m_cols + j];
        }
        k++;
    }

    PCAMatrix T_Matrix_U = Matrix_U.T_Matrix();

////////5
    PCAMatrix Matrix_Y = (T_Matrix_U * Matrix_X);

    for (int i = 0; i < pixels; i++)
    {
        imgPCA.m_data[i] = Matrix_Y.m_data[i];
    }
////////6
    double eigSum = 0;
    for (int i = 0; i < Matrix_S.m_cols; i++)
    {
        eigSum += eigvalues[i];
    }

    double pSum = 0;
    int pSum_index = head;
    for (int j = 0; j < P; j++)
    {
        pSum += eigvalues[pSum_index];
        pSum_index = index[pSum_index];
    }

    double varianceCR = pSum / eigSum;
    //afxmessage;

    PCAMatrix Matrix_result(Matrix_X.m_rows, Matrix_Y.m_cols);
    PCAMatrix Matrix_Ui(Matrix_X.m_rows, 1);
    PCAMatrix Matrix_Yi(1, Matrix_Y.m_cols);

    int IPCA_index = head;
    for (int i = 0; i < P; i++)
    {
        for (int j = 0; j < Matrix_X.m_rows; j++)
        {
            Matrix_Ui.m_data[j] = eigvectors[IPCA_index * Matrix_X.m_rows + j];
        }
        IPCA_index = index[IPCA_index];

        for (int j = 0; j < Matrix_Y.m_cols; j++)
        {
            Matrix_Yi.m_data[j] = Matrix_Y.m_data[i * Matrix_Y.m_cols + j];
        }

        PCAMatrix mtemp = Matrix_Ui * Matrix_Yi;
        Matrix_result += mtemp;
    }
////////7
    for (int i = 0; i < pixels; i++)
    {
        imgOut.m_data[i] = Matrix_result.m_data[i];
    }

    for (int count = 0; count < imgIn.m_xsize * imgIn.m_ysize * imgIn.m_rastercount; count++)
    {
        imgDif.m_data[count] = imgOut.m_data[count] - imgIn.m_data[count];
    }

    delete[]index;
    delete[]eigvalues;
    delete[]eigvectors;

    return TRUE;
}

BOOL CImageProcessingEx::RGB2IHS(CImageDataset &imgIn, CImageDataset &imgOut, CImageDataset &imgDif, double I_factor, double H_factor, double S_factor)
{
    if (imgIn.empty()) return FALSE;

    if (FALSE == imgOut.create(imgIn.m_xsize, imgIn.m_ysize, imgIn.m_rastercount)) return FALSE;

    if (FALSE == imgDif.create(imgIn.m_xsize, imgIn.m_ysize, imgIn.m_rastercount)) return FALSE;

    double *inputData = imgIn.m_data;
    double *outputData = imgOut.m_data;
    double *difData = imgDif.m_data;

    double R, G, B;
    double I, H, S;
    double theta;


    for (int row = 0; row < imgIn.m_ysize; row++)
    {
        for (int col = 0; col < imgIn.m_xsize; col++)
        {
            R = inputData[imgIn.m_xsize * row + col];
            G = inputData[imgIn.m_xsize * imgIn.m_ysize + imgIn.m_xsize * row + col];
            B = inputData[2 * imgIn.m_xsize * imgIn.m_ysize + imgIn.m_xsize * row + col];

            theta = acos((0.5 * ((R - G) + (R - B))) / (sqrt((R - G) * (R - G) + (R - B) * (G - B))) + (1 / DBL_MAX)) / PI * 180;
            I = (R + G + B) / 3;
            if (R + G + B == 0) S = 1;
            else S = 1 - (3 / (R + G + B)) * (R > G ? (G > B ? B : G) : (R > B ? B : R));
            H = (B <= G ? theta : 360 - theta);

            I *= I_factor;
            H *= H_factor;
            S *= S_factor;

            if (I > 255) I = 255;
            if (H > 360) H = 360;
            if (S > 1) S = 1;

            if (H >= 0 && H < 120)
            {
                R = I * (1 + S * cos(H / 180 * PI) / cos((60 - H) / 180 * PI));
                B = I * (1 - S);
                G = 3 * I - (R + B);
            } else if (H >= 120 && H < 240)
            {
                R = I * (1 - S);
                G = I * (1 + S * cos((H - 120) / 180 * PI) / cos((180 - H) / 180 * PI));
                B = 3 * I - (R + G);
            } else if (H >= 240 && H <= 360)
            {
                G = I * (1 - S);
                B = I * (1 + S * cos((H - 240) / 180 * PI) / cos((300 - H) / 180 * PI));
                R = 3 * I - (G + B);
            }


            outputData[imgIn.m_xsize * row + col] = R;
            outputData[imgIn.m_xsize * imgIn.m_ysize + imgIn.m_xsize * row + col] = G;
            outputData[2 * imgIn.m_xsize * imgIn.m_ysize + imgIn.m_xsize * row + col] = B;

            int difs = 1;
            difData[imgIn.m_xsize * row + col] = abs(difs * (R - inputData[imgIn.m_xsize * row + col]));
            difData[imgIn.m_xsize * imgIn.m_ysize + imgIn.m_xsize * row + col] = abs(
                    difs * (G - inputData[imgIn.m_xsize * imgIn.m_ysize + imgIn.m_xsize * row + col]));
            difData[2 * imgIn.m_xsize * imgIn.m_ysize + imgIn.m_xsize * row + col] = abs(
                    difs * (B - inputData[2 * imgIn.m_xsize * imgIn.m_ysize + imgIn.m_xsize * row + col]));
        }
    }
    return TRUE;
}


BOOL CImageProcessingEx::otsu(CImageDataset &imgIn, CImageDataset &imgOut, bool otsuSwitch, double thresholdValue)
{
    if (imgIn.empty()) return FALSE;

    if (FALSE == imgOut.create(imgIn.m_xsize, imgIn.m_ysize, 1)) return FALSE;

    const int LEVEL = 256;
    vector<double> imgInHist(LEVEL);

    CalFunction::calHist(imgIn, 0, imgInHist);

    double maxVariance = 0;
    double optimalThreshold = 0;
    double count1;
    double count2;
    double sum1;
    double sum2;

    if (otsuSwitch)
    {
        for (int threshold = 0; threshold < LEVEL; threshold++)
        {
            count1 = 0;
            count2 = 0;
            sum1 = 0;
            sum2 = 0;

            for (int i = 0; i <= threshold; i++)
            {
                count1 += imgInHist[i];
                sum1 += i * imgInHist[i];
            }

            for (int i = threshold + 1; i < LEVEL; i++)
            {
                count2 += imgInHist[i];
                sum2 += i * imgInHist[i];
            }

            double mean1 = (count1 == 0) ? 0.0 : sum1 / count1;
            double mean2 = (count2 == 0) ? 0.0 : sum2 / count2;

            double variance = count1 * count2 * pow(mean1 - mean2, 2);

            if (variance > maxVariance)
            {
                maxVariance = variance;
                optimalThreshold = threshold;
            }
        }
    } else
    {
        optimalThreshold = thresholdValue;
    }

    int pixels = imgIn.m_xsize * imgIn.m_ysize;
    double P_i;
    double totalMean = 0;
    double totalPixels = 0;
    double totalSum = 0;
    for (int i = 0; i < 256; ++i)
    {
        totalPixels += imgInHist[i];
        totalSum += i * imgInHist[i];
    }

    totalMean = totalSum / totalPixels;

    double totalVariance = 0;
    for (int i = 0; i < 256; ++i)
    {
        totalVariance += pow(double(i - totalMean), 2) * imgInHist[i];
    }

    for (int i = 0; i <= optimalThreshold; i++)
    {
        count1 += imgInHist[i];
        sum1 += i * imgInHist[i];
    }

    double partVariance = 0;
    partVariance = pow((totalMean * count1 - sum1), 2) / (count1 * (1 - count1));

    separability = partVariance / totalVariance;

    for (int row = 0; row < imgIn.m_ysize; row++)
    {
        for (int col = 0; col < imgIn.m_xsize; col++)
        {
            if (imgIn.m_data[row * imgIn.m_xsize + col] > optimalThreshold)
                imgOut.m_data[row * imgIn.m_xsize + col] = 255;
            else imgOut.m_data[row * imgIn.m_xsize + col] = 0;
        }
    }
    return TRUE;
}

BOOL CImageProcessingEx::normalThresholdprocessing(CImageDataset &imgIn, CImageDataset &imgOut)
{
    if(imgIn.empty()) return FALSE;

    if()
}








