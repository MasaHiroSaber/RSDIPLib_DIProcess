//
// Created by MasaHiroSaber on 2023/12/5.
//

#include "CalFunction.h"

void CalFunction::calHist(CImageDataset &imgIn, int band, vector<double> &hist)
{
    for (int row = 0; row < imgIn.m_ysize; row++)
    {
        for (int col = 0; col < imgIn.m_xsize; col++)
        {
            int index = band * imgIn.m_xsize * imgIn.m_ysize + row * imgIn.m_xsize + col;
            int hisValue = UINT8(imgIn.m_data[index]);
            hist[hisValue]++;
        }
    }
}

void CalFunction::calCumHist(vector<double> &hist, vector<double> &cumHist, int area, int LEVEL)
{
    cumHist[0] = hist[0] / area;
    for (int i = 1; i < LEVEL - 1; i++)
    {
        cumHist[i] = cumHist[i - 1] + hist[i] / area;
    }
    cumHist[255] = 1;
}

void CalFunction::zeroFill(CImageDataset &imgIn, CImageDataset &Temp, int kerLen)
{
    int border = kerLen / 2;

    const double *imgInput = imgIn.m_data;

    CImageDataset Fill;
    Fill.create(imgIn.m_xsize + 2 * border, imgIn.m_ysize + 2 * border, imgIn.m_rastercount);

    //0填充
    for (int band = 0; band < imgIn.m_rastercount; band++)
    {
        for (int row = border; row < Fill.m_ysize - border; row++)
        {
            for (int col = border; col < Fill.m_xsize - border; col++)
            {
                int indexTemp = band * Fill.m_ysize * Fill.m_xsize + row * Fill.m_xsize + col;
                int indexIn = band * imgIn.m_ysize * imgIn.m_xsize + (row - border) * imgIn.m_xsize + (col - border);
                Fill.m_data[indexTemp] = imgInput[indexIn];
            }
        }
    }
    Fill.duplicate(Temp);
}


void CalFunction::OneD_DFT(vector<complex<double>> &complexDFT, bool inverse)
{
    complex<double> Temp;
    vector<complex<double>> complex_copy;
    //vector<complex<double>> complex_copy;
    size_t complexSize = complexDFT.size();
    complex_copy.assign(complexDFT.begin(), complexDFT.end());

    for (size_t k = 0; k < complexSize; k++)
    {
        Temp = 0;
        for (size_t n = 0; n < complexSize; n++)
        {
            complex<double> calTemp(cos((2 * PI * k * n) / complexSize),
                                    ((inverse ? 1 : -1) * sin((2 * PI * k * n) / complexSize)));
            Temp += complex_copy[n] * calTemp;
        }
        complexDFT[k] = (inverse ? Temp / (double) complexSize : Temp);
    }
}


void CalFunction::visualDFT(complex<double> *imgDFTdata, double *imgOutput, int doWhat, bool trans_option, int rows, int cols, int bands)
{

    double (*trans_fun)(const complex<double> &x);//定义函数指针
    enum doWhat
    {
        FREQUENCY = 0, ANGLE = 1, ENERGY = 2
    };
    switch (doWhat)
    {
        case FREQUENCY://显示频谱
        default:
            trans_fun = CalFunction::amplitude;
            break;
        case ANGLE://显示相角
            trans_fun = CalFunction::angle;
            break;
        case ENERGY://显示能量
            trans_fun = CalFunction::energy;
            break;
    }

    //取模的同时找最大模和最小模
    for (int band = 0; band < bands; band++)
    {
        double max = 0;
        double min = DBL_MAX;
        for (int i = 0; i < rows * cols; i++)
        {
            int index = band * rows * cols + i;
            imgOutput[index] = trans_fun(imgDFTdata[index]);//转换为幅度，相角，或者功率

            //记录当前最大值和最小值
            if (imgOutput[index] > max)
                max = imgOutput[index];
            if (imgOutput[index] < min)
                min = imgOutput[index];
        }

        if (trans_option)
        {
            //对数变换
            double c = 255 / log(max + 1);
            for (int i = 0; i < rows * cols; i++)
            {
                imgOutput[band * rows * cols + i] =
                        c * log(imgOutput[band * rows * cols + i] + 1);
            }
        } else
        {
            //线性拉伸
            double extr_diff = max - min;//极差
            for (int i = 0; i < rows * cols; i++)
            {
                imgOutput[band * rows * cols + i] =
                        (imgOutput[band * rows * cols + i] - min) / (extr_diff * 255);
            }
        }
    }
}

double CalFunction::amplitude(const complex<double> &x)
{
    return sqrt(pow(x.imag(), 2) + pow(x.real(), 2));
}

double CalFunction::angle(const complex<double> &x)
{
    return atan(x.imag() / x.real());
}

double CalFunction::energy(const complex<double> &x)
{
    return x.imag() * x.imag() + x.real() * x.real();
}

//void CalFunction::typeHPF(complex<double> *imgHPFdata, int rows, int cols, int bands, int type, double distance)
//{
//    double (*filter)(int row, int col, double cut, int rows, int cols);
//    enum filter_type
//    {
//        IDEAL = 0, BUTTERWORTH = 1, GUASS = 2
//    };
//    switch (type)
//    {
//        case IDEAL://理想高通滤波器
//            filter = idealHPF;
//            break;
//        case BUTTERWORTH://butterworth
//            filter = butterworthHPF;
//            break;
//        case GUASS://高斯
//            filter = gaussHPF;
//            break;
//    }
//
//    //做点乘
//    for (int band = 0; band < bands; band++)
//    {
//        for (int row = 0; row < rows; ++row)
//        {
//            for (int col = 0; col < cols; ++col)
//            {//遍历全图
//                int index = band * rows * cols + row * cols + col;
//                imgHPFdata[index] = imgHPFdata[index] * filter(row, col, distance, rows, cols);
//            }
//        }
//    }
//}

double CalFunction::idealHPF(int row, int col, double distance, int rows, int cols)
{
    return sqrt(double(pow(row - rows / 2, 2)) + double(pow(col - (cols / 2), 2))) > distance ? 0 : 1;
}

double CalFunction::butterworthHPF(int row, int col, double distance, int rows, int cols)
{
    double pointDist = sqrt(double(pow(row - rows / 2, 2)) + double(pow(col - (cols / 2), 2)));
    return 1 / (1 + pow(distance / pointDist, 4));
}

double CalFunction::gaussHPF(int row, int col, double distance, int rows, int cols)
{
    double pointDist = sqrt(double(pow(row - rows / 2, 2)) + double(pow(col - (cols / 2), 2)));
    return 1 - exp(-pow(pointDist, 2) / (2 * pow(distance, 2)));
}

double CalFunction::wedgeFilter(int row, int col, double inclination, double angle, int rows, int cols)
{
    double col_mid = double(cols) / 2;
    double row_mid = double(rows) / 2;
    if (col == col_mid || col == 0) return 1;
    double ang[9];
    ang[0] = row / col;                    //左上角
    ang[1] = row / (col - col_mid);        //上中点
    ang[2] = row / (col - 2 * col_mid);    //右上角
    ang[3] = (row - row_mid) / col;        //左边中点
    ang[4] = (row - row_mid) / (col - col_mid);        //正中心
    ang[5] = (row - row_mid) / (col - 2 * col_mid);    //右边中点
    ang[6] = (row - 2 * row_mid) / col;                //左下角
    ang[7] = (row - 2 * row_mid) / (col - col_mid);        //下边中点
    ang[8] = (row - 2 * row_mid) / (col - 2 * col_mid);    //右下角

    for (int i = 0; i < 9; i++)
    {
        if (ang[i] < tan((inclination + 90 + angle / 2) / 180 * PI) &&
            (ang[i] > tan((inclination + 90 - angle / 2) / 180 * PI)))
            return 0;
    }
    return 1;
}

int CalFunction::insertSort(float *arr, int *link, int length)
{
    int head = 0;
    for (int i = 1; i < length; i++) //i在目标数组中顺序遍历
    {
        if (arr[i] > arr[head]) //如果当前位置的值小于head位置的值
        {
            link[i] = head;//将head连到当前位置的后面
            head = i;//将当前位置作为新的head
        } else //如果当前位置的值大于head位置的值
        {
            int j = head;//从head开始找
            while ((link[j] != -1) && (arr[link[j]] > arr[i]))
                j = link[j];//从head顺着link往后找，直到找到未连接的或者找到比当前值小的最大值为止

            if (link[j] == -1)//如果找到了未连接的
                link[j] = i;

            else //如果找到的是比当前值小的最大值
            {//则将当前值插到这个值后面

                link[i] = link[j];//相当于 curr->next = p->next;
                link[j] = i;      //相当于    p->next = curr;
            }
        }
    }
    return head;
}









