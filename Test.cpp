#include "iostream"
#include "cmath"
#define PI 3.1415926535
using namespace std;
int main()
{
    double R,G,B,I,H,S;
    R = 255;
    G = 0;
    B = 1;

    double theta = acos((0.5 * ((R - G) + (R - B))) / (sqrt((R - G) * (R - G) + (R - B) * (G - B)) + (1 / DBL_MAX))) /  PI * 180;
    cout << "cos" << cos(180 / 180 * PI) << endl;
    cout << "theta:" << theta << endl;
    I = (R + G + B) / 3;
    if (R + G + B == 0) S = 1;
    else S = 1 - (3 / (R + G + B)) * (R > G ? (G > B ? B : G) : (R > B ? B : R));

    H = (B <= G ? theta : 360 - theta);


    I *= 1;
    H *= 0.5;
    S *= 1;

    if (I > 255) I = 255;
    if (H > 360) H = 360;
    if (S > 1) S = 1;

    cout <<"I "<<I<<endl;
    cout <<"H "<<H<<endl;
    cout <<"S "<<S<<endl;

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

    cout <<"R "<<R<<endl;
    cout <<"G "<<G<<endl;
    cout <<"B "<<B<<endl;
}