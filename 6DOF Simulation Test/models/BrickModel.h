#include <string>
#include <iostream>
#include <cmath>
using namespace std;

class BrickModel{
    public:
    enum DataSet{Wings, Inertia, Measurements, Coefficients};
    BrickModel(){
        cout << "Standard US Brick" << endl;
    }
    double getData(DataSet category, int i){
        switch(category){
            case Wings: return wingData[i]; break;
            case Inertia: return inertiaData[i]; break;
            case Measurements: return objectData[i]; break;
            case Coefficients: return coefficients[i]; break;
            default: return 0;
        }
    }
    private:
    const double deg_rad{PI/180}, slug_kg{14.5939}, inch_m{0.0254}, sqft_sqm{0.092903}, ft_m{0.304878}, slugsqft_kgsqm{1.3558179619},
    wingData[4] = {0.33333*ft_m /*wing span*/, 0.66667*ft_m /*wing chord*/, 0.22222*sqft_sqm /*Aref*/, 0*sqft_sqm /*Wing Area*/},
    inertiaData[4] = {0.001894220*slugsqft_kgsqm /*Ixx*/, 0.006211019*slugsqft_kgsqm /*Iyy*/, 0.007194665*slugsqft_kgsqm /*Izz*/, 0*slugsqft_kgsqm /*Ixz*/},
    objectData[5] = {0.155404754*slug_kg /*mass*/, 8*inch_m /*length*/, 4*inch_m /*width*/, 2.25*inch_m /*height*/, 0.5 /*radius*/},
    coefficients[6] = {0.01 /*Cd*/, -1 /*Clp*/, 0 /*Clr*/, -1 /*Cmq*/, 0 /*Cnp*/, -1 /*Cnr*/};  
};