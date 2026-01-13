#include <string>
#include <iostream>
#include <cmath>
using namespace std;

class B52{
    public:
    B52(){
        cout << "B-52 Model loaded" << endl;
        mass = emptyMass + motorData[2];
        updateInertia();
        initialState();
    }
    void updateInertia(){ //Assume perfect cylinder
        inertiaData[0] = 0.5*mass*pow(radius, 2);
        inertiaData[1] = inertiaData[2] = 0.25*mass*pow(radius, 2) + (mass*pow(length, 2))/12;
    }
    void initialState(){
        modelData[0] = mass;
        for(int i = 0; i < 4; ++i){modelData[i+1] = inertiaData[i];}
        modelData[5] = length;
        modelData[6] = radius;
        for(int i = 0; i < 5; ++i){modelData[i+7] = wingData[i];}
        modelData[10] += wingData[5];

        modelData[12] = tau;
        modelData[13] = thrust;
        modelData[14] = C_DragBase;
        modelData[15] = C_LiftBase;
        modelData[16] = 0;
    }
    void newState(){
        modelData[0] = mass;
        for(int i = 0; i < 4; ++i){modelData[i+1] = inertiaData[i];}
        modelData[5] = length;
        modelData[6] = radius;
        modelData[13] = thrust;
        modelData[15] = C_LiftBase;
        modelData[16] = CG_CP_avg;
    }
    double transferData(int i){
        return modelData[i];
    }
    const double deg_rad{PI/180}, slug_kg{14.5939}, inch_m{0.0254}, sqft_sqm{0.092903}, ft_m{0.304878}, slugsqft_kgsqm{1.3558179619},
    wingData[6] = {185*ft_m/*wing span*/, 21.62*ft_m /*wing chord*/, 200*sqft_sqm /*Aileron Wing Area*/, 4000*sqft_sqm /*Main Wing Area*/, 200*sqft_sqm /*Rudder Wing Area*/, 1004*sqft_sqm/*Booster Wing Area*/},
    emptyMass = 83250;
    double mass{0}, length{160.9*ft_m}, radius{10*ft_m}, tau{0.5}, thrust{400000}, C_LiftBase{0.2}, C_DragBase{0.0119}, CG_CP_avg{31.7754},
        modelData[17],
        motorData[3] = {0/*Booster Fuel*/, 0/*Booster empty mass*/, 70000/*Sustainer Fuel*/},
        inertiaData[4] = {0 /*Ixx*/, 0 /*Iyy*/, 0 /*Izz*/, 0 /*Ixz*/},
        actuatorMaxAngles[3] = {12*deg_rad /*aileron*/, 28*deg_rad /*elevator*/, 28*deg_rad /*rudder*/};
};