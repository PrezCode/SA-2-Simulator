#include <string>
#include <iostream>
#include <cmath>
using namespace std;

class V755DSU{
    public:
    V755DSU(){
        cout << "V-755DSU Model Loaded" << endl;
        updateMass(0);
        updateInertia();
        initialState();
    }
    void updateMissile(double TOF, double dt){
        motor(TOF, dt);
        updateMass(TOF);
        if(TOF == 3.1){updateDimensions();}
        updateInertia();
        newState();
    }
    void motor(double TOF, double dt){  //Thrust level then fuel mass change
        if(TOF <= 3.0){thrust = 460000.0 - 130000.0*(TOF/3); motorData[0] -= 174*dt;}  //Boost main phase
            else if(TOF > 3.0 && TOF <= 3.5){thrust = 1084821.429*TOF*TOF + -7710625*TOF + 13700714.29; motorData[0] -= 174*dt;} //Boost transition
            else if(TOF > 3.5 && TOF <= 4.0){thrust = 35000.0*(TOF-3.5); motorData[2] -= 723.4/42*dt;}    //Sustainer throttle up
            else if(TOF > 4.0 && TOF <= 45.0){thrust = 35000.0; motorData[2] -= 723.4/42*dt;}    //Sustainer Max
            else if(TOF > 45.0 && TOF <= 45.5){thrust = -17500.0*(TOF-44.0) + 35000.0; motorData[2] -= 723.4/42*dt;}   //Sustainer transition
            else{thrust = 0.0;}   //Glide Phase 
    }
    void updateDimensions(){
        length = 8.215; 
        radius = 0.25;
        modelData[10] = wingData[3];
    }
    
    void updateMass(double TOF){
        if(TOF < 3.1){mass = emptyMass + motorData[1];}
            else{mass = emptyMass;}        
        mass += motorData[0] + motorData[2];
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
        modelData[16] = CG_CP_avg;
    }
    void newState(){
        modelData[0] = mass;
        for(int i = 0; i < 4; ++i){modelData[i+1] = inertiaData[i];}
        modelData[5] = length;
        modelData[6] = radius;
        modelData[13] = thrust;
        modelData[14] = C_DragBase;
        modelData[15] = C_LiftBase;
        modelData[16] = CG_CP_avg;
    }
    double transferData(int i){
        return modelData[i];
    }
    const double deg_rad{PI/180.0}, slug_kg{14.5939}, inch_m{0.0254}, sqft_sqm{0.092903}, ft_m{0.304878}, slugsqft_kgsqm{1.3558179619},
    wingData[6] = {1.691/*wing span*/, 1.0685 /*wing chord*/, 0.0442848 /*Aileron Wing Area*/, 1.2554875 /*Main Wing Area*/, 0.161876 /*Rudder Wing Area*/, 2.027676/*Booster Wing Area*/},
    emptyMass = 667.0;
    double mass{0}, length{10.778}, radius{0.327}, tau{0.1}, thrust{0}, C_LiftBase{0.01}, C_DragBase{0.27}, CG_CP_avg{0.8065},
        modelData[17],
        motorData[3] = {609.0/*Booster Fuel*/, 398.5/*Booster empty mass*/, 723.4/*Sustainer Fuel*/},
        inertiaData[4] = {0 /*Ixx*/, 0 /*Iyy*/, 0 /*Izz*/, 0 /*Ixz*/};
};