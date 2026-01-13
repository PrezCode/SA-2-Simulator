//Version 1, 8 November 2025
#include <iostream>
#include <string>
#include <cmath>
using namespace std;

class CommandGenerationSystem{
    public:
    CommandGenerationSystem(){}
    void generateSignal(double time, double Em, double Et, double Bm, double Bt, double Etdot, double Btdot){
        if(time <= 14){range_missile = 5000;}
            else if(time > 14 && time < 55){range_missile = 5000 + (730)*(time - 14);}
        //Calculate distance error between missile and target LOS in meters in the Azimuth and Elevation planes, function limites to +/-1200m
        hDeltaE = range_missile*(Em - Et);
            if(hDeltaE > 1200){hDeltaE = 1200;}
                else if(hDeltaE < -1200){hDeltaE = -1200;}
        hDeltaB = range_missile*(Bm - Bt);
            if(hDeltaB > 1200){hDeltaB = 1200;}
                else if(hDeltaB < -1200){hDeltaB = -1200;}
        //Calculate error correction signal to be applied to hDelta based on time
        if(time <= 14){hKD = 3.5;}
            else{hKD = 3.5 - 0.66*(time - 14);}
        //Combine hDelta and hKD to create new distance error signal
        hE = hDeltaE + hKD;
            if(hE > 800){hE = 800;}
                else if(hE < -800){hE = -800;}
        hB = hDeltaB + hKD;
            if(hB > 800){hB = 800;}
                else if(hB < -800){hB = -800;}
        //Determine new value of hE/hB based on smoothing formula
        if(abs(hE) <= 175){hSmoothE = hE;}
            else if(abs(hE) > 175){hSmoothE = 175 + (hE - 175)/6;}
        if(abs(hB) <= 175){hSmoothB = hB;}
            else if(abs(hB) > 175){hSmoothB = 175 + (hB - 175)/6;}
        //Determine dampening value to add to control signal
        hDampE = (2)/(1 + 0.22)*hE;
        hDampB = (2)/(1 + 0.22)*hB;
        //Determine dynamic compensaion error for motor burn time
        if(time > 42){Xt = 2372 + 19.7*time;}
            else{Xt = 1.34*(1770 + 14.7*time);}
        hDynamicE = Xt*Etdot;
        hDynamicB = Xt*Btdot*cos(Et);
        //Generate extra coefficient based on direct distance error.
        Kp = 2.2/4.5;
        //combines all calculations and then multiplies by a conversion factor to turn the value into voltage
        signalEpsilon = k*(hSmoothE + hDampE + hDynamicE + hV)*Kp;
        if(signalEpsilon > 38){signalEpsilon = 38;} 
            else if(signalEpsilon < -38){signalEpsilon = -38;}
        signalBeta = k*(hSmoothB + hDampB + hDynamicB)*Kp;    
        //Generate K1 and K2 channel signals for maneuver controls  
        gamma = (time - 6)*Btdot*sin(Et);
        K1 = (signalBeta*cos(PI/4 - gamma) + signalEpsilon*sin(PI/4 - gamma))/k*0.355;
        K2 = (signalEpsilon*cos(PI/4 - gamma) - signalBeta*sin(PI/4 - gamma))/k*0.355;
    }
    double getK1(){return K1;}
    double getK2(){return K2;}

    private:
    double signalBeta, signalEpsilon, k{0.0855 /*V/m*/}, Kp, Xt, hV{1.3898}, hB, hE, hSmoothB, hSmoothE, hDampE, hDampB, hDeltaE, hDeltaB, hKD, hDynamicE, hDynamicB, range_missile, gamma, K1, K2;
};