//Version 1, 20 January 2016
#pragma once
#include <iostream>
#include <string>
#include <cmath>
using namespace std;

class CommandGenerationSystem{
    public:
    CommandGenerationSystem(){}
    void generateSignal(double time, double Em, double Et, double Bm, double Bt, double Etdot, double Btdot){
        if(time <= 14.0){range_missile = 5000.0;}
            else if(time > 14.0 && time < 55.0){range_missile = 5000.0 + (775.0)*(time - 14.0);}
        //Calculate distance error between missile and target LOS in meters in the Azimuth and Elevation planes, function limites to +/-1200m
        hDeltaE = clamp(range_missile*(Em - Et), -1200.0, 1200.0);
        hDeltaB = clamp(range_missile*(Bm - Bt), -1200.0, 1200.0);
        //Calculate error correction signal to be applied to hDelta based on time
        if(time <= 14.0){hKD = 3.5;}
            else{hKD = 3.5 - 0.66*(time - 14.0);}
        //Combine hDelta and hKD to create new distance error signal
        hE = clamp(hDeltaE + hKD, -800.0, 800.0);
        hB = clamp(hDeltaB, -800.0, 800.0);
        //Determine new value of hE/hB based on smoothing formula
        if(abs(hE) <= 175.0){hSmoothE = hE;}
            else if(abs(hE) > 175.0){hSmoothE = 175.0 + (hE - 175.0)/6.0;}
        if(abs(hB) <= 175.0){hSmoothB = hB;}
            else if(abs(hB) > 175.0){hSmoothB = 175.0 + (hB - 175.0)/6.0;}
        //Determine dampening value to add to control signal
        hDampE = 2.0/(1 + 0.22)*hE;
        hDampB = 2.0/(1 + 0.22)*hB;
        //Determine dynamic compensaion error for motor burn time
        if(time <= 42){Xt = 1.34*(1770 + 14.7*time);}
            else{Xt = 3199.116;}
        hDynamicE = Xt*Etdot;
        hDynamicB = Xt*Btdot*cos(Et);
        //Generate extra coefficient based on direct distance error.
        Kp = 2.2/4.5;
        //combines all calculations and then multiplies by a conversion factor to turn the value into voltage
        lambdaEpsilon = clamp(k*(hSmoothE + hDampE + hDynamicE)*Kp + hV, -38.0, 38.0);
        lambdaBeta = clamp(k*(hSmoothB + hDampB + hDynamicB)*Kp, -170.0, 170.0);    
        //Generate K1 and K2 channel signals for maneuver controls  
        if(time >= 6.0){gamma = Btdot*sin(Et)*time - Btdot*sin(Et)*6.0;}
        clamp(gamma, -30.0*d2r, 30.0*d2r);
        K1 = lambdaEpsilon;
        K2 = lambdaBeta;

        //K1 = (lambdaBeta*cos(M_PI/4.0) + lambdaEpsilon*sin(M_PI/4.0));
        //K2 = (lambdaEpsilon*cos(M_PI/4.0) - lambdaBeta*sin(M_PI/4.0));
    }
    double rangeAmplifier(){return 3.6*(range_missile/1000.0);}

    double getK1(){return K1;}
    double getK2(){return K2;}
    double getRangeAmp(){return rangeAmplifier();}
    private:
    double d2r{M_PI/180.0}, lambdaBeta, lambdaEpsilon, k{0.0855 /*V/m*/}, Kp, Xt, hV{1.3898}, hB, hE, hSmoothB, hSmoothE, hDampE, hDampB, hDeltaE, hDeltaB, hKD, hDynamicE, hDynamicB, range_missile, gamma, K1{0.0}, K2{0.0};
    static double clamp(double x, double low, double hi) {return (x < low) ? low : (x > hi) ? hi : x;}
};