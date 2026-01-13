#include <iostream>
#include <cmath>
#include <string>


class B52Autopilot{
    public:
    enum TrimOption{Moment_Equilibrium, Steady_Glide};
    B52Autopilot(){
        VmagTarget = 200;
        heightTarget = 10000*ft_m;
        flightPathTarget = 0;
        rollTarget = 0;
        pitchTarget = 0;
        yawTarget = 0;
        pTarget = 0;
        qTarget = 0;
        rTarget = 0;
        aileronLimit = 45;
        elevatorLimit = 45;
        rudderLimit = 45;
        
    }
    void autopilotGenerate(double* aeroData, double dt, TrimOption option){
        rollRate = aeroData[0];
        pitchRate = aeroData[1];
        yawRate = aeroData[2];
        rollAngle = wrapPi(aeroData[3]);
        pitchAngle = wrapPi(aeroData[4]);
        yawAngle = wrapPi(aeroData[5]);
        delA = aeroData[6];
        delE = aeroData[7];
        delR = aeroData[8];
        AoA = aeroData[9];
        sideSlip = aeroData[10];
        Vmag = aeroData[11];
        height = aeroData[12];
        switch(option){
            case Moment_Equilibrium: J = pow(rollRate, 2) + pow(pitchRate, 2) + pow(yawRate, 2);
            case Steady_Glide: J = 0;
        }
    }
    double autopilotCommands(int i){
        return commands[i];
    }
    double wrapPi(double a){
        while (a >  PI) a -= 2.0*PI;
        while (a < -PI) a += 2.0*PI;
        return a;
    }
    private:
    double commands[3] = {0,0,0}, Pq, deg_rad{PI/180}, rad_deg{1/deg_rad}, AoA, ft_m{0.304878}, J,
    rollAngle, pitchAngle, yawAngle, 
    rollRate, pitchRate, yawRate, 
    rollTarget, pitchTarget, yawTarget, 
    pTarget, qTarget, rTarget,
    delA, delE, delR,
    aileronLimit, elevatorLimit, rudderLimit,
    VmagTarget, heightTarget, flightPathTarget,
    sideSlip, Vmag, height;
};