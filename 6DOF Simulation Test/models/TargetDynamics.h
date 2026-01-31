//Version 1.01, 19 January 2026
#pragma once
#include <iostream>
#include <string>
#include <cmath>
using namespace std;

class Target{
    public:
    Target(double bearing_deg, double range_km, double heading_deg, double speed_mps, double altitude_km){   //Generate initial target information
        Vsum = speed_mps;
        beta = bearing_deg*(M_PI/180.0); //heading to target in rad
        alpha = heading_deg*(M_PI/180.0); //target heading in rad
        epsilon = atan2(altitude_km, range_km); 
        range_ground = range_km*1000.0;
        height = altitude_km*1000.0;         
        X = range_ground*cos(beta);
        Y = range_ground*sin(beta);
        Vx = speed_mps*cos(alpha);
        Vy = speed_mps*sin(alpha);      
    }
    void move(double dt){X += Vx*dt; Y += Vy*dt;}//move target along XY-plane
    double getBeta(){return beta;}
    double getEpsilon(){return epsilon;}
    double getX(){return X;}
    double getY(){return Y;}
    double getZ(){return height;}
    private:
    double range_ground, height, X, Y, Vsum, Vx, Vy, Vz, beta, epsilon, alpha;
};