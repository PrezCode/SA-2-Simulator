//Version 1, 4 November 2025
#include <iostream>
#include <string>
#include <cmath>
using namespace std;

class CoordinateDeterminationSystem{
    public:
    CoordinateDeterminationSystem(double zero){
        rangeInitial = zero; rangeFinal = zero; rangeRate = zero;
        betaInitial = zero; betaFinal = zero; betaRate = zero;
        epsilonInitial = zero; epsilonFinal = zero; epsilonRate = zero;
    }
    //Trackfile calculations
    void initial(double X, double Y, double Z){
        rangeInitial = sqrt(X*X + Y*Y + Z*Z); //Create first slant range reference value
        betaInitial = atan2(Y, X);  //Create first beta reference value
        epsilonInitial = atan2(Z, sqrt(X*X + Y*Y)); //Create first epsilon reference value
        polarInitialX = rangeInitial*cos(betaInitial);
        polarInitialY = rangeInitial*sin(betaInitial);
        heightInitial = rangeInitial*sin(epsilonInitial);
    }
    void final(double X, double Y, double Z){
        rangeFinal = sqrt(X*X + Y*Y + Z*Z); //Create second slant range reference value
        betaFinal = atan2(Y, X);    //Create second beta reference value
        epsilonFinal = atan2(Z, sqrt(X*X + Y*Y));   //Create second epsilon reference value
        polarFinalX = rangeFinal*cos(betaFinal);
        polarFinalY = rangeFinal*sin(betaFinal);
        heightFinal = rangeFinal*sin(epsilonFinal);
    }
    void rates(double dt){
        rangeRate = (rangeFinal - rangeInitial)/dt;
        epsilonRate = (epsilonFinal - epsilonInitial)/dt;
        betaRate = (betaFinal - betaInitial)/dt;
                
        Vxy = abs(sqrt(pow(rangeFinal*betaRate, 2) + pow(rangeFinal*epsilonRate, 2) + pow(rangeRate, 2)));
    }
    void WEZ(){        
        /*if(Vxy <= 640){    //Rmax for Three-Point Method
            if(heightFinal < 1000){Rmax = 15000;}
                else if(heightFinal >= 1000 && heightFinal <= 30000){Rmax = -0.0000070689*pow(heightFinal, 2) + 1.3565*heightFinal + 20824.795;}
                else{Rmax = 0;}
        }else{Rmax = 0;}*/
        if(Vxy <= 640){    //Rmax for Three-Point Method
            if(heightFinal < 1000){Rmax = 15000;}
                else if(heightFinal >= 1000 && heightFinal <= 17500){Rmax = 21473.96*pow(1.0000393, heightFinal);}
                else if(heightFinal > 17500 && heightFinal <= 25000){Rmax = 43000;}
                else{Rmax = 0;}
        }else{Rmax = 0;}
        alpha = atan2(polarFinalY - polarInitialY, polarFinalX - polarInitialX);
        theta = PI/2 - abs(alpha - (betaFinal + PI));
        range_cross = rangeFinal*cos(epsilonFinal)*cos(theta);
        TimeToIntercept = (sqrt(pow(rangeFinal*cos(epsilonFinal), 2) - pow(range_cross, 2)) - sqrt(pow(Rmax, 2) - pow(range_cross, 2)))/Vxy;
    }
    void Fuze(double height){
        if(height >= 5000){fuzeSensitivity = 300;}
            else{fuzeSensitivity = 100;}
    }
    double TTI(){return TimeToIntercept;}
    double getRange(){return rangeFinal;}
    double getGroundRange(){return sqrt(pow(polarFinalX, 2) + pow(polarFinalY, 2));}
    double crossRange(){return range_cross;}
    double maxRange(){return Rmax;}
    double getE(){return epsilonFinal;}
    double getB(){return betaFinal;}
    double getEdot(){return epsilonRate;}
    double getBdot(){return betaRate;}
    double getRangeRate(){return rangeRate;}
    double getFuze(){return fuzeSensitivity;}
    double getHeight(){return heightFinal;}
    private:
    double epsilonInitial, epsilonFinal, epsilonRate, betaInitial, betaFinal, betaRate, rangeInitial, rangeFinal, rangeRate, fuzeSensitivity;
    double polarInitialX, polarFinalX, polarInitialY, polarFinalY, alpha, theta, range_cross, TimeToIntercept, Vxy, Rmax, heightInitial, heightFinal;
};