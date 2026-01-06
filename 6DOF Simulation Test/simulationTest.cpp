#include <iostream>
#include <cmath>
#include <fstream>
#include "dynamics.h"
#include <chrono>
using namespace std;

int main(){
    auto start = chrono::system_clock::now();
    int stepCounter{0};
    double  time{0}, dt{0.001}, simEnd{60}, deg_rad{M_PI/180}, slug_kg{14.5939}, inch_m{0.0254}, sqft_sqm{0.092903}, ft_m{0.304878};
    double mass{0.155404754*slug_kg}, 
            length{8*inch_m}, width{4*inch_m}, height{2.25*inch_m}, radius{6*inch_m}, Aref{0.1963495*sqft_sqm},
            wingSpan{0.33333*ft_m}, wingChord{0.66667*ft_m},
            coefficients[6] = {0 /*Cd*/, -1 /*Clp*/, 0 /*Clr*/, -1 /*Cmq*/, 0 /*Cnp*/, -1 /*Cnr*/},
            u{1/pow(10,8)}, v{0}, w{0},
            p{10*deg_rad}, q{20*deg_rad}, r{30*deg_rad},
            phi{0}, theta{0}, psi{0},
            x{0}, y{0}, z{-30000*ft_m};
    Dynamics testObject(mass, length, width, height, radius, Dynamics::Shape::RectPrism);
    ofstream state ("states.txt");

    testObject.initialConditions(u, v, w, p, q, r, phi, theta, psi, x, y, z);

    state << time << " ";
    for(int i = 0; i < 12; i++){state << testObject.getInstance(i) << " ";}
    for(int i = 0; i < 4; i++){state << testObject.getAirData(i) << " ";}
    state << endl;
    
    while(stepCounter < simEnd/dt){
        time += dt;
        testObject.iterate(dt, Aref, wingSpan, wingChord, coefficients);

        state << time << " ";
        for(int i = 0; i < 12; i++){state << testObject.getInstance(i) << " ";}
        for(int i = 0; i < 4; i++){state << testObject.getAirData(i) << " ";}
        state << endl;
        if(testObject.getInstance(11) > 0 && stepCounter > 0){break;}
        stepCounter++;
    }
    auto end = chrono::system_clock::now();
    chrono::duration<double> runtime = end - start;
    cout << "Real Runtime: " << runtime.count() << " seconds" << endl;
    cout << "Sim Runtime: " << time << " seconds" << endl;
    return 0;
    //system("PAUSE");
}