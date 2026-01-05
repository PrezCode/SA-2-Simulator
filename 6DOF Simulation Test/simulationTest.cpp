#include <iostream>
#include <fstream>
#include "dynamics.h"
#include <chrono>
using namespace std;

int main(){
    auto start = chrono::system_clock::now();
    int aircraft{1}, missile{2}, stepCounter{0};
    double mass{0.024}, length{0}, radius{0.0254/4}, time{0}, dt{0.001}, simEnd{60}, Cd{0.5}, Aref{M_PI*pow(radius, 2)};
    dynamics testObject(mass, length, radius, aircraft);
    ofstream state ("states.txt");

    testObject.initialConditions(0.001, 0, 0, 0, 0, 0, 0, -M_PI/2, 0, 0, 0, -10000);
    state << time << " ";
    for(int i = 0; i < 12; i++){state << testObject.getInstance(i) << " ";}
    state << endl;
    
    while(stepCounter < simEnd/dt){
        time += dt;
        testObject.iterate(dt, Cd, Aref);
        state << time << " ";
        for(int i = 0; i < 12; i++){state << testObject.getInstance(i) << " ";}
        state << endl;
        if(testObject.getInstance(11) > 0){break;}
        stepCounter++;
    }
    auto end = chrono::system_clock::now();
    chrono::duration<double> runtime = end - start;
    cout << "Runtime: " << runtime.count() << " seconds" << endl;
    cout << "Elapsed Time: " << time << " seconds" << endl;
    cout << "Final Vx: " << testObject.getInstance(0) << endl;
    cout << "Final Vy: " << testObject.getInstance(1) << endl;
    cout << "Final Vz: " << testObject.getInstance(2) << endl;
    system("PAUSE");
}