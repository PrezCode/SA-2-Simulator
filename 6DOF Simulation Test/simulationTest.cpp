#include <iostream>
#include <cmath>
#include <fstream>
#include "computers/dynamics.h"
#include "models/V-755DSU/V-755DSUModel.h"
#include "models/V-755DSU/V755Autopilot.h"
#include <chrono>
using namespace std;

struct object{
    double CMD[3]{}, model[17]{};
};

int main(){
    auto start = chrono::system_clock::now();
    V755DSU missile1;
    V755Autopilot autopilot755;
    ofstream state_Missile("MATLAB/states_Missile.txt");
    object interceptor1;
    int stepCounter{0};
    const double dt{0.001}, simEnd{60.0}, deg_rad{PI/180.0}, ft_m{0.304878};
    double time{0}, initialState[15] = {0 /*u*/, 0 /*v*/, 0 /*w*/, 
                    0*deg_rad /*p*/, 0*deg_rad /*q*/, 0*deg_rad /*r*/,
                    0*deg_rad /*phi*/, 48*deg_rad /*theta*/, 0*deg_rad /*psi*/, 
                    0 /*x*/, 0 /*y*/, 0 /*z*/, 
                    0 /*deltaAileron*/, 0 /*deltaElevator*/, 0/*deltaRudder*/};
    for(int i = 0; i < 17; ++i){interceptor1.model[i] = missile1.transferData(i);}
    Dynamics testMissile(Dynamics::Atmosphere::ON, interceptor1.model);
    testMissile.initialConditions(initialState);
    //Initial Data Points
    state_Missile << time << "\t"; 
        for(int i = 0; i < 15; i++){state_Missile << testMissile.getSpecificInstance(i) << " ";}
        for(int i = 0; i < 8; i++){state_Missile << testMissile.getAirData(i) << " ";}
        state_Missile << endl;
    while(stepCounter < simEnd/dt){       
        //if(stepCounter % 100 == 0){cout << "Fb[z]: " << testMissile.getQuickData() << " at " << time << " sec" << endl;}
        time += dt;
        //Simulation Sequence
        testMissile.iterate(dt);
            //Control sequence
        double apInput[10] = {testMissile.getSpecificInstance(3), testMissile.getSpecificInstance(4), testMissile.getSpecificInstance(5),
            testMissile.getDerivatives(3), testMissile.getDerivatives(4), testMissile.getDerivatives(5), 
            testMissile.getMeasuredAccel(Dynamics::y), testMissile.getMeasuredAccel(Dynamics::z), testMissile.getAirData(0), testMissile.getSpecificInstance(6)};
        /*Desired control inputs go here*/ 
        autopilot755.autopilotGenerate(apInput, dt, time);
        for(int i = 0; i < 3; ++i){interceptor1.CMD[i] = autopilot755.autopilotCommands(i);}
        testMissile.steeringCommand(interceptor1.CMD);
            //Update missile based on TOF
        missile1.updateMissile(time, dt);
        for(int i = 0; i < 17; ++i){interceptor1.model[i] = missile1.transferData(i);}
        testMissile.updateObject(interceptor1.model);
        //Data Retrieval
        state_Missile << time << "\t";
        for(int i = 0; i < 15; i++){state_Missile << testMissile.getSpecificInstance(i) << " ";}
        for(int i = 0; i < 8; i++){state_Missile << testMissile.getAirData(i) << " ";}
        state_Missile << endl;        
        if((-testMissile.getSpecificInstance(11) < 0) && time > 1){cout << "[Missile-Ground Impact]" << endl; break;}
        if(testMissile.getAirData(1) > 5  && time > 1){cout << "[Missile Exceeded Mach 5]" << endl; break;}
        stepCounter++;
    }
    auto end = chrono::system_clock::now();
    chrono::duration<double> runtime = end - start;
    cout << "Real Runtime: " << runtime.count() << " seconds" << endl;
    cout << "Sim Runtime: " << time << " seconds" << endl;
    system("PAUSE");
}