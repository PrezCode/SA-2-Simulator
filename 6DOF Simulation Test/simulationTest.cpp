#include <iostream>
#include <cmath>
#include <fstream>
#include "computers/dynamics.h"
#include "models/V-755DSU/V-755DSUModel.h"
#include "models/V-755DSU/V755Autopilot.h"
#include "computers/CommandGenerationSystem.h"
#include "computers/CoordinateDeterminationSystem.h"
#include <chrono>
using namespace std;


int main(){
    V755DSU missile;
    V755Autopilot autopilot755;
    CommandGenerationSystem signals;
    CoordinateDeterminationSystem track(0), beacon(0);
    ofstream state_Missile("MATLAB/states_Missile.txt");
    auto start = chrono::system_clock::now();
    int stepCounter{0};
    const double dt{0.001}, simEnd{60}, deg_rad{PI/180}, ft_m{0.304878};
    double time{0}, missileData[17], targetData[17], steerCMD[3], aeroData[11],
        initialStateMissile[15] = {0 /*u*/, 0 /*v*/, 0 /*w*/, 
                            0*deg_rad /*p*/, 0*deg_rad /*q*/, 0*deg_rad /*r*/,
                            0*deg_rad /*phi*/, 60*deg_rad /*theta*/, 0*deg_rad /*psi*/, 
                            0 /*x*/, 0 /*y*/, 0 /*z*/, 
                            0 /*deltaAileron*/, 0 /*deltaElevator*/, 0/*deltaRudder*/};
    for(int i = 0; i < 17; ++i){missileData[i] = missile.transferData(i);}
    Dynamics testMissile(Dynamics::Atmosphere::ON, missileData);
    testMissile.initialConditions(initialStateMissile);
    //Initial Data Points
    state_Missile << time << "\t"; 
        for(int i = 0; i < 15; i++){state_Missile << testMissile.getSpecificInstance(i) << " ";}
        for(int i = 0; i < 8; i++){state_Missile << testMissile.getAirData(i) << " ";}
        state_Missile << endl;
    
    while(stepCounter < simEnd/dt){       
        time += dt;
        //Simulation Sequence
        track.initial(25000, 0, 10000);
        track.final(25000, 0, 10000);
        track.rates(dt);

        beacon.initial(testMissile.getSpecificInstance(9), testMissile.getSpecificInstance(10), testMissile.getSpecificInstance(11));
        testMissile.iterate(dt);
        beacon.final(testMissile.getSpecificInstance(9), testMissile.getSpecificInstance(10), testMissile.getSpecificInstance(11));
        beacon.rates(dt);
        signals.generateSignal(time, beacon.getE(), track.getE(), beacon.getB(), track.getB(), track.getEdot(), track.getBdot());
        /*
        if(time > 7 && stepCounter % 100 == 0){
            autopilot755.autopilotGenerate(
            testMissile.getAirData(8), testMissile.getAirData(0), testMissile.getSpecificInstance(4), 
            testMissile.getSpecificInstance(5), testMissile.getDerivatives(10), testMissile.getDerivatives(11), 
            signals.getK1(), signals.getK2(), dt);
            for(int i = 0; i < 3; ++i){steerCMD[i] = autopilot755.autopilotCommands(i);}
            testMissile.steeringCommand(steerCMD);
        }
        */
        missile.updateMissile(time, dt);
        for(int i = 0; i < 17; ++i){missileData[i] = missile.transferData(i);}
        testMissile.updateObject(missileData);
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