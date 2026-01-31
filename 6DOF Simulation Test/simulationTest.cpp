#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <chrono>
#include "computers/MissileDynamics copy.h"
#include "models/V-755DSU/V-755DSUModel.h"
#include "models/V-755DSU/AP755U_anglesOnly.h"
#include "computers/CoordinateDeterminationSystem.h"
#include "computers/CommandGenerationSystem.h"
#include "models/TargetDynamics.h"
struct object{
    double CMD[3]{}, model[20]{};
}; 

int main(){
    auto start = chrono::system_clock::now();
    CommandGenerationSystem P16_1;
    CoordinateDeterminationSystem track, launch1;
    Target drone(0.0, 50.0, 180.0, 400.0, 18.0);
    V755DSU missile1;
    AP755U autopilot755;
    ofstream state_Missile("MATLAB/states_Missile.txt");
    object interceptor1;
    bool trackreset = true, missilereset = true;
    int stepCounter{0};
    const double dt{0.001}, simEnd{60.0}, deg_rad{M_PI/180.0}, ft_m{0.304878};
    double time{0.0}, initialState[15] = {0.0 /*u*/, 0.0 /*v*/, 0.0 /*w*/, 
                    0.0*deg_rad /*p*/, 0.0*deg_rad /*q*/, 0.0*deg_rad /*r*/,
                    0.0*deg_rad /*phi*/, 30.0*deg_rad /*theta*/, 0.0*deg_rad /*psi*/, 
                    0.0 /*x*/, 0.0 /*y*/, 0.0 /*z*/, 
                    0.0 /*deltaAileron*/, 0.0 /*deltaElevator*/, 0.0 /*deltaRudder*/};
    for(int i = 0; i < 20; ++i){interceptor1.model[i] = missile1.transferData(i);}
    MissileDynamics missile(MissileDynamics::Atmosphere::ON, interceptor1.model);
    missile.initialConditions(initialState);
    //Initial Data Points
    state_Missile << time << "\t"; 
        for(int i = 0; i < 15; i++){state_Missile << missile.getInstance(i) << " ";}
        for(int i = 0; i < 7; i++){state_Missile << missile.getAirData(i) << " ";}
        state_Missile << drone.getX() << " " << drone.getY() << " " << drone.getZ() << " ";
        for(int i = 0; i <3; i++){state_Missile << autopilot755.getData(i) << " ";}
        state_Missile << endl;
        double command{0.0};
    while(stepCounter < simEnd/dt){
        time += dt;

        if(trackreset == true){track.initial(drone.getX(), drone.getY(), -drone.getZ()); trackreset = false;}
        drone.move(dt);
        if(stepCounter % 65 == 0){track.final(drone.getX(), drone.getY(), -drone.getZ()); trackreset = true;}
        //track.initial(25000.0, 0.0, -10000.0);
        //track.final(25000.0, 0.0, -10000.0);
        track.rates(dt);
            //Control sequence
        double apInput[11] = {
            missile.getInstance(6), //Phi 
            missile.getInstance(4), //Q
            missile.getInstance(5), //R
            missile.getAirData(8),          //Dynamic Pressure N/m^2
            missile.getMeasuredAccel(MissileDynamics::y), //Gs in OY axis
            missile.getMeasuredAccel(MissileDynamics::z), //Gs in OZ axis
            missile.getInstance(7),    //Theta 
            missile.getInstance(8),    //Psi
            missile.getInstance(12),    //Aileron Angle
            missile.getInstance(13),    //Elevator Angle
            missile.getInstance(14)     //Rudder Angle
        };

        autopilot755.setAutopilot(P16_1.getRangeAmp(),  P16_1.getK1(), P16_1.getK2());
        
        autopilot755.autopilotGenerate(apInput, dt, time);
        for(int i = 0; i < 3; ++i){interceptor1.CMD[i] = autopilot755.getAPCommands(i);}   
        if(time > 0.5){missile.steeringCommand(interceptor1.CMD);}

        if(missilereset == true){launch1.initial(missile.getInstance(9), missile.getInstance(10), missile.getInstance(11)); missilereset = false;}
        missile.iterate(dt);
        if(stepCounter % 65 == 0){launch1.final(missile.getInstance(9), missile.getInstance(10), missile.getInstance(11)); missilereset = true;}
        launch1.rates(dt);

        if(stepCounter % 22 == 0){P16_1.generateSignal(time, launch1.getE(), track.getE(), launch1.getB(), track.getB(), track.getEdot(), track.getBdot());}
        //Update missile based on TOF
        missile1.updateMissile(time, dt);
        for(int i = 0; i < 17; ++i){interceptor1.model[i] = missile1.transferData(i);}
        missile.updateObject(interceptor1.model);

        //Data Retrieval
        state_Missile << time << "\t";
        for(int i = 0; i < 15; i++){state_Missile << missile.getInstance(i) << " ";}
        for(int i = 0; i < 7; i++){state_Missile << missile.getAirData(i) << " ";}
        state_Missile << drone.getX() << " " << drone.getY() << " " << drone.getZ() << " ";
        for(int i = 0; i <3; i++){state_Missile << autopilot755.getData(i) << " ";}
        state_Missile << endl;
        if(launch1.getRange() > track.getRange()){break;}
        if((-missile.getInstance(11) < 0) && time > 1){cout << "[Missile-Ground Impact]" << endl; break;}
        if(missile.getAirData(1) > 5  && time > 1){cout << "[Missile Exceeded Mach 5]" << endl; break;}
        stepCounter++;
    }
    auto end = chrono::system_clock::now();
    chrono::duration<double> runtime = end - start;
    cout << "Real Runtime: " << runtime.count() << " seconds" << endl;
    cout << "Sim Runtime: " << time << " seconds" << endl;
    system("PAUSE");
}