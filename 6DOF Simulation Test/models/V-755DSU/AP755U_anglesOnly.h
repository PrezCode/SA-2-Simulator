#pragma once
#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
struct Commands{double OX0_deg, OX_deg, OY0_deg, OY_deg , OZ0_deg, OZ_deg, OX_rad, OY_rad, OZ_rad, OX_pqfb_deg;};
struct Gains{double dg, ddg, du, ddu, idu, sg, dsg, kpt, dkpt;};
struct Axis{double measured{0.0}, initial{0.0}, dError{0.0}, Error0{0.0}, Error{0.0}, phi_deg{0.0}, wRateMeasured{0.0},
            wRate0{0.0}, wRateFiltered{0.0}, wAccel{0.0}, wDamp{0.0}, g{0.0}, gAccel{0.0}, g0{0.0}, gDelta{0.0}, feedback_deg{0.0}, feedbackPq{0.0};};
enum Subscripts{A = 0, E = 1, R = 2};
class AP755U{
    public:
    Commands cmd;
    Gains K;
    Axis OX, OY, OZ;
    AP755U(){ 
        //Autopilot Settings
        cmd.OX0_deg = 0.0;   //Force roll angle of 45 deg
        cmd.OY0_deg = 0.0;      //Set OY actuator to 0, conventionally the elevator
        cmd.OZ0_deg = 0.0;      //Set OZ actuator to 0, conventionally the rudder
        OX.measured = cmd.OX0_deg;
        K.dg = 0.7;             //Angular rate gain, deg/deg/sec
        K.du = 13.0;            //Lateral acceleration gain, deg/g
        K.sg = 2.8;             //Initial actuator gain of ailerons, command deg per measured deg
        K.kpt = 3.3533;            //Initial actuator gain of rudders, command deg per measured deg
        deltaMaxA = 12.0*d2r;   //Aileron Limit
        deltaMaxE = 28.0*d2r;   //Elevator Limit
        deltaMaxR = 28.0*d2r;   //Rudder Limit
    }
    void setAutopilot(double range, double pitchcmd, double yawcmd){
        if(pitchcmd < 0){cmd.OY0_deg = (28.0/(147 + range))*(pitchcmd - range);}
            else{cmd.OY0_deg = (28.0/(147 + range))*(pitchcmd + range);}
        //if(pitchcmd == pitchcmd0){pitchcmd = 0.0;}
        //pitchcmd0 = pitchcmd;
        if(yawcmd < 0){cmd.OZ0_deg = -(28.0/(93 + range))*(yawcmd - range);}
            else{cmd.OZ0_deg = -(28.0/(93 + range))*(yawcmd + range);}
        //if(yawcmd == yawcmd0){yawcmd = 0.0;}
        //yawcmd0 = yawcmd;
    }
    void autopilotGenerate(double* apInput, double dt, double time){
        OX.measured = apInput[0]*r2d;          //phi in deg
        OY.wRateMeasured = threshold(clamp(apInput[1]*r2d, -45.0, 45.0), 0.3);  //q in deg/s, max measurement 45*/s
        OZ.wRateMeasured = threshold(clamp(apInput[2]*r2d, -45.0, 45.0), 0.3);  //r in deg/s, max measurement 45*/s
        double Pq = apInput[3]*Pa_kgf;                          //Dynamic Pressure converted to kg/cm2
        OY.g = clamp(apInput[4], -8.0, 8.0);                    //Gs in missile vertical axis 
        OZ.g = clamp(apInput[5], -8.0, 8.0);                    //Gs in missile horizontal axis
        OY.measured = apInput[6]*r2d;   //theta in deg
        OZ.measured = apInput[7]*r2d;   //psi in deg
        del[A] = apInput[8]*r2d;   //Measured aileron angle
        del[E] = apInput[9]*r2d;        //Measured elevator angle
        del[R] = apInput[10]*r2d;       //Measured rudder angle

        US14Gyro(dt);
        US2Sensor(dt);
        C12FreeGyro(dt, Pq);
        US17SignalCombinator(dt, Pq);
        deltaCMD[A] = clamp(cmd.OX_rad, -deltaMaxA, deltaMaxA);
        if(time > 6.0){
            deltaCMD[E] = clamp(cmd.OY_rad, -deltaMaxE, deltaMaxE);
            deltaCMD[R] = clamp(cmd.OZ_rad, -deltaMaxR, deltaMaxR);
        }else{deltaCMD[E] = deltaCMD[R] = 0.0;}
    }
    double getData(int i){
        switch(i){
            case 0: return cmd.OY_deg; break;
            case 1: return OY.gAccel; break;
            case 2: return OY.wDamp; break;
            default: return 0;
        }
    }
    void US14Gyro(double dt){    //Gyro to dampen roll rate in OY and OZ axes
        double tau = (4.0E-6)*7700.0;        //Time constant of forming cell
        const double J = 3.4E-6;
        const double b = 4.4E-4;
        const double k = 0.018;
        K.ddg = K.dg*tau;                           //Derivative gain
        OY.wRateFiltered = J*OY.wAccel + b*OY.wRateMeasured + k*OY.wRateMeasured*dt;
        OZ.wRateFiltered = J*OZ.wAccel + b*OZ.wRateMeasured + k*OZ.wRateMeasured*dt;
        OY.wAccel = (OY.wRateFiltered - OY.wRate0)/dt;      //Measured acceleration
        OZ.wAccel = (OZ.wRateFiltered - OZ.wRate0)/dt;
        OY.wDamp = K.dg*OY.wRateMeasured + K.ddg*OY.wAccel; //Actuation signal to be subtracted
        OZ.wDamp = K.dg*OZ.wRateMeasured + K.ddg*OZ.wAccel;
        OY.wRate0 = OY.wRateFiltered;                       //Reset initial rate for next iteraiton
        OZ.wRate0 = OZ.wRateFiltered;
    }
    void US2Sensor(double dt){         //Linear Acceleration Sensor to correct acceleration in OY and OZ axes
        double tau = (99.0E-6)*(22600.0);   //Time constant of forming cell
        K.idu = K.du/tau;                  //Integral gain
        OY.gAccel = K.du*OY.g + K.idu*OY.g*dt;    //PID shaped command to add
        OZ.gAccel = K.du*OZ.g + K.idu*OZ.g*dt;
        OY.g0 = OY.g;                           //Reset g measurement for next iteration
        OZ.g0 = OZ.g;
    }
    double D7Integrator(Subscripts ER, double dt, double Pq){
        //double tau = (99.0E-6)*(22600.0);
        if(Pq >= 1.6){K.kpt = 5.36567/Pq + 0.0035414;}
            else{K.kpt = 3.3533;}
        //K.dkpt = K.kpt*tau;
        switch(ER){
            case E:
                cmd.OY_deg = K.kpt*cmd.OY0_deg;
                /*OY.Error = cmd.OY0_deg - OY.measured;
                OY.dError = (OY.measured - OY.initial)/dt;
                OY.initial = OY.measured;
                cmd.OY_deg = (K.kpt*OY.Error + K.dkpt*OY.dError);*/
                return clamp(cmd.OY_deg + OY.gAccel, -30.0, 30.0); break;
            case R:
                cmd.OZ_deg = K.kpt*cmd.OZ0_deg;
                /*OZ.Error = cmd.OZ0_deg - OZ.measured;
                OZ.dError = (OZ.measured - OZ.initial)/dt;
                OZ.initial = OZ.measured;
                cmd.OZ_deg = (K.kpt*OZ.Error + K.dkpt*OZ.dError);
                cmd.OZ_deg = cmd.OZ0_deg;*/
                return clamp(cmd.OZ_deg + OZ.gAccel, -30.0, 30.0); break;
            default: return 0;
        }
    }
    void C12FreeGyro(double dt, double Pq){  //Bank Angle Sensor
        double PqVoltage = C5_2Sensor(Pq);
        double tau = (4.0E-6)*20000.0;         //Time constant of forming cell
        K.sg = 2.8*(PqVoltage/26.0);
        K.dsg = K.sg*tau;                               //Gain of derivative
        OX.Error = clamp((cmd.OX0_deg - OX.measured), -25.0, 25.0);
        OX.dError = (OX.initial - OX.measured)/dt;          //Change in roll rate
        cmd.OX_deg = (K.sg*OX.Error + K.dsg*OX.dError); //Command roll amount
        OX.initial = OX.measured;
    }
    double C5_2Sensor(double Pq){   //Automatically sets voltage delivered to C12FreeGyro based on airspeed
        const double Rpotmax = 1050.0;  //C-5-2 pot
        const double percentMin = 0.01;
        const double percentMax = 0.997;
        const double RIII_14 = 2500.0;   //Pot for trim, max 5kOhms
        const double RIII_7 = 470.0;     //Constant resistor for trim
        const double RIII_18 = 124.0;    //Load resistor
        const double Vsource = 26.0;
        const double Pqmin = 0.3;
        const double Pqmax = 2.3;
        double alpha;
        if(Pq < 0.3){alpha = 0.01;}
            else if(Pq >= 0.3 && Pq <= 2.3){alpha = percentMin + (percentMax - percentMin)*((Pq - Pqmin)/(Pqmax - Pqmin));}
            else{alpha = 0.997;}
        double Rtop = alpha*Rpotmax; 
        double Rbot = (1.0 - alpha)*Rpotmax; 
        double Rdown = 1.0/(1.0/(Rbot + RIII_18) + 1.0/(RIII_14 + RIII_7));
        double Vdelivered = Vsource*(Rdown/(Rtop + Rdown));
        return Vdelivered;
    }
    void US17SignalCombinator(double dt, double Pq){    //Combines Signals to generate stability feedback and final actuator commands
        double formingCell = 486.0/3500.0;  //Approx 13.8% of command
        cmd.OX_rad = (cmd.OX_deg - OX.feedback_deg - OX.feedbackPq)*d2r;
        OX.feedback_deg = cmd.OX_rad*formingCell*r2d;
        OX.feedbackPq = K.sg*del[A];
        
        cmd.OY_rad = (D7Integrator(E, dt, Pq) - OY.wDamp  - OY.feedback_deg - OY.feedbackPq)*d2r; 
        OY.feedback_deg = cmd.OY_rad*formingCell*r2d;
        OY.feedbackPq = K.kpt*del[E];
        cmd.OZ_rad = (D7Integrator(R, dt, Pq) - OZ.wDamp  - OZ.feedback_deg - OZ.feedbackPq)*d2r; 
        OZ.feedback_deg = cmd.OZ_rad*formingCell*r2d;
        OZ.feedbackPq = K.kpt*del[R];
    }
    double getAPCommands(int i){return deltaCMD[i];}
    private:
    const double  g{9.80665};
    double deltaCMD[3]{0,0,0}, d2r{M_PI/180.0}, r2d{180.0/M_PI}, Pa_kgf{1.0/98066.5}, deltaMaxA, deltaMaxE, deltaMaxR, temp{0.0}, del[3]{}, pitchcmd0{0.0}, yawcmd0{0.0};
    static double clamp(double x, double low, double hi) {return (x < low) ? low : (x > hi) ? hi : x;}
    static double threshold(double x, double threshold){return (x > threshold) ? x : (x < -threshold) ? x : 0;}
    static double forceNeg(double x){
        if(x > 0){return -x;}
        return x;
    }
    static double wrapToPi(double a){
        while(a >  M_PI){a -= 2.0*M_PI;}
        while(a < -M_PI){a += 2.0*M_PI;}
        return a;
    }
};