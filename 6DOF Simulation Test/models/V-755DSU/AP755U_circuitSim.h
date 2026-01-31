#pragma once
#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
struct Commands{double OX0_deg, OX_deg, OY0_deg, OZ0_deg, OX_rad, OY_rad, OZ_rad, OX_pqfb_deg, OY_signal, OZ_signal;};
struct Gains{double dg, ddg, du, ddu, idu, sg, dsg, kpt, rkpt;};
struct Axis{double measured{0.0}, initial{0.0}, dError{0.0}, Error0{0.0}, Error{0.0}, phi_deg{0.0},
            wRate0{0.0}, wRateFiltered{0.0}, wAccel{0.0}, wDamp{0.0}, g{0.0}, gAccel{0.0}, g0{0.0}, gDelta{0.0}, feedback_deg{0.0}, feedbackPq{0.0};};
enum Subscripts{A = 0, E = 1, R = 2};
class AP755U{
    public:
    Commands cmd;
    Gains K;
    Axis OX, OY, OZ;
    AP755U(){ 
        //Autopilot Settings
        cmd.OX0_deg = 45.0;   //Force roll angle of 45 deg
        cmd.OY0_deg = 0.0;      //Set OY actuator to 0, conventionally the elevator
        cmd.OZ0_deg = 0.0;      //Set OZ actuator to 0, conventionally the rudder
        OX.measured = cmd.OX0_deg;
        K.dg = 0.7;             //Angular rate gain, deg/deg/sec
        K.du = 13.0;            //Lateral acceleration gain, deg/g
        K.sg = 2.8;             //Initial actuator gain of ailerons, deg/deg
        K.kpt = 100.0;          //Initial gain of actuator track, deg/amp-turn
        K.rkpt = 83.5;          //Gain of Aileron actuator track, deg/amp-turn
        deltaMaxA = 12.0*d2r;   //Aileron Limit
        deltaMaxE = 28.0*d2r;   //Elevator Limit
        deltaMaxR = 28.0*d2r;   //Rudder Limit
    }
    void setAutopilot(double rollcmd, double K1, double K2){
        cmd.OY0_deg = 0.25*(K1/2870)*K.kpt;
        cmd.OZ0_deg = 0.25*(K2/2870)*K.kpt;
        cmd.OX0_deg = rollcmd;
    }
    void autopilotGenerate(double* apInput, double dt, double time){
        OX.measured = apInput[0]*r2d;          //phi in deg
        OY.wRateFiltered = threshold(clamp(apInput[1]*r2d, -45.0, 45.0), 0.3);  //q in deg/s, max measurement 45*/s
        OZ.wRateFiltered = threshold(clamp(apInput[2]*r2d, -45.0, 45.0), 0.3);  //r in deg/s, max measurement 45*/s
        double Pq = apInput[3]*Pa_kgf;                          //Dynamic Pressure converted to kg/cm2
        OY.g = clamp(apInput[4], -8.0, 8.0);                    //Gs in missile vertical axis 
        OZ.g = clamp(apInput[5], -8.0, 8.0);                    //Gs in missile horizontal axis
        double AileronAngle = apInput[6]*r2d;
        double K1angle = apInput[7]*r2d;
        double K2angle = apInput[8]*r2d;

        AC5Sensor(Pq);
        US14Gyro(dt);
        US2Sensor(dt);
        C12FreeGyro(dt, Pq);
        AC34a2(AileronAngle);        
        US17SignalCombinator(dt);
        deltaCMD[A] = clamp(cmd.OX_rad, -deltaMaxA, deltaMaxA);
        if(time > 6.0){
            //deltaCMD[E] = clamp(cmd.OY_rad, -deltaMaxE, deltaMaxE);
            //deltaCMD[R] = clamp(cmd.OZ_rad, -deltaMaxR, deltaMaxR);
        }else{deltaCMD[E] = deltaCMD[R] = 0.0;}
    }
    double getData(int i){
        switch(i){
            case 0: return cmd.OX_deg; break;
            case 1: return OX.feedback_deg; break;
            case 2: return cmd.OX_pqfb_deg; break;
            default: return 0;
        }
    }
    void US14Gyro(double dt){    //Gyro to dampen roll rate in OY and OZ axes
        double tau = 4.0*pow(10,-6)*(10250);        //Time constant of forming cell
        K.ddg = K.dg*tau;                           //Derivative gain
        OY.wAccel = (OY.wRateFiltered - OY.wRate0)/dt;      //Measured acceleration
        OZ.wAccel = (OZ.wRateFiltered - OZ.wRate0)/dt;
        OY.wDamp = forceNeg(K.dg*OY.wRateFiltered + K.ddg*OY.wAccel); //Actuation signal to be subtracted
        OZ.wDamp = forceNeg(K.dg*OZ.wRateFiltered + K.ddg*OZ.wAccel);
        OY.wRate0 = OY.wRateFiltered;                       //Reset initial rate for next iteraiton
        OZ.wRate0 = OZ.wRateFiltered;
    }
    void US2Sensor(double dt){         //Linear Acceleration Sensor to correct acceleration in OY and OZ axes
        double tau = 99.0*pow(10, -6)*(1000.0);   //Time constant of forming cell
        K.ddu = K.du*tau;                       //Derivatve gain
        K.idu = K.du*dt;                        //Integral gain
        OY.gDelta = (OY.g - OY.g0)/dt;          //Change in acceleration
        OZ.gDelta = (OZ.g - OZ.g0)/dt;
        OY.gAccel = K.du*OY.g + K.idu*OY.g*dt + K.ddu*OY.gDelta;    //PID shaped command to add
        OZ.gAccel = K.du*OZ.g + K.idu*OZ.g*dt + K.ddu*OZ.gDelta;
        OY.g0 = OY.g;                           //Reset g measurement for next iteration
        OZ.g0 = OZ.g;
    }
    void C12FreeGyro(double dt, double Pq){  //Bank Angle Sensor
        double PqVoltage = C5_2Sensor(Pq);
        temp = PqVoltage;
        K.sg = 2.8*(PqVoltage/26.0);
        OX.Error = PqVoltage*1.04*clamp(cmd.OX0_deg - OX.measured, -25.0, 25.0);    //Measured roll error signal, max +/-25 deg, 1.04V/deg
        OX.dError = (OX.Error - OX.Error0)/dt; //Change in voltage measured proportional to angular rate
        OX.Error0 = OX.Error;
        cmd.OX_deg = (OX.Error/20510.0 + 0.000004*OX.dError);    //Input current to US-17
        /*
        //double tau = 4.0*pow(10, -6)*(20000.0);         //Time constant of forming cell
        OX.phiError = clamp(cmd.roll0_deg - OX.phi, -25.0, 25.0);
        K.dsg = K.sg*tau;                               //Gain of derivative
        OX.dPhi_degps = (OX.phi0 - OX.phi)/dt;          //Change in roll rate
        cmd.roll_deg = (K.sg*OX.phi_deg + K.dsg*OX.dPhi_degps); //Command roll amount
        OX.phi0 = OX.phi;
        */
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
    void AC34a2(double delE){    //Aileron motion feedback system
        const double Rload = 5100.0;
        const double Vsource = 26.0;
        double delEchange = (delE0 - delE);
        delE0 = delE;
        cmd.OX_pqfb_deg = Vsource*(delEchange/13.0)/Rload;
    }
    void AC5Sensor(double Pq){                      //Dynamic Pressure Sensor and gain control
        //Gains change hyperbolically a function of Pq 
        //K(Pq) = A/Pq + B
        //A = (Kmin - Kmax)/(Pmin^-1 - Pmax^-1), Kmin is K at Pmin
        //B = Kmin - A/Pmin
        if(Pq >= 1.6 && Pq <= 5.1){K.kpt = 160.16914/Pq + 0.105714;}
            else{K.kpt = 100.0;}
    }
    void US17SignalCombinator(double dt){    //Combines Signals to generate stability feedback and final actuator commands
        double formingCell = 486.0/3500.0;  //Approx 13.8% of command
        
        cmd.OX_rad = (K.rkpt*threshold(0.25*(cmd.OX_deg - OX.feedback_deg - cmd.OX_pqfb_deg), 0.007))*d2r;
        OX.feedback_deg = FeedbackLoop(cmd.OX_rad*r2d/K.sg/K.rkpt, A, dt);
        //cmd.OY_rad = (clamp(cmd.ch1_deg + OY.gAccel, -30.0, 30.0) + OY.wDamp  - OY.feedback)*d2r; 
        //cmd.OZ_rad = (clamp(cmd.ch2_deg + OZ.gAccel, -30.0, 30.0) + OZ.wDamp  - OZ.feedback)*d2r;
        //OY.feedback = (clamp(cmd.ch1_deg + OZ.gAccel, -30.0, 30.0) + OY.wDamp)*formingCell;
        //OZ.feedback = (clamp(cmd.ch2_deg + OY.gAccel, -30.0, 30.0) + OZ.wDamp)*formingCell;
        
    }
    double FeedbackLoop(double commandCurrent, Subscripts AER, double dt){
        static double cap = 0.000004;
        static double Rfc = 3500.0;
        static double Rsense = 486.0;
        static double Iout0 = 0.007*Rsense;
        Vfb[AER][1] = commandCurrent*Rsense + Iout0;
        double Vdelta = (Vfb[AER][1] - Vfb[AER][0])/dt;
        Vfb[AER][0] = Vfb[AER][1];
        double feedbackCurrent = clamp(commandCurrent*(Rsense/Rfc) + cap*Vdelta, -commandCurrent, commandCurrent);
        return feedbackCurrent;
    }
    double getAPCommands(int i){return deltaCMD[i];}
    private:
    const double  g{9.80665};
    double deltaCMD[3]{0,0,0}, d2r{M_PI/180.0}, r2d{180.0/M_PI}, Pa_kgf{1.0/98066.5}, Vfb[3][2]{0.007*486.0}, deltaMaxA, deltaMaxE, deltaMaxR, delE0{0.0}, rollFeedback{0.0}, temp{0.0};
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