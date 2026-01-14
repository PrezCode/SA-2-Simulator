#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>

struct Commands{double pdot_rps, qdot_rps, rdot_rps, p, q, r, nz_g, ny_g, phi_rad;};
struct Gains{double az0, ay0, pq, dq, pr, dr, az, ay, iz, iy, pphi, iphi, pp, dp;};
enum Subscripts{A = 0, E = 1, R = 2};
class V755Autopilot{
    public:
    Commands cmd;
    Gains K;
    V755Autopilot(){
        //Autopilot Settings
        t_inhibit = 6.0;
        t_ramp = 2.0;
        tau = 0.1; 
        Wq = Wr = 1.0/(4.0*tau);
        Waz = Way = Wq/3.0;
        cmd.ny_g = 0.0;
        cmd.nz_g = 0.0; 
        cmd.phi_rad = 0.0*deg_rad;
        K.ay0 = 0.3/g;
        K.az0 = 0.3/g;
        K.iz = 0.0;
        K.iy = 0.0;
        K.pq = 2.5;
        K.pr = 2.5;
        K.dq = 0.2;
        K.dr = 0.2;
        K.pphi = 2.0;
        K.iphi = 0.0;
        K.pp = 2.0;
        K.dp = 0.2;
        V0_mps = 750.0;
        Vmin_mps = 200.0;
        Von = 250.0;
        Vfull = 400.0;
        deltaMaxA = 12.0/180.0*PI;
        deltaMaxE = 28.0/180.0*PI;
        deltaMaxR = 28.0/180.0*PI;
    }
    void setAutopilot(double ny_g, double nz_g){
        cmd.ny_g = ny_g;
        cmd.nz_g = nz_g;
    }
    void autopilotGenerate(double* apInput, double dt, double time){
        double p = apInput[0];
        double q = apInput[1];
        double r = apInput[2];
        double pdot = apInput[3];
        double qdot = apInput[4];
        double rdot = apInput[5];
        double ny_meas_g = apInput[6]/g;
        double nz_meas_g = -apInput[7]/g;
        double Vmag_mps = apInput[8];
        double phi_rad = apInput[9];
        S = clamp((time - t_inhibit)/t_ramp, 0, 1);
        if(time < t_inhibit){Et = 0.0;}   
            else if((t_inhibit <= time) && (time <= t_inhibit + t_ramp)){Et = S*S*(3.0 - 2.0*S);}
            else{Et = 1.0;}
        Ev = calculateEv(Vmag_mps);
        Efactor = Et*Ev;
        gainSchedule(Vmag_mps);
        if(time < t_inhibit || Efactor <= 0.0){
            iaz = iay = 0.0;
            cmd.q = cmd.r = 0.0;
            cmd.qdot_rps = cmd.rdot_rps = 0.0;
            deltaCMD[A] = deltaCMD[E] = deltaCMD[R] = 0.0;
            return;
        }
        rollAngleCorrection(phi_rad, dt, p, pdot);
        accelerationError(time, dt, ny_meas_g, nz_meas_g, phi_rad);
        rateTracking(p, q, r);
        computeCommands(p, q, r, pdot, qdot, rdot);
        deltaCMD[A] = clamp(deltaCMD[A], -deltaMaxA, deltaMaxA);
        deltaCMD[E] = clamp(deltaCMD[E], -deltaMaxE, deltaMaxE);
        deltaCMD[R] = clamp(deltaCMD[R], -deltaMaxR, deltaMaxR);
    }
    void rollAngleCorrection(double phi_rad, double dt, double p, double pdot){
        double phiError = wrapToPi(cmd.phi_rad - phi_rad);
        iphi = clamp(iphi + Efactor*phiError*dt, -iphiMax, iphiMax);
        cmd.p = K.pphi*phiError + K.iphi*iphi;
        cmd.pdot_rps = Wr*(cmd.p - p);
        deltaCMD[A] = (K.pp*(cmd.p - p) + K.dp*(cmd.pdot_rps - pdot))*Efactor;
    }
    void accelerationError(double time, double dt, double ny_meas_g, double nz_meas_g, double phi_rad){
        iaz = clamp(iaz + Efactor*(cmd.nz_g - nz_meas_g)*dt, -iazMax, iazMax);
        iay = clamp(iay + Efactor*(cmd.ny_g - ny_meas_g)*dt, -iayMax, iayMax);
        cmd.q = K.az*(cmd.nz_g - nz_meas_g) + K.iz*iaz;
        cmd.r = K.ay*(cmd.ny_g - ny_meas_g) + K.iy*iay;
    }
    void computeCommands(double p, double q, double r, double pdot, double qdot, double rdot){
        deltaCMD[E] = (K.pq*(cmd.q - q) + K.dq*(cmd.qdot_rps - qdot))*Efactor;
        deltaCMD[R] = (K.pr*(cmd.r - r) + K.dr*(cmd.rdot_rps - rdot))*Efactor;
    }
    void rateTracking(double p, double q, double r){
        cmd.qdot_rps = Wq*(cmd.q - q);
        cmd.rdot_rps = Wr*(cmd.r - r);
    }
    void gainSchedule(double Vmag_mps){
        const double denom = std::max(Vmag_mps, Vmin_mps);
        const double product = clamp(V0_mps/denom, 0.5, 3.0);
        K.az = K.az0*product;
        K.ay = K.ay0*product;
    }
    double autopilotCommands(int i){
        return deltaCMD[i];
    }
    private:
    double deltaCMD[3]{0,0,0}, tau, deg_rad{PI/180}, g{9.81}, Wq{0}, Wr{0}, Wp{0}, Waz{0}, Way{0}, 
    Efactor{0.0}, Et{0.0}, Ev{0.0}, Vmin_mps, V0_mps, iaz{0.0}, iay{0.0}, iazMax{5.0}, iayMax{5.0}, iphi{0.0}, iphiMax{0.5},
    S, t_inhibit, t_ramp, Von, Vfull, deltaMaxA, deltaMaxE, deltaMaxR;
    static double clamp(double x, double lo, double hi) {
        return (x < lo) ? lo : (x > hi) ? hi : x;
    }
    static double wrapToPi(double a){
        while(a >  PI){a -= 2.0*PI;}
        while(a < -PI){a += 2.0*PI;}
        return a;
    }
    double calculateEv(double Vmag){
        if(Vmag <= Von){return 0.0;}
        if(Vmag >= Vfull){return 1.0;}
        return (Vmag - Von)/(Vfull - Von);
    }
};