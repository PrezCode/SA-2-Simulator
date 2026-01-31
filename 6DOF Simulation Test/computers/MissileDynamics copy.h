//6-DOF Physics Simulation for Missiles
#pragma once
#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
class MissileDynamics{
    public:
    enum ModelType{Missile, Aircraft};
    enum Atmosphere{ON, OFF};
    enum Variables{x = 0, y = 1, z = 2, xz = 3, u = 0, v = 1, w = 2, p = 0, q = 1, r = 2, L = 0, M = 1, N = 2};
    enum DataType{Mass = 0, Ix = 1, Iy = 2, Iz = 3, Ixz = 4, Length = 5, Radius = 6, Span = 7, Chord = 8, AileronARef = 9, 
        WingARef = 10, RudderARef = 11, RudderTau = 12, Thrust = 13, CDrag = 14, CLift = 15, CGCP = 16, AileronTau = 17, DelAMax = 18, DelERMax = 19};
    MissileDynamics(Atmosphere setting, double* modelData){
        if(setting == ON){atmosphereON = true;}
            else if(setting == OFF){atmosphereON = false;}
        mass = modelData[Mass];
        I[x] = modelData[Ix];
        I[y] = modelData[Iy];
        I[z] = modelData[Iz];
        I[xz] = modelData[Ixz];
        length = modelData[Length];
        radius = modelData[Radius];
        b = modelData[Span];
        c = modelData[Chord];
        Aref_aileron = modelData[AileronARef];
        Aref_wing = modelData[WingARef];
        Aref_rudder = modelData[RudderARef];
        tau_rudder = modelData[RudderTau];
        tau_aileron = modelData[AileronTau];
        thrust = modelData[Thrust];
        C_DragBase = modelData[CDrag];
        volume = M_PI*pow(radius, 2)*length;
        Aref = M_PI*pow(radius, 2);
        C_LiftBase = modelData[CLift];
        CGCPavg = modelData[CGCP];
        deltaAmax = modelData[DelAMax];
        deltaERmax = modelData[DelERMax];
        diameter = 2.0*radius;
    }
    void updateObject(double* modelData){
        thrust = modelData[Thrust];
        mass = modelData[Mass];
        I[x] = modelData[Ix];
        I[y] = modelData[Iy];
        I[z] = modelData[Iz];
        I[xz] = modelData[Ixz];
        length = modelData[Length];
        radius = modelData[Radius];
        volume = M_PI*pow(radius, 2)*length;
        Aref = M_PI*pow(radius, 2);
        diameter = 2.0*radius;
    }
    void initialConditions(double* state0){
        for(int i = 0; i < 3; i++){
            instance[i] = V[i] = state0[i];
            instance[i+3] = W[i] = state0[i+3];
            instance[i+6] = eulerAngles[i] = state0[i+6];
            instance[i+9] = position[i] = state0[i+9];
            instance[i+12] = deltaCMD[i] = state0[i+12];
        }
        initializeDCM();
        computeAirData();
        estimateCoeffs();
    }
    void steeringCommand(double* command){
        deltaCMD[p] = command[p];
        deltaCMD[q] = command[q];
        deltaCMD[r] = command[r];
    }
    void externalForces(){
        double dynamicPressure = qBar();
        estimateCoeffs();
        dragForce = C_Drag*dynamicPressure*Aref;
        liftForce = C_LiftPitch*dynamicPressure*Aref_wing;
        sideForce = C_LiftYaw*dynamicPressure*Aref_wing;
        Fb[x] = -(Cm_wbT[0][0]*dragForce + Cm_wbT[0][1]*sideForce + Cm_wbT[0][2]*liftForce) + thrust;
        Fb[y] = -(Cm_wbT[1][0]*dragForce + Cm_wbT[1][1]*sideForce + Cm_wbT[1][2]*liftForce);
        Fb[z] = -(Cm_wbT[2][0]*dragForce + Cm_wbT[2][1]*sideForce + Cm_wbT[2][2]*liftForce);
        OZ_meas = (Fb[y]/mass + g*Cm_bn[1][2])/g;   //Felt lateral acceleration in the OZ axis, traditionally the East/Y axis of the body
        OY_meas = (Fb[z]/mass + g*Cm_bn[2][2])/g;   //Felt lateral acceleration in the OY axis, traditionall the Down/Z axis of the body
        moments[L] = MomentCoefficient(L)*dynamicPressure*(diameter + 0.164)*Aref;
        moments[M] = MomentCoefficient(M)*dynamicPressure*length*Aref_wing;
        moments[N] = MomentCoefficient(N)*dynamicPressure*length*Aref_wing;
    }
    double getMeasuredAccel(Variables yz){
        switch(yz){
            case z: return OZ_meas; break;
            case y: return OY_meas; break;
        }
        return 0;
    }
    void derivatives(){
        //Velocity Accelerations
        instanceDeriv[0] = Fb[x]/mass + g*Cm_bn[0][2] - V[w]*W[q] + V[v]*W[r];
        instanceDeriv[1] = Fb[y]/mass + g*Cm_bn[1][2] - V[u]*W[r] + V[w]*W[p];
        instanceDeriv[2] = Fb[z]/mass + g*Cm_bn[2][2] - V[v]*W[p] + V[u]*W[q];
        //Rotation Accelerations
        instanceDeriv[3] = (I[xz]*(I[x] - I[y] + I[z])*W[p]*W[q] - (I[z]*(I[z] - I[y]) + pow(I[xz], 2))*W[q]*W[r] + I[z]*moments[L] + I[xz]*moments[N])/(I[x]*I[z] - pow(I[xz], 2));
        instanceDeriv[4] = ((I[z] - I[x])*W[r]*W[p] - I[xz]*(pow(W[p], 2) - pow(W[r], 2)) + moments[M])/I[y];
        instanceDeriv[5] = (-I[xz]*(I[x] - I[y] + I[z])*W[q]*W[r] + (I[x]*(I[x] - I[y]) + pow(I[xz], 2))*W[p]*W[q] + I[xz]*moments[L] + I[x]*moments[N])/(I[x]*I[z] - pow(I[xz], 2));
        //Euler Angle rates
        instanceDeriv[6] = (W[p] + sin(eulerAngles[p])*tan(eulerAngles[q])*W[q] + cos(eulerAngles[p])*tan(eulerAngles[q])*W[r]);
        instanceDeriv[7] = (cos(eulerAngles[p])*W[q] - sin(eulerAngles[p])*W[r]);
        instanceDeriv[8] = ((sin(eulerAngles[p])/cos(eulerAngles[q]))*W[q] + (cos(eulerAngles[p])/cos(eulerAngles[q]))*W[r]); 
        //XYZ Motion
        instanceDeriv[9]  = Cm_bnT[0][0]*V[u] + Cm_bnT[0][1]*V[v] + Cm_bnT[0][2]*V[w];
        instanceDeriv[10] = Cm_bnT[1][0]*V[u] + Cm_bnT[1][1]*V[v] + Cm_bnT[1][2]*V[w];
        instanceDeriv[11] = Cm_bnT[2][0]*V[u] + Cm_bnT[2][1]*V[v] + Cm_bnT[2][2]*V[w];
        //Actuator Angle Rates
        instanceDeriv[12] = clamp((deltaCMD[0] - instance[12])/tau_aileron, -250.0*d2r, 250.0*d2r);    //aileron
        instanceDeriv[13] = clamp((deltaCMD[1] - instance[13])/tau_rudder, -250.0*d2r, 250.0*d2r);    //elevator
        instanceDeriv[14] = clamp((deltaCMD[2] - instance[14])/tau_rudder, -250.0*d2r, 250.0*d2r);    //rudder
    }
    void computeInstance(){
        computeAirData();
        cosineMatrix(eulerAngles[p], eulerAngles[q], eulerAngles[r], AoA, sideslip);
        externalForces();
        derivatives();
    }
    void iterate(double dt){
        RK4(dt);
        unpackInstance();
        computeAirData();
    }
    double getInstance(int i){return instance[i];}
    double getDerivatives(int i){return instanceDeriv[i];}
    double getAirData(int i){
        switch(i){
            case 0: return Vmag; break;
            case 1: return mach; break;
            case 2: return AoA; break;
            case 3: return sideslip; break;
            case 4: return C_Drag; break;
            case 5: return C_LiftPitch; break;
            case 6: return C_LiftYaw; break;
            case 7: return rho; break;
            case 8: return 0.5*rho*pow(Vmag, 2); break;
            default: return 0; break;
        }
    }
    double getQuickData(){return Cm_bnT[1][0]*Cm_bn[0][1];}
    private:
    bool atmosphereON;
    const double g{9.80665};
    double Cm_bn[3][3], Cm_bnT[3][3],Cm_wb[3][3], Cm_wbT[3][3],   //Cosine matrices for body-normal, normal-body, wind-body
            I[4],       //Inertia, XX, YY, ZZ, XZ
            position[3],//XYZ vector
            Fb[3],      //force vector XYZ
            V[3],       //velocity vectors u, v, w
            W[3],       //rotation vectors p, q, r
            eulerAngles[3],  //roll: phi, pitch: theta, yaw: psi 
            moments[3],    //Moment vectors L, M, N
            deltaCMD[3] = {0,0,0},   //Actuator command values: Aileron, Elevator, Rudder
            instance[15], instanceDeriv[15],  //Arrays to store simulation state and rates for integration
            Aref, Aref_aileron, Aref_wing, Aref_rudder, C_DragBase, C_LiftBase, C_Drag, C_LiftPitch, C_LiftYaw, diameter{0}, tau_aileron{}, deltaAmax{}, deltaERmax{},
            CnBeta{0}, Cmq{0}, Clp{0}, Cnr{0}, CmAoA{0}, CmBeta{-2.0}, CnAoA{-2.0}, volume, CndelR{0}, CmdelE{0}, CldelA{0}, d2r{M_PI/180.0}, C_Nalpha{0}, K_D{0}, OZ_meas{0.0}, OY_meas{0.0},
            mass, AoA, AoAdelE, AoAdelR, sideslip, height, soundBarrier, mach, dragForce, sideForce, liftForce, Vmag{0}, tau_rudder{1}, thrust{0}, length, radius, b, c, CGCPavg, rho;
    //Various calculators
    void computeAirData(){
        Vmag = sqrt(V[u]*V[u] + V[v]*V[v] + V[w]*V[w]);
        if(soundBarrier > 0.0){mach = Vmag/soundBarrier;}
            else{mach = 0.0;}
        if(V[u] == 0.0){AoA = AoAdelE = 0.0;}
            else{
                AoA = atan2(V[w], V[u]);
            }
        if(Vmag == 0.0){sideslip = asin(0);}
            else{
                double const Vxz = sqrt(V[u]*V[u] + V[w]*V[w]);
                sideslip = atan2(V[v], Vxz);
            }
        AoAdelE = AoA + instance[13];
        AoAdelR = sideslip + instance[14];
    }
    void RK4(double h){
        double y0[15], k[4][15];
        for(int i = 0; i < 15; ++i){y0[i] = instance[i];}
        for(int i = 0; i < 4; ++i){            
            unpackInstance();
            computeAirData();
            cosineMatrix(eulerAngles[p], eulerAngles[q], eulerAngles[r], AoA, sideslip);
            externalForces();
            derivatives();
            for(int j = 0; j < 15; ++j){
                k[i][j] = instanceDeriv[j];
                if(i < 2){instance[j] = y0[j] + 0.5*h*k[i][j];}
                    else if(i == 2){instance[j] = y0[j] + h*k[i][j];}                
            }
        }
        for (int i = 0; i < 15; i++){instance[i] = y0[i] + h*(k[0][i] + 2.0*k[1][i] + 2.0*k[2][i] + k[3][i])/6.0;}
    }
    void initializeDCM(){
        computeAirData();
        // in Dynamics constructor after member initialization
        for (int i=0;i<3;i++){
            for (int j=0;j<3;j++){
                Cm_wb[i][j]  = 0.0;
                Cm_wbT[i][j] = 0.0;
                Cm_bn[i][j]  = 0.0;
                Cm_bnT[i][j] = 0.0;
            }
            Cm_wb[i][i]  = 1.0;
            Cm_wbT[i][i] = 1.0;
            Cm_bn[i][i]  = 1.0;
            Cm_bnT[i][i] = 1.0;
        }
    }
    void cosineMatrix(double phi /*roll*/, double theta /*pitch*/, double psi /*yaw*/, double alpha /*AoA*/, double beta /*sideslip*/){
        //NED -> Body DCM 
        Cm_bn[0][0] = cos(theta)*cos(psi);
        Cm_bn[0][1] = cos(theta)*sin(psi);
        Cm_bn[0][2] = -sin(theta);
        Cm_bn[1][0] = -cos(phi)*sin(psi) + sin(phi)*sin(theta)*cos(psi);
        Cm_bn[1][1] = cos(phi)*cos(psi) + sin(phi)*sin(theta)*sin(psi);
        Cm_bn[1][2] = sin(phi)*cos(theta);
        Cm_bn[2][0] = sin(phi)*sin(psi) + cos(phi)*sin(theta)*cos(psi);
        Cm_bn[2][1] = -sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi);
        Cm_bn[2][2] = cos(phi)*cos(theta);
        //NED -> Body DCM Transpose aka Body -> NED DCM
        Cm_bnT[0][0] = Cm_bn[0][0];
        Cm_bnT[0][1] = Cm_bn[1][0];
        Cm_bnT[0][2] = Cm_bn[2][0];
        Cm_bnT[1][0] = Cm_bn[0][1];
        Cm_bnT[1][1] = Cm_bn[1][1];
        Cm_bnT[1][2] = Cm_bn[2][1];
        Cm_bnT[2][0] = Cm_bn[0][2];
        Cm_bnT[2][1] = Cm_bn[1][2];
        Cm_bnT[2][2] = Cm_bn[2][2];
        //Body -> Wind DCM
        Cm_wb[0][0] = cos(alpha)*cos(beta);
        Cm_wb[0][1] = sin(beta);
        Cm_wb[0][2] = sin(alpha)*cos(beta);
        Cm_wb[1][0] = -cos(alpha)*sin(beta);
        Cm_wb[1][1] = cos(beta);
        Cm_wb[1][2] = -sin(alpha)*sin(beta);
        Cm_wb[2][0] = -sin(alpha);
        Cm_wb[2][1] = 0;
        Cm_wb[2][2] = cos(alpha);
        //Body->Wind DCM Transpose
        Cm_wbT[0][0] = Cm_wb[0][0];
        Cm_wbT[0][1] = Cm_wb[1][0];
        Cm_wbT[0][2] = Cm_wb[2][0];
        Cm_wbT[1][0] = Cm_wb[0][1];
        Cm_wbT[1][1] = Cm_wb[1][1];
        Cm_wbT[1][2] = Cm_wb[2][1];
        Cm_wbT[2][0] = Cm_wb[0][2];
        Cm_wbT[2][1] = Cm_wb[1][2];
        Cm_wbT[2][2] = Cm_wb[2][2];
    }
    double qBar(){
        if(atmosphereON == true){
            double tempurature, pressure;
            if(height < 11000.0){  //calculations from NASA: https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html
                tempurature = 15.04 - 0.00649*(height);
                pressure = 101.29*pow(((tempurature + 273.1) / 288.08), 5.256);
            }else if(height >= 11000.0 && height < 25000.0){
                tempurature = -56.46;
                pressure = 22.65*exp(1.73 - 0.000157*height);
            }else{
                tempurature = -131.21 + .00299*height;
                pressure = 2.488*pow((tempurature + 273.1) / 216.6, -11.388);
            }
            soundBarrier = sqrt(1.4*287.05*(tempurature + 273.1));
            rho = pressure/(0.2869*(tempurature + 273.1));
            return 0.5*rho*pow(Vmag, 2);
        }else{soundBarrier = 0.0; return 0;}
    }
    double MomentCoefficient(Variables axis){
        //Moment Coefficients on body        
        if(mach <= 1.2){
            CldelA = 0.078395 + 0.004599*mach - 0.002012*mach*mach;
            Clp = -1.1339*CldelA;
            CndelR = CmdelE = 12.33 + 70.20*mach - 33.43*mach*mach;
            CnBeta = -CndelR;
            CmAoA = -CndelR;
            Cnr = Cmq = -12.456*CndelR;
        }else{
            CldelA = 0.2061*pow(mach, -1.413);
            Clp = -1.1339*CldelA;
            CndelR = CmdelE = 50.8*pow(mach, -1.40);
            CnBeta = -CndelR;
            CmAoA = -CndelR; 
            Cnr = Cmq = -12.456*CndelR;
        }
        // static stability
        if (CmAoA > 0) CmAoA = -std::abs(CmAoA);
        if (CnBeta < 0) CnBeta =  std::abs(CnBeta);
        // damping
        if (Cmq > 0) Cmq = -std::abs(Cmq);
        if (Cnr > 0) Cnr = -std::abs(Cnr);
        if (Clp > 0) Clp = -std::abs(Clp);
        
        if(Vmag > 0){
            switch(axis){
                case L: return Clp*(W[p]*diameter)/(2.0*Vmag) + CldelA*instance[12]; break;
                case M: return CmAoA*AoA + CmBeta*sideslip + Cmq*(W[q]*diameter)/(2.0*Vmag) + CmdelE*instance[13]; break;
                case N: return CnBeta*sideslip + CnAoA*AoA + Cnr*(W[r]*diameter)/(2.0*Vmag) + CndelR*instance[14]; break;
            }
        }
        return 0;
    }
    void estimateCoeffs(){  
        const double ARwings = b*b/Aref_wing;
        const double ARrudder = 1.072*1.072/Aref_rudder;
        const double betaM = std::max(sqrt(std::max(1 - mach*mach, 1e-9)), 0.3);
        const double cosLamW = cos(55.0*d2r);
        const double cosLamR = cos(45.0*d2r);
        auto safe_sqrt_pos = [](double value) {
            constexpr double eps = 1e-9;
            return sqrt(std::max(value, eps));
        };
        const double Mn_w = mach * cosLamW;
        const double Mn_r = mach * cosLamR;
        // Wing contribution a_w(M)
        double a_w = 0.0;
        if (Mn_w <= 1.0) {
            // subsonic finite-AR (with betaM safely clamped)
            const double denom = 2.0 + sqrt(4.0 + std::pow(ARwings*betaM, 2.0));
            a_w = (2.0*M_PI*ARwings)/denom * cosLamW;
        } else {
            // supersonic-normal (thin-airfoil with finite-AR factor)
            const double denom = safe_sqrt_pos(Mn_w*Mn_w - 1.0);
            a_w = (4.0/denom) * (ARwings/(ARwings + 2.0)) * cosLamW;
        }
        // Rudder contribution a_r(M)
        double a_r = 0.0;
        if (Mn_r <= 1.0) {
            const double denom = 2.0 + sqrt(4.0 + std::pow(ARrudder*betaM, 2.0));
            a_r = (2.0*M_PI*ARrudder)/denom * cosLamR;
        } else {
            const double denom = safe_sqrt_pos(Mn_r*Mn_r - 1.0);
            a_r = (4.0/denom) * (ARrudder/(ARrudder + 2.0)) * cosLamR;
        }
        C_Nalpha = 2.0 + 0.90*(Aref_wing/Aref)*a_w + 0.85*(Aref_rudder/Aref)*a_r;
        K_D = 1.1*(1 + 0.15*std::max(mach - 1, 0.0));
        const double C_Nlin = C_Nalpha*AoA;
        const double C_Nhigh = 0.5*C_Nalpha*sin(2*AoA);
        const double alpha_b = 15.0 * d2r;
        const double dAlpha  =  5.0 * d2r;
        const double weight  = 0.5 * (1.0 + tanh((std::abs(AoA) - alpha_b) / dAlpha));
        const double C_Nmag = (1 - weight)*C_Nlin + weight*C_Nhigh;
        //Drag Coefficient
        C_Drag = C_DragBase + C_DragBase*exp(-pow((mach - 1.0)/0.25, 2))
                + 0.5*C_DragBase*(1.0 - exp(-(mach - 1)/0.9)) + 0.02*(mach - 1.0)*H(mach - 1.0)
                + K_D*sin(AoA)*sin(AoA);
        //Lift/Sideforce Coefficient
        C_LiftPitch = C_Nmag*cos(sideslip);
        C_LiftYaw = -C_Nmag*sin(sideslip);
    }
    void unpackInstance(){
        for(int i = 0; i < 3; i++){
            V[i] = instance[i];
            W[i] = instance[i+3];
            eulerAngles[i] = instance[i+6];
            position[i] = instance[i+9];
        }
        height = -position[z];
    }
    double H(double machValue){
        if(machValue < 1){return 0;}
        return 1;
    }
    static double wrapToPi(double a){
        while(a >  M_PI){a -= 2.0*M_PI;}
        while(a < -M_PI){a += 2.0*M_PI;}
        return a;
    }
    static double clamp(double x, double low, double hi) {return (x < low) ? low : (x > hi) ? hi : x;}
};