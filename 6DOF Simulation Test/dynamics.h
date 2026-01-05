#include <iostream>
#include <cmath>
#include <string>
using namespace std;

class dynamics{
    public:
    dynamics(double m, double length, double radius, int object){
        mass = m;
        AoA = sideslip = 0;
        for(int i = 0; i < 3; i++){
            instance[i] = V[i] = 0;
            instance[i+3] = W[i] = 0;
            instance[i+6] = eulerAngles[i] = 0;
            instance[i+9] = position[i] = 0;    
            moments[i] = 0;
        }
        switch(object){
            case 1: I[x] = I[y] = I[z] = 0.4*mass*pow(radius, 2); break;
            case 2: 
                I[x] = mass*pow(radius, 2)/2;
                I[y] = I[z] = mass*pow(radius, 2)/4 + mass*pow(length, 2)/12;
                break;
            default: break;
        }
    }
    void initialConditions(double i = 0, double j = 0, double k = 0, double roll = 0, double pitch = 0, double yaw = 0, double phi = 0, double theta = 0, double psi = 0, double Px = 0, double Py = 0, double Pz = 0){
        instance[0] = V[u] = i;
        instance[1] = V[v] = j; 
        instance[2] = V[w] = k;
        instance[3] = W[p] = roll; 
        instance[4] = W[q] = pitch; 
        instance[5] = W[r] = yaw;
        instance[6] = eulerAngles[p] = phi; 
        instance[7] = eulerAngles[q] = theta; 
        instance[8] = eulerAngles[r] = psi;
        instance[9] = position[x] = Px; 
        instance[10] = position[y] = Py; 
        instance[11] = position[z] = Pz;
    }
    double atmosphericDensity(){
        double tempurature, pressure;
        if(height < 11000){  //calculations from NASA: https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html
            tempurature = 15.04 - 0.00649*(height);
            pressure = 101.29*pow(((tempurature + 273.1) / 288.08), 5.256);
        }else if(height >= 11000 && height < 25000){
            tempurature = -56.46;
            pressure = 22.65*exp(1.73 - 0.000157*height);
        }else{
            tempurature = -131.21 + .00299*height;
            pressure = 2.488*pow((tempurature + 273.1) / 216.6, -11.388);
        }
        return pressure/(0.2869*(tempurature + 273.1));
    }
    void eulerKinematics(){
        instanceDeriv[6] = (W[p] + sin(eulerAngles[p])*tan(eulerAngles[q])*W[q] + cos(eulerAngles[p])*tan(eulerAngles[q])*W[r]);
        instanceDeriv[7] = (cos(eulerAngles[p])*W[q] - sin(eulerAngles[p])*W[r]);
        instanceDeriv[8] = (sin(eulerAngles[p])/cos(eulerAngles[q])*W[q] + cos(eulerAngles[p])/cos(eulerAngles[q])*W[r]); 
        if(V[u] == 0){AoA = atan(0);}
            else{AoA = atan2(-V[w], V[u]);}
        if(sqrt(V[u]*V[u] + V[v]*V[v] + V[w]*V[w]) == 0){sideslip = asin(0);}
            else{sideslip = asin(V[v]/sqrt(V[u]*V[u] + V[v]*V[v] + V[w]*V[w]));}
    }
    void accelerations(){
        instanceDeriv[0] = Fb[x]/mass - g*Cm_bn[0][2] - V[w]*W[q] + V[v]*W[r];
        instanceDeriv[1] = Fb[y]/mass - g*Cm_bn[1][2] - V[u]*W[r] + V[w]*W[p];
        instanceDeriv[2] = Fb[z]/mass - g*Cm_bn[2][2] - V[v]*W[p] + V[u]*W[q];
        instanceDeriv[3] = (W[q]*W[r]*((I[y] - I[z])/I[x]) + moments[L]/I[x]);
        instanceDeriv[4] = (W[p]*W[r]*((I[z] - I[x])/I[y]) + moments[M]/I[y]);
        instanceDeriv[5] = (W[p]*W[q]*((I[x] - I[y])/I[z]) + moments[N]/I[z]);
    }
    void externalForces(double Cd, double Aref){
        dragForce = 0.5*Cd*atmosphericDensity()*pow(sqrt(V[x]*V[x] + V[y]*V[y] + V[z]*V[z]), 2)*Aref;
        sideForce = liftForce = 0;
        Fb[x] = -(Cm_wb[0][0]*dragForce + Cm_wb[1][0]*sideForce + Cm_wb[2][0]*liftForce);
        Fb[y] = -(Cm_wb[0][1]*dragForce + Cm_wb[1][1]*sideForce + Cm_wb[2][1]*liftForce);
        Fb[z] = -(Cm_wb[0][2]*dragForce + Cm_wb[1][2]*sideForce + Cm_wb[2][2]*liftForce);
    }
    void motion(){
            instanceDeriv[9] = Cm_bn[0][0]*V[x] + Cm_bn[0][1]*V[y] + Cm_bn[0][2]*V[z];
            instanceDeriv[10] = Cm_bn[1][0]*V[x] + Cm_bn[1][1]*V[y] + Cm_bn[1][2]*V[z];
            instanceDeriv[11] = Cm_bn[2][0]*V[x] + Cm_bn[2][1]*V[y] + Cm_bn[2][2]*V[z];
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
    void iterate(double dt, double Cd, double Aref){
        cosineMatrix(eulerAngles[p], eulerAngles[q], eulerAngles[r], AoA, sideslip);
        externalForces(Cd, Aref);
        accelerations();
        eulerKinematics();
        motion();
        simpleEuler(dt);
        unpackInstance();
    }
    double getInstance(int i){
        return instance[i];
    }
    private:
    double Cm_bn[3][3], Cm_wb[3][3],   //Cosine matrices for body-normal, normal-body, wind-body
            I[3],       //Inertia
            position[3],//XYZ vector
            Fb[3],      //force vector XYZ
            V[3], Vdot[3],//velocity vectors u, v, w
            W[3], Wdot[3],//rotation vectors p, q, r
            eulerAngles[3], eulerRates[3],  //roll: phi, pitch: theta, yaw: psi 
            moments[3],    //Moment vectors L, M, N
            instance[12], instanceDeriv[12],  //Arrays to store simulation state and rates for integration
            g{9.81}, mass, AoA, sideslip, height, soundBarrier, dragForce, sideForce, liftForce;
    int x{0}, y{1}, z{2}, u{0}, v{1}, w{2}, p{0}, q{1}, r{2}, L{0}, M{1}, N{2}; //variables for array elements for ease of readability
    //Various calculators
    void simpleEuler(double h){
        for(int i = 0; i < 12; i++){
            instance[i] += instanceDeriv[i]*h;
        }
    }
    void cosineMatrix(double phi /*roll*/, double theta /*pitch*/, double psi /*yaw*/, double alpha /*AoA*/, double beta /*sideslip*/){
        //Body-Normal matrix calculations
        Cm_bn[0][0] = cos(theta)*cos(psi);
        Cm_bn[0][1] = cos(theta)*sin(psi);
        Cm_bn[0][2] = -sin(theta);
        Cm_bn[1][0] = -cos(phi)*sin(psi) + sin(phi)*sin(theta)*cos(psi);
        Cm_bn[1][1] = cos(phi)*cos(psi) + sin(phi)*sin(theta)*sin(psi);
        Cm_bn[1][2] = sin(phi)*cos(theta);
        Cm_bn[2][0] = sin(phi)*sin(psi) + cos(phi)*sin(theta)*cos(psi);
        Cm_bn[2][1] = -sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi);
        Cm_bn[2][2] = cos(phi)*cos(theta);
        //Wind-Body matrix calculations
        Cm_wb[0][0] = cos(alpha)*cos(beta);
        Cm_wb[0][1] = sin(beta);
        Cm_wb[0][2] = sin(alpha)*cos(beta);
        Cm_wb[1][0] = -cos(alpha)*sin(beta);
        Cm_wb[1][1] = cos(beta);
        Cm_wb[1][2] = sin(alpha)*sin(beta);
        Cm_wb[2][0] = -sin(alpha);
        Cm_wb[2][1] = 0;
        Cm_wb[2][2] = cos(alpha);
    }
};