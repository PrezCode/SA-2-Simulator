#include <iostream>
#include <cmath>
#include <string>

class Dynamics{
    public:
    enum Shape{Sphere, Cylinder, RectPrism};
    Dynamics(double m, double length, double width, double height_construct, double radius, Shape object){
        mass = m;
        AoA = sideslip = 0;
        for(int i = 0; i < 3; i++){
            instance[i] = V[i] = 0;
            instance[i+3] = W[i] = 0;
            instance[i+6] = eulerAngles[i] = 0;
            instance[i+9] = position[i] = 0;    
            moments[i] = Fb[i] = 0;
        }
        switch(object){
            case Sphere: I[x] = I[y] = I[z] = 0.4*mass*pow(radius, 2); I[xz] = 0; break;
            case Cylinder: 
                I[x] = mass*pow(radius, 2)/2;
                I[y] = I[z] = mass*pow(radius, 2)/4 + mass*pow(length, 2)/12;
                I[xz] = 0;
                break;
            case RectPrism:
                I[x] = mass*(pow(width, 2) + pow(height_construct, 2))/12;
                I[y] = mass*(pow(length, 2) + pow(height_construct, 2))/12;
                I[z] = mass*(pow(length, 2) + pow(width, 2))/12;
                I[xz] = 0;
                break;
        }
    }
    void initialConditions(double i, double j, double k, double roll, double pitch, double yaw, double phi, double theta, double psi, double Px, double Py, double Pz){
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
        Vmag = sqrt(i*i + j*j + k*k);
    }
    void externalForces(double Aref, double b, double c, double* coefficients){
        dragForce = coefficients[Cd]*qBar()*Aref;
        sideForce = liftForce = 0;
        Fb[x] = -(Cm_wb[0][0]*dragForce + Cm_wb[1][0]*sideForce + Cm_wb[2][0]*liftForce);
        Fb[y] = -(Cm_wb[0][1]*dragForce + Cm_wb[1][1]*sideForce + Cm_wb[2][1]*liftForce);
        Fb[z] = -(Cm_wb[0][2]*dragForce + Cm_wb[1][2]*sideForce + Cm_wb[2][2]*liftForce);
        moments[L] = MomentCoefficient(L, b, c, coefficients)*qBar()*Aref*b;
        moments[M] = MomentCoefficient(M, b, c, coefficients)*qBar()*Aref*c;
        moments[N] = MomentCoefficient(N, b, c, coefficients)*qBar()*Aref*b;
    }
    double MomentCoefficient(int i, double span, double chord, double* coefficients){
        switch(i){
            case 0: return coefficients[Clp]*W[p]*span/(2*Vmag) + coefficients[Clr]*W[r]*span;
            case 1: return coefficients[Cmq]*W[q]*chord/(2*Vmag);
            case 2: return coefficients[Cnp]*W[p]*span/(2*Vmag) + coefficients[Cnr]*W[r]*span/(2*Vmag);
        }
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
        /*
        instanceDeriv[3] = (W[q]*W[r]*((I[y] - I[z])/I[x]) + moments[L]/I[x]);
        instanceDeriv[4] = (W[p]*W[r]*((I[z] - I[x])/I[y]) + moments[M]/I[y]);
        instanceDeriv[5] = (W[p]*W[q]*((I[x] - I[y])/I[z]) + moments[N]/I[z]);
        */
        //Euler Angle rates
        instanceDeriv[6] = (W[p] + sin(eulerAngles[p])*tan(eulerAngles[q])*W[q] + cos(eulerAngles[p])*tan(eulerAngles[q])*W[r]);
        instanceDeriv[7] = (cos(eulerAngles[p])*W[q] - sin(eulerAngles[p])*W[r]);
        instanceDeriv[8] = (sin(eulerAngles[p])/cos(eulerAngles[q])*W[q] + cos(eulerAngles[p])/cos(eulerAngles[q])*W[r]); 
        //XYZ Motion
        instanceDeriv[9] = Cm_bn[0][0]*V[u] + Cm_bn[1][0]*V[v] + Cm_bn[2][0]*V[w];
        instanceDeriv[10] = Cm_bn[0][1]*V[u] + Cm_bn[1][1]*V[v] + Cm_bn[2][1]*V[w];
        instanceDeriv[11] = Cm_bn[0][2]*V[u] + Cm_bn[1][2]*V[v] + Cm_bn[2][2]*V[w];
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
    void iterate(double dt, double Aref, double wingSpan, double wingChord, double* coefficients){
        cosineMatrix(eulerAngles[p], eulerAngles[q], eulerAngles[r], AoA, sideslip);
        externalForces(Aref, wingSpan, wingChord, coefficients);
        derivatives();
        //simpleEuler(dt);
        RK4(dt);
        computeAirData();
        unpackInstance();
    }
    double getInstance(int i){
        return instance[i];
    }
    double getAirData(int i){
        switch(i){
            case 0: return Vmag;
            case 1: return mach;
            case 2: return AoA;
            case 3: return sideslip;
            default: return 0;
        }
    }
    private:
    double Cm_bn[3][3], Cm_wb[3][3],   //Cosine matrices for body-normal, normal-body, wind-body
            I[4],       //Inertia, XX, YY, ZZ, XZ
            position[3],//XYZ vector
            Fb[3],      //force vector XYZ
            V[3],       //velocity vectors u, v, w
            W[3],       //rotation vectors p, q, r
            eulerAngles[3], eulerRates[3],  //roll: phi, pitch: theta, yaw: psi 
            moments[3],    //Moment vectors L, M, N
            instance[12], instanceDeriv[12],  //Arrays to store simulation state and rates for integration
            g{9.81}, mass, AoA, sideslip, height, soundBarrier, mach, dragForce, sideForce, liftForce, Vmag;
    //variables for array elements for ease of readability
    int x{0}, y{1}, z{2}, xz{3},
        u{0}, v{1}, w{2},
        p{0}, q{1}, r{2},
        L{0}, M{1}, N{2},
        Cd{0}, Clp{1}, Clr{2}, Cmq{3}, Cnp{4}, Cnr{5}; 
    //Various calculators
    void computeAirData(){
        Vmag = sqrt(V[u]*V[u] + V[v]*V[v] + V[w]*V[w]);
        mach = Vmag/soundBarrier;
        if(V[u] == 0){AoA = atan(0);}
            else{AoA = atan2(V[w], V[u]);}
        if(Vmag == 0){sideslip = asin(0);}
            else{sideslip = asin(V[v]/Vmag);}
    }
    void simpleEuler(double h){
        for(int i = 0; i < 12; i++){
            instance[i] += instanceDeriv[i]*h;
        }
    }
    void RK4(double h){
        for(int i = 0; i < 12; i++){
            double k1 = instanceDeriv[i];
            double k2 = instanceDeriv[i] + h*k1/2;
            double k3 = instanceDeriv[i] + h*k2/2;
            double k4 = instanceDeriv[i] + h*k3;
            instance[i] += h*(k1 + 2*k2 + 2*k3 + k4)/6;
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
    double qBar(){
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
        soundBarrier = sqrt(1.4*287.05*(tempurature + 273.1));
        return 0.5*(pressure/(0.2869*(tempurature + 273.1)))*pow(Vmag, 2);
    }
};