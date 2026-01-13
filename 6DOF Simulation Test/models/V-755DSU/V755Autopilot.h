#include <algorithm>
#include <cmath>

class V755Autopilot {
public:
    V755Autopilot() = default;

    // Returns commanded actuator angles (radians):
    // commands[0] unused (kept for your 3-channel convention)
    // commands[1] Channel I (e.g., pitch)
    // commands[2] Channel II (e.g., yaw)
    double autopilotCommands(int i) const { return commands[i]; }

    // Main update
    //
    // Inputs:
    //  rho            kg/m^3
    //  V              true airspeed (m/s)  <-- recommended (avoids Mach/a ambiguity)
    //  pitchRate_rad_s, yawRate_rad_s   body rates (rad/s)
    //  gY_m_s2, gZ_m_s2                  measured linear accel components (m/s^2)
    //  Vpk1, Vpk2                        radio-command "applied voltage" into R1-8 path (V)
    //  dt                               timestep (s)
    //
    // Notes:
    //  - If your sim provides accel already in g's, pass g*9.80665 or change conversion below.
    //  - If Vpk1/Vpk2 are NOT actual volts across ~2870 ohm, you must adapt the V->I mapping.
    void autopilotGenerate(double rho, double V,
                           double pitchRate_rad_s, double yawRate_rad_s,
                           double gY_m_s2, double gZ_m_s2,
                           double Vpk1, double Vpk2,
                           double dt)
    {
        // Guard dt
        if (dt <= 0.0) return;

        constexpr double deg2rad = PI / 180.0;
        constexpr double rad2deg = 180.0 / PI;

        // ===== Constants from your docs =====
        constexpr double Kdg_deg_per_degps = 0.7;   // deg/(deg/s)
        constexpr double Kxu_deg_per_g     = 13.0;  // deg/g

        // US-17 "slope" spec: 0.25 mA / (mA*turn)
        constexpr double us17_slope_mA_per_mAturn = 0.25;
        constexpr double N_cmd_turns = 800.0;

        // Radio-command injection resistor (your earlier assumption)
        constexpr double R1_8_ohm = 2870.0;

        // AC-5 schedule reference gain (from your item-23 derived mapping usage)
        constexpr double Kref_deg_per_mA = 79.4;

        // AC-5 working range (kgf/cm^2) from your text
        constexpr double Pq_min_kgcm2 = 1.6;
        constexpr double Pq_max_kgcm2 = 5.0;

        // Surface/actuator limit for channels I/II
        constexpr double surf_lim_deg = 28.0;

        // Simple actuator dynamics (tunable)
        constexpr double tau_act_s = 0.10;        // first-order lag time constant
        constexpr double rate_lim_deg_s = 200.0;  // actuator deflection rate limit

        // ===== 1) Dynamic pressure =====
        // q [Pa] = 0.5 * rho * V^2
        const double q_Pa = 0.5 * rho * V * V;

        // Convert to kgf/cm^2
        double Pq_kgcm2 = q_Pa / 98066.5;

        // Clamp to stated AC-5 operating range
        Pq_kgcm2 = std::clamp(Pq_kgcm2, Pq_min_kgcm2, Pq_max_kgcm2);

        // AC-5 hyperbolic schedule (per your description)
        const double Kp_deg_per_mA = Kref_deg_per_mA * (Pq_min_kgcm2 / Pq_kgcm2);

        // ===== 2) Unit conversions =====
        const double pitchRate_deg_s = pitchRate_rad_s * rad2deg;
        const double yawRate_deg_s   = yawRate_rad_s   * rad2deg;

        const double gY = gY_m_s2 / 9.80665;  // convert m/s^2 to g
        const double gZ = gZ_m_s2 / 9.80665;

        // ===== 3) Radio command: V -> I_in(mA) -> US-17 dI_out(mA) =====
        auto V_to_Iin_mA = [&](double Vin) -> double {
            return (Vin / R1_8_ohm) * 1000.0; // mA
        };

        const double Iin1_mA = V_to_Iin_mA(Vpk1);
        const double Iin2_mA = V_to_Iin_mA(Vpk2);

        // Ampere-turns in (mA*turn) = N * Iin(mA)
        const double NI1_mAturn = N_cmd_turns * Iin1_mA;
        const double NI2_mAturn = N_cmd_turns * Iin2_mA;

        // US-17 incremental output current (mA)
        // dIout = slope * (N * Iin)
        const double dIout1_mA = us17_slope_mA_per_mAturn * NI1_mAturn;
        const double dIout2_mA = us17_slope_mA_per_mAturn * NI2_mAturn;

        // ===== 4) Build raw commanded deflection in degrees =====
        // Structure: + radio command term + accel term - rate damping term
        double cmd1_deg = (Kp_deg_per_mA * dIout1_mA) + (Kxu_deg_per_g * gZ) - (Kdg_deg_per_degps * pitchRate_deg_s);
        double cmd2_deg = (Kp_deg_per_mA * dIout2_mA) + (Kxu_deg_per_g * gY) - (Kdg_deg_per_degps * yawRate_deg_s);

        // Hard command limit (mechanical stop)
        cmd1_deg = std::clamp(cmd1_deg, -surf_lim_deg, surf_lim_deg);
        cmd2_deg = std::clamp(cmd2_deg, -surf_lim_deg, surf_lim_deg);

        // ===== 5) Actuator dynamics: first-order lag + rate limit =====
        // Track actuator states in degrees internally
        act1_deg = stepActuator(act1_deg, cmd1_deg, dt, tau_act_s, rate_lim_deg_s);
        act2_deg = stepActuator(act2_deg, cmd2_deg, dt, tau_act_s, rate_lim_deg_s);

        // ===== 6) Output in radians =====
        commands[0] = 0.0;
        commands[1] = act1_deg * deg2rad;
        commands[2] = act2_deg * deg2rad;
    }

private:
    // Actuator state (deg)
    double act1_deg {0.0};
    double act2_deg {0.0};

    // Output commands (rad)
    double commands[3] {0.0, 0.0, 0.0};

    // One step of a rate-limited first-order actuator
    static double stepActuator(double y_deg, double u_deg, double dt,
                              double tau, double rateLim_deg_s)
    {
        // First-order target rate (deg/s)
        const double dy_des_deg_s = (u_deg - y_deg) / std::max(tau, 1e-6);

        // Apply rate limit
        const double dy_lim_deg_s = std::clamp(dy_des_deg_s, -rateLim_deg_s, rateLim_deg_s);

        // Integrate
        return y_deg + dy_lim_deg_s * dt;
    }
};

/*
#include <iostream>
#include <cmath>
#include <string>

class V755Autopilot{
    public:
    V755Autopilot(){}
    void autopilotGenerate(double rho, double mach, double pitchRate, double yawRate, double gY, double gZ, double K1, double K2){
        double Kdg = 0.7;
        double Kxu = 13;
        K1input = 0.25*800*((K1/2870)*1000);
        K2input = 0.25*800*((K2/2870)*1000);
        Pq = rho*pow(mach, 2)/2;
        double Kref = 79.4;
        double KPq = Kref*(1.6/Pq);
        
        commands[0] = 0;
        commands[1] = (KPq*K1input + Kxu*gZ/g - Kdg*pitchRate/deg_rad)*deg_rad;
        if(commands[1] > 28*deg_rad){commands[1] = 28*deg_rad;}
            else if(commands[1] < -28*deg_rad){commands[1] = -28*deg_rad;}
        commands[2] = (KPq*K2input + Kxu*gY/g - Kdg*yawRate)*deg_rad;
        if(commands[2] > 28*deg_rad){commands[2] = 28*deg_rad;}
            else if(commands[2] < -28*deg_rad){commands[2] = -28*deg_rad;}
    }
    double autopilotCommands(int i){
        return commands[i];
    }
    
    double wrapPi(double a){
        while (a >  M_PI) a -= 2.0*M_PI;
        while (a < -M_PI) a += 2.0*M_PI;
        return a;
    }
    private:
    double commands[3] = {0,0,0}, Pq, deg_rad{M_PI/180}, K1input, K2input, g{9.81};
};
*/
