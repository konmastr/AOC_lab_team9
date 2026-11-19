// Lab 3: LQR balance of the pendulum, TORQUE command.   (student sketch: complete the TODOs)
//
// Given: the settings below, the energy swing-up helper (SwingUp, the Lab 2 law), the prints, the serial keys
//        ('1', '2', 'u', 'h') and the CSV logging.
// Yours (the TODOs, in this order):
//   Task 3.1   the Lab 3 model / linearization / discretization TODOs of PendulumModel.cpp (your library)
//   Task 3.2   myLqrGain(): the discrete LQR gain from the Riccati recursion iterated to convergence (DARE)
//   Task 3.2a  the objects and the global variables
//   Task 3.2b  updateLqr(): the discrete model, the call of myLqrGain(), the closed-loop eigenvalue check
//   Task 3.2c  setup(): start-up of the motor, the pendulum, the swing-up helper and the LQR
//   Task 3.2d  loop(): the control loop with the LQR law and its torque limit
// Enter your own rig-card calibration, mechanical parameters and stored alignment values below.
//
//   theta = 0 : hanging down      theta = pi : upright      (positive theta = positive torque)
//
// The control task is selected at runtime (keys '1' and '2' in the Serial Monitor):
//   EQUILIBRIUM 1  HANGING (stable)   : the pendulum is already down, so there is no swing-up; the LQR holds
//                      it at theta = 0 and rejects disturbances.
//   EQUILIBRIUM 2  UPRIGHT (unstable) : the LQR balances the pendulum at theta = pi. It gets there with
//                      USE_SWING_UP 1 (key 'u'): the energy swing-up (given) pumps it up from the bottom;
//                      USE_SWING_UP 0: the motor stays at zero torque until YOU bring the pendulum by hand
//                      close to the top (pendulum.updateMode() == BALANCE), as in Lab 2, Part 1.
//
// LQR: on the model linearized at the chosen equilibrium, state x = [e, omega], e = theta - theta_eq:
//                       x[k+1] = Ad*x[k] + Bd*tau[k]          (zero-order hold, Ts = control period)
//                       cost   = sum( Q_ANGLE*e^2 + Q_SPEED*omega^2 + R*tau^2 )
//                       tau    = -(K[0]*e + K[1]*omega)
//                 K comes from the discrete algebraic Riccati equation (DARE):
//                       S = R + B'PB
//                       P <- Q + A'PA - (A'PB)(A'PB)'/S          (iterate until P stops changing)
//                       K = B'PA / S
//                 YOU compute K on the Arduino (Task 3.2). At start-up the sketch checks the closed loop:
//                 the eigenvalues of Ad - Bd*K must be inside the unit circle.
//
// Procedure:
//   1. Fill in the calibration results and the pendulum parameters below (Lab 1). Paste the STORED_ lines.
//   2. Let the pendulum hang at rest, open the Serial Monitor (115200) and send 's'. Send 'x' to stop (0 V).
//      Keys: '1'/'2' choose the equilibrium, 'u' the swing-up helper, 'p' pauses the log, 'h' help.
//
// Log (CSV): t_s,theta_rad,omega_rad_s,torque_Nm,voltage_V,mode   (mode 0 = swing-up, 1 = balance)
//
// WARNING: the pendulum moves fast and can hit things. Keep the area clear.

#include <FOC_lar.h>
#include <Pendulum_lar.h>

#define EQUILIBRIUM 1  // initial equilibrium: 1 HANGING (stable), 2 UPRIGHT (switch at runtime)
#define USE_SWING_UP 1 // upright: 1 energy swing-up + LQR, 0 LQR only (lift by hand) (key 'u')

// Motor definitions

#define MOTOR_POLE_PAIRS 7
#define MOTOR_SUPPLY_VOLTAGE 12.0f
#define MOTOR_VOLTAGE_LIMIT 3.0f // limits the torque: max torque = TORQUE_CONSTANT * this
#define SENSOR_ALIGN_VOLTAGE 1.0f

#define STORED_ZERO_ELECTRIC_ANGLE NOT_SET // paste the values printed by motor_AS5600
#define STORED_SENSOR_DIRECTION Direction::UNKNOWN

// Calibration results (Lab 1)

#define TORQUE_CONSTANT NOT_SET // N*m/V    <- REPLACE
#define MOTOR_KE NOT_SET        // V*s/rad  <- REPLACE (low-speed Ke)
#define USE_BACK_EMF_COMP 1     // 1: the torque command includes the back-EMF compensation

// Pendulum parameters (Lab 1)

#define ROD_MASS NOT_SET         // kg   <- REPLACE
#define POINT_MASS NOT_SET       // kg   <- REPLACE
#define PENDULUM_LENGTH NOT_SET  // m    <- REPLACE
#define VISCOUS_FRICTION NOT_SET // N*m*s/rad  <- REPLACE (ring-down)
#define COULOMB_FRICTION NOT_SET // N*m, tau_c  <- REPLACE (ring-down)
#define USE_COULOMB_COMP 1       // 1: add tau_c*sign(omega) to the torque (the LQR model has no Coulomb term)

// Control

#define CONTROL_PERIOD_S 0.01f
#define VELOCITY_ALPHA 0.5f
#define RUN_TIME_S 3600.0f

#define CAPTURE_ANGLE 0.4f
#define RELEASE_ANGLE 1.0f

// Energy-based swing-up (given: the same law you implement in Lab 2)
#define SWING_ENERGY_BAND 0.10f  // energy error, as a fraction of E_up, over which the torque ramps to full
#define SWING_ENERGY_RATIO 1.03f // target energy / energy at the top
#define SWING_KICK_SPEED 0.3f    // below this |omega| (rad/s) the last direction is kept

// LQR cost:  R = LQR_R_SCALE / max_torque^2   (a torque of max_torque costs as much as LQR_R_SCALE)
#define LQR_Q_ANGLE 10.0f // weight of phi^2          (1/rad^2 order)
#define LQR_Q_SPEED 0.1f  // weight of omega^2
#define LQR_R_SCALE 1.0f  // larger = gentler torque

// Discrete LQR gain by iterating the Riccati equation. Returns the number of iterations (0 = no convergence).
int myLqrGain(const float Ad[2][2], const float Bd[2], float q_angle, float q_speed, float r, float K[2])
{
    float P[2][2] = {{q_angle, 0.0f}, {0.0f, q_speed}}; // start with P = Q

    for (int it = 1; it <= 20000; it++)
    {
        // TODO (Task 3.2): one iteration of the Riccati recursion (handout, section 1.3).
        //   T = P*Bd (2-vector),  S = r + Bd'*T (scalar),  L = Ad'*T (2-vector, = A'PB)
        //   Pn = Q + Ad'*P*Ad - L*L'/S       (Q = diag(q_angle, q_speed))
        //   P <- Pn
        //   if the largest change of P is below 1e-6 * max(P[0][0], 1):
        //       recompute T and S with the final P, set K = (T'*Ad)/S  (K[0], K[1]) and return it
        // Remove the next line when you are done.
        (void)Ad;
        (void)Bd;
        (void)r;
        (void)K;
    }
    return 0;
}

// GIVEN: energy-based swing-up (same law as the Lab 2 implementation)

struct SwingUp
{
    const PendulumModel *model = 0;
    float max_torque = 0.0f;
    float direction = 1.0f; // push direction, kept when the pendulum is almost at rest

    void begin(const PendulumModel &m, float max_torque_nm)
    {
        model = &m;
        max_torque = max_torque_nm;
        direction = 1.0f;
    }

    void reset()
    {
        direction = 1.0f;
    }

    float torque(float theta, float omega)
    {
        if (fabsf(omega) > SWING_KICK_SPEED)
        {
            direction = pendulum_sign(omega);
        }
        float e_up = model->energyUpright();
        float e_error = (SWING_ENERGY_RATIO * e_up - model->energy(theta, omega)) / (SWING_ENERGY_BAND * e_up);
        return pendulum_clamp(e_error, -1.0f, 1.0f) * max_torque * direction;
    }
};

// TODO (Task 3.2a): declare the objects and the global variables of the sketch. The rest of the sketch uses these
// names, so it compiles only once they exist (look at how Lab 2 declares foc and the pendulum object):
//   1. foc: the motor, an FOC_lar object built with MOTOR_POLE_PAIRS and the pins of the SimpleFOCMini v1.0 on an
//      Arduino UNO, in the order IN1, IN2, IN3, EN = 11, 10, 9, 8.
//   2. pendulum: the pendulum, a Pendulum_lar object (model, measured state, timing, log, swing-up/balance switch).
//   3. swing: the energy swing-up helper, a SwingUp object (the struct above).
//   4. max_torque: a float, the largest torque the motor can apply (N*m). It is set in setup(); start it at 0.
//   5. K: an array of 2 floats, the LQR gain used by tau = -(K[0]*e + K[1]*omega); start it at {0, 0}.
//   6. equilibrium: an int, the equilibrium in use, 1 hanging or 2 upright; start it at EQUILIBRIUM
//      (the keys '1' and '2' change it at runtime).
//   7. lqr_iterations: an int, the number of iterations of your Riccati solver (0 = it did not converge); start
//      it at 0.
bool use_swing_up = USE_SWING_UP; // upright: true = energy swing-up + LQR, false = LQR only (lift by hand)

// Recomputes K for the selected equilibrium and prints the closed-loop check
void updateLqr(void)
{
    float r_torque = LQR_R_SCALE / (max_torque * max_torque);

    // TODO (Task 3.2b): the discrete model x[k+1] = Ad*x[k] + Bd*tau[k] of the selected equilibrium, with the torque
    // held over each control period (zero-order hold, Ts = CONTROL_PERIOD_S). Fill Ad and Bd with the library
    // (your Lab 3 TODOs of PendulumModel.cpp):
    //   - upright (equilibrium == 2): pendulum.discretizeUpright(Ts, Ad, Bd); the state is [phi, omega], phi = theta - pi
    //   - hanging (equilibrium == 1): pendulum.discretizeAt(0.0f, Ts, Ad, Bd); the state is [theta, omega]
    // Check in the print below: Ad[1][0] is positive at the top and negative at the bottom.
    float Ad[2][2] = {{0.0f, 0.0f}, {0.0f, 0.0f}}; // <- fill (the prints below use Ad and Bd)
    float Bd[2] = {0.0f, 0.0f};

    // TODO (Task 3.2b): the LQR gain. Call your myLqrGain(Ad, Bd, q_angle, q_speed, r, K) with the weights
    // LQR_Q_ANGLE, LQR_Q_SPEED and r_torque (above), so that it writes the gain into the global K, and store the
    // number of iterations it returns in lqr_iterations (0 = it did not converge: the sketch stops in setup()).
    (void)r_torque;     // <- remove once you use r_torque
    lqr_iterations = 0; // <- replace
    if (lqr_iterations == 0)
    {
        Serial.println(F("#ERROR: your Riccati iteration did not converge (or is not implemented)"));
        return;
    }

    Serial.print(F("#Ad = ["));
    Serial.print(Ad[0][0], 5);
    Serial.print(F(", "));
    Serial.print(Ad[0][1], 5);
    Serial.print(F("; "));
    Serial.print(Ad[1][0], 5);
    Serial.print(F(", "));
    Serial.print(Ad[1][1], 5);
    Serial.print(F("]  Bd = ["));
    Serial.print(Bd[0], 6);
    Serial.print(F(", "));
    Serial.print(Bd[1], 4);
    Serial.println(F("]"));
    Serial.print(F("#your K = ["));
    Serial.print(K[0], 5);
    Serial.print(F(", "));
    Serial.print(K[1], 5);
    Serial.print(F("]  ("));
    Serial.print(lqr_iterations);
    Serial.println(F(" iterations)"));

    // TODO (Task 3.2b): check the closed loop x[k+1] = (Ad - Bd*K) x[k]. It is stable when both eigenvalues of the
    // 2x2 matrix M = Ad - Bd*K are inside the unit circle (|eigenvalue| < 1).
    //   1. M[i][j] = Ad[i][j] - Bd[i] * K[j]   (four numbers M00, M01, M10, M11)
    //   2. The eigenvalues solve  lambda^2 - tr*lambda + det = 0,  with tr = M00 + M11, det = M00*M11 - M01*M10.
    //      disc = tr^2 - 4*det:
    //      - disc >= 0: two real eigenvalues (tr + sqrt(disc))/2 and (tr - sqrt(disc))/2: take the larger
    //        absolute value (fabsf(), fmaxf(), sqrtf());
    //      - disc < 0: a complex pair with the same magnitude, |lambda|^2 = det: take sqrtf(fabsf(det)).
    //   3. Store the largest magnitude in lambda: the print below reports "stable" when lambda < 1.
    float lambda = NAN; // <- replace
    Serial.print(F("#closed loop: max |eigenvalue| = "));
    Serial.print(lambda, 4);
    Serial.println(lambda < 1.0f ? F("  (stable)") : F("  (UNSTABLE!)"));
}

void printState(void)
{
    Serial.print(equilibrium == 1 ? F("#equilibrium 1: hanging") : F("#equilibrium 2: upright"));
    if (equilibrium == 2)
    {
        Serial.print(use_swing_up ? F(", swing-up + LQR") : F(", LQR only: lift the pendulum by hand"));
    }
    Serial.println();
}

void printHelp(void)
{
    Serial.println(F("#Serial commands:"));
    pendulum.printHelp();
    Serial.println(F("1 : equilibrium 1: hanging, LQR holds theta = 0"));
    Serial.println(F("2 : equilibrium 2: upright, LQR balances theta = pi"));
    Serial.println(F("u : upright: swing-up helper on / off"));
    Serial.println(F("***"));
}

// Applies the keys read by Pendulum_lar::update()
void handleCommands(void)
{
    char c;
    while ((c = pendulum.readCommand()) != 0)
    {
        switch (c)
        {
        case '1':
        case '2':
            equilibrium = c - '0';
            swing.reset();
            pendulum.resetMode(Pendulum_lar::SWING_UP);
            updateLqr();
            printState();
            break;

        case 'u':
            use_swing_up = !use_swing_up;
            swing.reset();
            pendulum.resetMode(Pendulum_lar::SWING_UP);
            printState();
            break;

        case 'h':
            printHelp();
            break;

        default:
            break;
        }
    }
}

void setup()
{
    Serial.begin(115200);
    Serial.println();
    Serial.println(F("#Lab 3: LQR balance (torque)"));

    // TODO (Task 3.2c): start the motor and the pendulum, in this order (use the defines at the top; Lab 2 did it):
    //   1. The motor with the stored alignment of Lab 1, so that the pendulum does not move:
    //        foc.begin(supply voltage, voltage limit, align voltage, stored zero angle, stored direction)
    //      It returns false if the motor or the sensor does not answer: then stop here with  while (1);
    //   2. The calibration of Lab 1: foc.setTorqueConstant(k) and foc.setBackEmfConstant(Ke); Ke only when
    //      USE_BACK_EMF_COMP is 1 (otherwise 0: no back-EMF compensation).
    //   3. The pendulum: pendulum.begin(foc, control period, velocity filter alpha, run time).
    //   4. Its model: pendulum.setParameters(rod mass, point mass, length, b, tau_c). b is VISCOUS_FRICTION with the
    //      back-EMF compensation; without it the motor adds the electrical damping TORQUE_CONSTANT * MOTOR_KE.
    //   5. max_torque = the largest torque the motor can apply: foc.getMaxTorque() (= k * voltage limit).

    // given: the model and the torque limit (needs your step 4 and 5)
    Serial.print(F("#G = "));
    Serial.print(pendulum.gravityTorque() * 1000.0f, 2);
    Serial.print(F(" mN*m, max torque = "));
    Serial.print(max_torque * 1000.0f, 2);
    Serial.print(F(" mN*m, unstable pole = "));
    Serial.print(pendulum.naturalFrequency(), 2);
    Serial.println(F(" rad/s"));

    // TODO (Task 3.2c): the rest of the start-up:
    //   6. pendulum.setCoulombCompensation(USE_COULOMB_COMP != 0) and pendulum.setCaptureAngles(CAPTURE_ANGLE, RELEASE_ANGLE)
    //      (the thresholds of pendulum.updateMode()).
    //   7. The swing-up helper: swing.begin(pendulum, max_torque).
    //   8. The LQR gain of the initial equilibrium: updateLqr().
    //   9. If your Riccati solver did not converge (lqr_iterations == 0), stop here with  while (1);

    // given: what the sketch prints at power-up
    printState();
    printHelp();
}

void loop()
{
    // TODO (Task 3.2d): the control loop. loop() runs thousands of times per second; never use delay() here.
    //   1. foc.update(): the FOC algorithm. It must run on EVERY call of loop(), first.
    //   2. bool control_tick = pendulum.update(): reads the keys 's', 'x', 'p', keeps the run time and the timing,
    //      and returns true once per control period (while a run is active) with a fresh pendulum.theta, pendulum.omega.
    handleCommands(); // given: '1', '2', 'u', 'h'
    // TODO (Task 3.2d), continued:
    //   3. If a run has just started (pendulum.started()): swing.reset() and pendulum.resetMode(Pendulum_lar::SWING_UP).
    //   4. If there is no new sample (!control_tick), return: nothing to compute until the next control period.
    //   5. The angle error e to the equilibrium, wrapped to [-pi, pi) with pendulum_wrap_angle():
    //      - hanging (equilibrium == 1): e = theta - 0;
    //      - upright (equilibrium == 2): first call pendulum.updateMode(). While it returns Pendulum_lar::SWING_UP the
    //        pendulum is far from the top: apply the swing-up torque swing.torque(pendulum.theta, pendulum.omega) if
    //        use_swing_up is true, otherwise 0 (you lift the pendulum by hand), log with pendulum.log(0) and return.
    //        Otherwise (BALANCE) e = theta - pi (PENDULUM_LAR_PI).
    //   6. The LQR law tau = -(K[0]*e + K[1]*omega), limited to +/- max_torque with pendulum_clamp(); apply it with
    //      pendulum.applyTorque() and log the sample with pendulum.log(1) (mode 1 = balance).
}
