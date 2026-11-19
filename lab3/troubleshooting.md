# Lab 3 — Troubleshooting

- **DARE does not converge / closed-loop eigenvalue is outside the unit circle:** check the units and signs of `J`, `G`, `b`, the equilibrium selected in `EQUILIBRIUM`, and the sign convention for encoder/torque. Do not run the controller until the upright or hanging matrix is stable as appropriate.
- **Upright LQR does nothing from the bottom:** with `USE_SWING_UP 0` (key `u`) the motor stays at zero torque until `rig.updateMode()` reports `BALANCE`: lift the pendulum by hand close to the top. Press `u` for the energy swing-up.
- **The sketch does not compile (`foc`, `rig`, `K` ... not declared):** complete Task 3.2a, the global declarations, first.
- **LQR does not swing up from hanging:** this is expected when `tau_max < G`. LQR is a local balancing law, not the Lab 2 energy controller. Use the energy swing-up (`USE_SWING_UP 1`, key `u`) or bring the pendulum near the top by hand.
- **Motor/encoder jumps at startup:** redo the alignment procedure with the shaft free, then paste the stored zero angle and sensor direction from the alignment sketch.
- **Torque is lower than expected:** check that back-EMF compensation is enabled consistently, `TORQUE_CONSTANT` and `MOTOR_VOLTAGE_LIMIT` are measured correctly, and the supply voltage does not sag.
