| Item | Symbol / sketch setting | Value | Unit |
|---|---|---:|---|
| Pole pairs | `MOTOR_POLE_PAIRS` | 7 |  |
| Stored zero electrical angle | `STORED_ZERO_ELECTRIC_ANGLE` | 4.65717  | rad |
| Stored sensor direction | `STORED_SENSOR_DIRECTION` | CW |  |
| Torque constant | `TORQUE_CONSTANT`, $k$ | 0.009340 | N·m/V |
| Back-EMF constant, low-speed estimate | `MOTOR_KE`, $K_e$ | 0.02839 | V·s/rad |
| Rod mass | `ROD_MASS`, $m_r$ | 0.02 | kg |
| Point mass | `POINT_MASS`, $m_p$ | 0.05 | kg |
| Length, axis to point-mass centre | `PENDULUM_LENGTH`, $l$ | 0.09 | m |
| Inertia | $J$ | 0.459e-03 | kg·m² |
| Gravity-torque coefficient | $G$ | 0.052974 | N·m |
| Natural frequency | $\sqrt{G/J}$ | 10.743 | rad/s |
| Predicted period | $2\pi\sqrt{J/G}$ | 0.585 | s |
| Mechanical viscous friction | `VISCOUS_FRICTION`, $b$ | 1.3147e-05 | N·m·s/rad |
| Mechanical Coulomb friction | `COULOMB_FRICTION`, $\tau_c$ | 0.5263e-03 | N·m |
| Voltage limit planned for controller labs | `MOTOR_VOLTAGE_LIMIT` | **3 (suggested)** | V |
| Estimated maximum torque at 3 V | $\tau_{max}=3k$ | 0.02802 | N·m |
| Torque-to-gravity ratio | $\tau_{max}/G$ | 0.5289 | — |