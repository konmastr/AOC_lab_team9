# Lab 2 — Model calculations (handout §2)

Inputs: Lab 1 rig card (`../lab1/rig_card.md`). Design settings: handout §1.3–1.4 and §5.2.

| Input | Symbol | Value | Unit |
|---|---|---:|---|
| Inertia | $J$ | 0.459e-03 | kg·m² |
| Gravity-torque coefficient | $G$ | 0.052974 | N·m |
| Viscous friction | $b$ | 1.3147e-05 | N·m·s/rad |
| Torque constant | $k$ | 0.009340 | N·m/V |
| Voltage limit | $U_{\max}$ | 3.0 | V |

## 1. Model

| Quantity | Formula | Value | Unit |
|---|---|---:|---|
| Natural frequency | $\omega_0=\sqrt{G/J}$ | 10.743 | rad/s |
| Max torque | $\tau_{\max}=kU_{\max}$ | 0.02802 | N·m |
| Energy at the top | $E_{up}=2G$ | 0.10595 | J |

## 2. PID gains (torque gains, $\zeta = 0.7$)

$k_p=J\omega_n^2-q$, $k_d=2\zeta\omega_nJ-b$, $k_i=\ell\,k_p$, with $q=+G$ hanging, $q=-G$ upright.

| | Hanging | Upright | Unit |
|---|---:|---:|---|
| $\omega_n$ | 25 | 15 | rad/s |
| $\ell = k_i/k_p$ | 1.0 | 0.5 | 1/s |
| $k_p$ | **0.23390** | **0.15625** | N·m/rad |
| $k_d$ | **0.016052** | **0.0096259** | N·m·s/rad |
| $k_i$ | **0.23390** | **0.078125** | N·m/(rad·s) |

Hanging condition $\omega_n^2>G/J$: $625 > 115.4$ ✔ ($k_p>0$).

Sketch defines:

```cpp
#define PID_WN_HANGING 25.0f
#define PID_WN_UPRIGHT 15.0f
#define PID_ZETA 0.7f
#define PID_KI_OVER_KP_HANGING 1.0f
#define PID_KI_OVER_KP_UPRIGHT 0.5f

#define SWING_ENERGY_RATIO 1.03f
#define SWING_ENERGY_BAND 0.10f
#define SWING_KICK_SPEED 0.3f
```

`printGains()` at power-up should print these $k_p$, $k_i$, $k_d$ (up to float rounding).

## 3. Torque limit

| Quantity | Formula | Value |
|---|---|---:|
| Torque ratio | $\tau_{\max}/G$ | 0.529 (< 1: cannot lift directly) |
| Largest recoverable angle from the top (at rest) | $\varphi_c=\arcsin(\tau_{\max}/G)$ | **0.557 rad = 31.9°** |

Part 1, Test B: release the pendulum at rest clearly inside ±32° of the top (aim for ≲ 10–15°, see 5.3).

## 4. Swing-up (Part 2)

$\omega_{top}=\sqrt{4G(\rho-1)/J}$ with $\rho=1.03$: **3.72 rad/s** < 8 rad/s capture speed ✔

## 5. Extra checks (not in the handout)

1. **Closed-loop poles.** PD (design): hanging $-17.5\pm17.85j$, upright $-10.5\pm10.71j$.
   With $k_i$ (third-order $A_{PID}=Js^3+(b+k_d)s^2+(q+k_p)s+k_i$): hanging $-17.07\pm17.45j,\ -0.86$; upright $-10.09\pm10.33j,\ -0.82$.
   The slow real pole (~0.8 rad/s, time constant ~1.2 s) is the "slow tail" of the integral term.
2. **Routh bound on $k_i$:** stable iff $k_i<(b+k_d)(q+k_p)/J$ → hanging $k_i<10.0$, upright $k_i<2.17$. Chosen values are far inside ✔
3. **Capture vs. recoverable angle:** `CAPTURE_ANGLE` 0.4 rad < $\varphi_c$ = 0.557 rad ✔
4. **Saturation angle of the PD at rest**, $\tau_{\max}/k_p$: hanging 0.120 rad (6.9°), upright 0.179 rad (10.3°).
   Beyond these angles the output is at the limit and the response is no longer the linear design. At capture
   (0.4 rad) the upright PID saturates — this is why Part 2 needs the clamp and the anti-windup.
5. **Sensitivity to $k$.** The torque fit with offset (Lab 1 handout, Task 2.7) gives $k = 0.008939$ instead of
   the rig-card $k_0 = 0.009340$ (−4.5%). With it: $\tau_{\max} = 0.02682$ N·m, $\tau_{\max}/G = 0.506$,
   $\varphi_c = 0.531$ rad (30.4°). The gains above do not depend on $k$ (they are torque gains).
