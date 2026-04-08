import parameters as p
import numpy as np

def compute_u_vortex(height_cyl):
    """Return U_vortex for a given cylinder height (m)."""
    if height_cyl <= 0:
        raise ValueError("height_cyl must be > 0")

    r_landeau = ((3 / 4) * p.R_cyl**2 * height_cyl) ** (1 / 3)
    phi_local = (p.m_0 / p.rho_p) / (np.pi * p.R_cyl**2 * height_cyl)
    rho_clouds_local = phi_local * p.rho_p + (1 - phi_local) * p.rho_released

    return np.sqrt(
        p.a2_Landeau
        * (rho_clouds_local - p.rho_infinity)
        / rho_clouds_local
        * p.g
        * r_landeau
        + p.a3_Landeau * p.g * r_landeau
    )

print(compute_u_vortex(0.027))

def solve_h_cyl_for_u_target(u_target, h_min=1e-10, h_max=100.0, tol=1e-10, max_iter=200):
    """Solve U_vortex(H_cyl)=u_target by bisection in [h_min, h_max]."""
    if u_target <= 0:
        raise ValueError("u_target must be > 0")
    if h_min <= 0 or h_max <= h_min:
        raise ValueError("Use bounds such that 0 < h_min < h_max")

    def residual(height):
        return compute_u_vortex(height) - u_target

    f_min = residual(h_min)
    f_max = residual(h_max)

    if f_min == 0:
        return h_min
    if f_max == 0:
        return h_max
    if f_min * f_max > 0:
        raise ValueError(
            "Target U_vortex is outside solver bounds. "
            "Increase h_max or reduce h_min."
        )

    left, right = h_min, h_max
    for _ in range(max_iter):
        mid = 0.5 * (left + right)
        f_mid = residual(mid)

        if abs(f_mid) < tol or (right - left) < tol:
            return mid

        if f_min * f_mid < 0:
            right = mid
            f_max = f_mid
        else:
            left = mid
            f_min = f_mid

    raise RuntimeError("Bisection did not converge within max_iter")



target_u_vortex_values = [0.5, 0.21, 0.1, 0.05, 0.01]

print("Target U_vortex (m/s) | H_cyl (m) | abs(U(H)-U_target)")
print("-" * 58)

for target_u in target_u_vortex_values:
    try:
        h_solution = solve_h_cyl_for_u_target(target_u)
        check_error = abs(compute_u_vortex(h_solution) - target_u)
        print(f"{target_u:>20.6f} | {h_solution:>9.6f} | {check_error:.3e}")
    except ValueError as exc:
        print(f"{target_u:>20.6f} | {'N/A':>9} | {exc}")
