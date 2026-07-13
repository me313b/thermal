"""correlations.py - correlations shared between the lumped solver
(app.py) and the zonal solver (zonal.py), so a cross-check between them
can never diverge on a correlation choice again (red-team F5).

Single source of truth for the internal water-side Nusselt number:
Hausen (laminar, entry-corrected) below Re = 2300, Gnielinski above
Re = 3000, with a linear blend across the 2300-3000 transition band so
the film coefficient is continuous through transition (no step).
"""
import math


def water_nu(Re, Pr, d_i, L):
    """Internal forced-convection Nusselt number for water in a round
    tube. Continuous across the laminar-turbulent transition.

    Re  Reynolds number (4 mdot / (pi mu d_i))
    Pr  Prandtl number
    d_i tube inner diameter [m]
    L   heated length [m]
    returns (Nu, regime)
    """
    def nu_lam(Re_):
        gz = (d_i / L) * Re_ * Pr
        return 3.66 + 0.0668 * gz / (1.0 + 0.04 * gz ** (2.0 / 3.0))

    def nu_turb(Re_):
        f = (0.790 * math.log(max(Re_, 3000.0)) - 1.64) ** -2
        return ((f / 8.0) * (Re_ - 1000.0) * Pr
                / (1.0 + 12.7 * math.sqrt(f / 8.0)
                   * (Pr ** (2.0 / 3.0) - 1.0)))

    if Re <= 0:
        return 3.66, "no flow"
    if Re < 2300.0:
        return nu_lam(Re), "laminar"
    if Re < 3000.0:
        w = (Re - 2300.0) / 700.0
        return (1.0 - w) * nu_lam(2300.0) + w * nu_turb(3000.0), \
            "transitional"
    return nu_turb(Re), "turbulent"
