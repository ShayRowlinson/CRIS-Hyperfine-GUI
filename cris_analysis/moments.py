"""Nuclear moments from hyperfine constants, relative to a reference isotope.

    mu = mu_ref * (A * I) / (A_ref * I_ref)
    Q  = Q_ref * B / B_ref

Inputs can be plain floats or ``uncertainties.ufloat`` values; errors propagate.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


def parse_spin(value) -> float:
    """Turn ``"9/2"``, ``4.5`` or ``"4.5"`` into a float."""
    if isinstance(value, (int, float)):
        return float(value)
    return float(Fraction(str(value).strip()))


def magnetic_moment(A, I, A_ref, I_ref, mu_ref):
    if A_ref == 0 or I_ref == 0:
        raise ZeroDivisionError("Reference A and I must be non-zero to scale the magnetic moment")
    return mu_ref * (A * I) / (A_ref * I_ref)


def quadrupole_moment(B, B_ref, Q_ref):
    if B_ref == 0:
        raise ZeroDivisionError("Reference B must be non-zero to scale the quadrupole moment")
    return Q_ref * B / B_ref


@dataclass
class Moments:
    mu_l: object | None = None
    mu_u: object | None = None
    Q_l: object | None = None
    Q_u: object | None = None


def compute_moments(scan: dict, reference: dict) -> Moments:
    """Moments from both atomic states.

    ``scan`` needs I, A_l, A_u, B_l, B_u. ``reference`` needs I, mu, Q, A_l, A_u, B_l, B_u.
    A moment is left as None when its reference constant is zero, instead of v3's
    behaviour of refusing to compute anything (which blocked mu whenever B_ref was 0,
    e.g. for an I = 1/2 reference).
    """
    I, I_ref = parse_spin(scan["I"]), parse_spin(reference["I"])
    out = Moments()
    for state in ("l", "u"):
        A_ref, B_ref = reference[f"A_{state}"], reference[f"B_{state}"]
        if A_ref != 0 and I_ref != 0:
            setattr(out, f"mu_{state}", magnetic_moment(scan[f"A_{state}"], I, A_ref, I_ref, reference["mu"]))
        if B_ref != 0:
            setattr(out, f"Q_{state}", quadrupole_moment(scan[f"B_{state}"], B_ref, reference["Q"]))
    return out
