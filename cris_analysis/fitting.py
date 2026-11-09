"""Hyperfine-structure fitting with SATLAS2.

Fits counts per bunch against frequency with one or more hyperfine components (ground
state, isomers, contaminants) plus a polynomial background.

Differences from v3, which used the deprecated ``satlas2.interface.HFSModel``:

* v4 builds ``satlas2.HFS`` directly. The interface wrapper silently added a second,
  fixed background (its ``background_params``) on top of the fitted one, so the fitted
  background reported by v3 was offset by the starting value. A, B and centroid are
  unaffected.
* Component names are unique, so two components with the same spin no longer collide.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import satlas2

SOURCE_NAME = "Fit"


@dataclass
class HFSComponent:
    """Starting values and constraints for one hyperfine spectrum."""

    I: float
    J_l: float
    J_u: float
    A_l: float = 0.0
    A_u: float = 0.0
    B_l: float = 0.0
    B_u: float = 0.0
    centroid: float = 0.0
    scale: float = 1.0
    fwhm_g: float = 50.0
    fwhm_l: float = 50.0
    racah: bool = True
    fix: set[str] = field(default_factory=set)  # any of "Al", "Au", "Bl", "Bu", "centroid", "FWHMG", ...
    A_ratio: float | None = None  # if set, Au = A_ratio * Al
    B_ratio: float | None = None  # if set, Bu = B_ratio * Bl
    name: str | None = None


@dataclass
class ParameterValue:
    value: float
    uncertainty: float | None
    vary: bool


@dataclass
class FitResult:
    parameters: dict[str, dict[str, ParameterValue]]
    chisqr: float
    redchi: float
    nfree: int
    ndata: int
    report: str
    evaluate: Callable[[np.ndarray], np.ndarray]
    uncertainties_estimated: bool = True

    def table(self) -> pd.DataFrame:
        rows = [
            {"component": comp, "parameter": name, "value": p.value, "uncertainty": p.uncertainty, "vary": p.vary}
            for comp, params in self.parameters.items()
            for name, p in params.items()
        ]
        return pd.DataFrame(rows)


def _build_model(comp: HFSComponent, name: str) -> satlas2.HFS:
    model = satlas2.HFS(
        comp.I,
        [comp.J_l, comp.J_u],
        A=[comp.A_l, comp.A_u],
        B=[comp.B_l, comp.B_u],
        C=[0, 0],
        df=comp.centroid,
        fwhmg=comp.fwhm_g,
        fwhml=comp.fwhm_l,
        scale=comp.scale,
        racah=comp.racah,
        name=name,
    )
    for key in ("Cl", "Cu"):
        model.params[key].vary = False
    # B is undefined for J < 1 (and A for J = 0). Leaving it free makes the covariance
    # matrix singular, which is why v3 often reported "uncertainties could not be estimated".
    for state, J in (("l", comp.J_l), ("u", comp.J_u)):
        if J < 1:
            model.params["B" + state].value = 0.0
            model.params["B" + state].vary = False
        if J == 0:
            model.params["A" + state].value = 0.0
            model.params["A" + state].vary = False
    for key in comp.fix:
        if key not in model.params:
            raise KeyError(f"Unknown parameter '{key}' for component {name}")
        model.params[key].vary = False
    for upper, lower, ratio in (("Au", "Al", comp.A_ratio), ("Bu", "Bl", comp.B_ratio)):
        if ratio is None:
            continue
        if model.params[lower].value == 0 and not model.params[lower].vary:
            model.params[upper].vary = False
        else:
            model.params[upper].expr = f"{ratio}*{SOURCE_NAME}___{name}___{lower}"
            model.params[upper].vary = False
    return model


def spectrum_arrays(spectrum: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """x, counts per bunch and its error, as v3 fitted them."""
    bunches = spectrum["bunches"].to_numpy(dtype=float)
    return (
        spectrum["x"].to_numpy(dtype=float),
        spectrum["y"].to_numpy(dtype=float) / bunches,
        spectrum["yerr"].to_numpy(dtype=float) / bunches,
    )


def fit_spectrum(
    spectrum: pd.DataFrame,
    components: list[HFSComponent],
    background: list[float] | float = 0.001,
    **fit_kwargs,
) -> FitResult:
    """Fit ``components`` plus a polynomial background to a binned spectrum.

    ``background`` holds the starting polynomial coefficients, highest order first
    (a single number is a flat background). Extra keyword arguments go to
    ``satlas2.Fitter.fit``.
    """
    if not components:
        raise ValueError("Give at least one HFSComponent")
    x, y, yerr = spectrum_arrays(spectrum)
    source = satlas2.Source(x, y, yerr, name=SOURCE_NAME)

    models = {}
    for i, comp in enumerate(components, start=1):
        name = comp.name or f"hfs{i}"
        models[name] = _build_model(comp, name)
        source.addModel(models[name])
    bkg = satlas2.Polynomial(np.atleast_1d(background).tolist(), name="bkg")
    source.addModel(bkg)

    fitter = satlas2.Fitter()
    fitter.addSource(source)
    fitter.fit(**fit_kwargs)
    report = fitter.reportFit()

    params = {
        name: {k: ParameterValue(p.value, p.unc, p.vary) for k, p in model.params.items()}
        for name, model in {**models, "bkg": bkg}.items()
    }
    return FitResult(
        parameters=params,
        chisqr=float(fitter.chisqr),
        redchi=float(fitter.redchi),
        nfree=int(fitter.nfree),
        ndata=int(fitter.ndata),
        report=report,
        evaluate=source.evaluate,
        uncertainties_estimated="uncertainties could not be estimated" not in report,
    )
