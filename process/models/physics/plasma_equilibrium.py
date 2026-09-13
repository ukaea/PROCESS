import os
from pathlib import Path
from typing import Any

import numpy as np
import scipy as sp

import process.core.constants as constants
from process.core.data_structure.base import DataStructure
from process.core.model import Model
from process.models.physics.profiles import NeProfile, TeProfile, PlasmaProfileShapeType

_veqpy: Any = None


def _import_veqpy():
    """Import veqpy on demand; site-packages install needs a writable numba cache dir."""
    global _veqpy
    if _veqpy is None:
        numba_cache = Path(__file__).resolve().parent / ".numba_cache"
        numba_cache.mkdir(exist_ok=True)
        os.environ.setdefault("NUMBA_CACHE_DIR", str(numba_cache))
        import veqpy as veq

        _veqpy = veq
    return _veqpy


class PlasmaEquilibrium(Model):
    def __init__(self):
        self.eq = None
        self.te_profile = TeProfile()
        self.ne_profile = NeProfile()
        self.ne_axis = 0.0
        self.te_axis = 0.0

    def _bind_profile_data(self):
        self.te_profile.data = self.data
        self.ne_profile.data = self.data

    def run(self):
        """PlasmaEquilibrium is invoked on demand, not in the main model loop."""

    def output(self):
        """PlasmaEquilibrium has no dedicated output block."""

    @staticmethod
    def solve_equilibrium_veqpy(
        pres_thermal_plasma_array,
        rho_pres_array,
        j_toroidal_array,
        rho_j_array,
        current_plasma_total_ampere,
        rminor,
        rmajor,
        kappa,
        delta,
        b_toroidal_rmajor,
    ):
        veq = _import_veqpy()
        nrho = 51
        rho = np.linspace(0.0, 1.0, nrho)
        p = np.interp(rho, rho_pres_array, pres_thermal_plasma_array)
        jtor = np.interp(rho, rho_j_array, j_toroidal_array)

        # PJ1 + coordinate=rho + nodes=uniform：psin 由源积分得到，不能再开 active psin 剖面
        topology = veq.KernelTopology(
            h_count=3,
            v_count=0,
            kappa_count=6,
            psin_count=0,
            F_count=0,
            c_counts=(),
            s_counts=(3,),
            Nr=16,
            Nt=16,
            route="PJ1",
            coordinate="rho",
            nodes="uniform",
            constraint="ip",
            sample_count=nrho,
        )

        kernel = veq.build(
            topology=topology,
            recipe=veq.KernelRecipe(backend="numba"),
            config=veq.KernelConfig(initial="cold"),
        )
        boundary = veq.KernelBoundary(
            a=rminor,
            R0=rmajor,
            Z0=0.0,
            B0=b_toroidal_rmajor,
            ka=kappa,
            s_offsets=(float(np.arcsin(np.clip(delta, -1.0, 1.0))),),
        )
        source = veq.KernelSource(p=p, jtor=jtor, Ip=float(current_plasma_total_ampere))
        result = kernel.solve(boundary=boundary, source=source)
        if not result.success:
            raise RuntimeError(
                f"veqpy equilibrium solve failed: residual={result.raw_norm:.3e}, nfev={result.nfev}"
            )
        return kernel.build_equilibrium()

    @staticmethod
    def get_volume_average(eq, rho_profile, profile):
        """Volume average of a 1D profile f(rho) using equilibrium geometry."""
        profile_on_grid = np.interp(eq.rho, rho_profile, profile)
        integrand = profile_on_grid[:, None] * eq.R * eq.J
        volume = eq.grid.integrate(eq.R * eq.J)
        return float(eq.grid.integrate(integrand) / volume)

    @staticmethod
    def get_volume_average_2(eq, profile):
        fun = np.asarray(profile, dtype=np.float64) * eq.R * eq.J
        vol = eq.grid.integrate(eq.R * eq.J)
        return float(eq.grid.integrate(fun) / vol)

    def pres_profile(
        self,
        alphan,
        alphat,
        tbeta,
        ne_ped,
        ne_sep,
        ne_axis,
        rho_nd_ped,
        f_pres_ie,
        te_ped,
        te_sep,
        te_axis,
        rho_temp_ped,
    ):
        rho = np.linspace(0.0, 1.0, 101)
        self.ne_profile.profile_y = np.zeros_like(rho)
        self.te_profile.profile_y = np.zeros_like(rho)
        self.ne_profile.calculate_profile_y(
            rho=rho,
            radius_plasma_pedestal_density_norm=rho_nd_ped,
            n0=ne_axis,
            nped=ne_ped,
            nsep=ne_sep,
            alphan=alphan,
        )

        self.te_profile.calculate_profile_y(
            rho=rho,
            radius_plasma_pedestal_temp_norm=rho_temp_ped,
            t0=te_axis,
            temp_plasma_pedestal_kev=te_ped,
            temp_plasma_separatrix_kev=te_sep,
            alphat=alphat,
            tbeta=tbeta,
        )
        pres_e_profile = (
            self.ne_profile.profile_y
            * self.te_profile.profile_y
            * constants.KILOELECTRON_VOLT
        )
        pres_i_profile = (
            pres_e_profile * f_pres_ie
        )
        return rho, (pres_e_profile + pres_i_profile)

    def iterate_equilibrium(self, ne_axis, te_axis, f_pres_ie):
        rho = np.linspace(0.0, 1.0, 101)
        j_toroidal_array = (1.0 - rho**2) ** self.data.physics.alphaj
        rho, pres_profile = self.pres_profile(
            alphan=self.data.physics.alphan,
            alphat=self.data.physics.alphat,
            tbeta=self.data.physics.tbeta,
            ne_ped=self.data.physics.nd_plasma_pedestal_electron,
            ne_sep=self.data.physics.nd_plasma_separatrix_electron,
            ne_axis=ne_axis,
            rho_nd_ped=self.data.physics.radius_plasma_pedestal_density_norm,
            f_pres_ie=f_pres_ie,
            te_ped=self.data.physics.temp_plasma_pedestal_kev,
            te_sep=self.data.physics.temp_plasma_separatrix_kev,
            te_axis=te_axis,
            rho_temp_ped=self.data.physics.radius_plasma_pedestal_temp_norm,
        )
        eq = self.solve_equilibrium_veqpy(
            pres_thermal_plasma_array=pres_profile,
            rho_pres_array=rho,
            j_toroidal_array=j_toroidal_array,
            rho_j_array=rho,
            current_plasma_total_ampere=self.data.physics.plasma_current,
            rminor=self.data.physics.rminor,
            rmajor=self.data.physics.rmajor,
            kappa=self.data.physics.kappa,
            delta=self.data.physics.triang,
            b_toroidal_rmajor=self.data.physics.b_plasma_toroidal_on_axis,
        )
        self.eq = eq
        ne_vol_avg = self.get_volume_average(eq, rho, self.ne_profile.profile_y)
        te_vol_avg = self.get_volume_average(eq, rho, self.te_profile.profile_y)
        return ne_vol_avg, te_vol_avg, eq

    def solve_axis_for_volume_averages(
        self,
        f_pres_ie,
        ne_vol_avg_target=None,
        te_vol_avg_target=None,
        ne_axis=None,
        te_axis=None,
        tol=1.0e-3,
        max_iter=25,
    ):
        """Iterate ``ne_axis`` and ``te_axis`` until volume averages match targets.

        Returns
        -------
        tuple[float, float, float, float, object]
            (ne_axis, te_axis, ne_vol_avg, te_vol_avg, eq)
        """
        self._bind_profile_data()
        physics = self.data.physics
        eq = None
        if ne_vol_avg_target is None:
            ne_vol_avg_target = physics.nd_plasma_electrons_vol_avg
        if te_vol_avg_target is None:
            te_vol_avg_target = physics.temp_plasma_electron_vol_avg_kev

        if ne_axis is None:
            ne_axis = NeProfile.ncore(
                physics.radius_plasma_pedestal_density_norm,
                physics.nd_plasma_pedestal_electron,
                physics.nd_plasma_separatrix_electron,
                ne_vol_avg_target,
                physics.alphan,
            )
        if te_axis is None:
            te_axis = TeProfile.tcore(
                physics.radius_plasma_pedestal_temp_norm,
                physics.temp_plasma_pedestal_kev,
                physics.temp_plasma_separatrix_kev,
                te_vol_avg_target,
                physics.alphat,
                physics.tbeta,
            )

        for _ in range(max_iter):
            ne_vol_avg, te_vol_avg, eq = self.iterate_equilibrium(ne_axis, te_axis, f_pres_ie)
            ne_err = abs(ne_vol_avg - ne_vol_avg_target) / ne_vol_avg_target
            te_err = abs(te_vol_avg - te_vol_avg_target) / te_vol_avg_target
            if ne_err <= tol and te_err <= tol:
                self.ne_axis = ne_axis
                self.te_axis = te_axis
                self.eq = eq
                return ne_axis, te_axis, ne_vol_avg, te_vol_avg, eq

            if ne_vol_avg > 0.0:
                ne_axis *= ne_vol_avg_target / ne_vol_avg
            if te_vol_avg > 0.0:
                te_axis *= te_vol_avg_target / te_vol_avg

        raise RuntimeError(
            "Failed to converge ne_axis/te_axis to target volume averages: "
            f"ne={ne_vol_avg:.4e} (target {ne_vol_avg_target:.4e}), "
            f"te={te_vol_avg:.4f} (target {te_vol_avg_target:.4f})"
        )

    def calculate_ind_plasma_internal_norm(self, eq, b_poloidal_avg):
        r_t = eq.surface_fields[2]   # ∂R/∂θ
        z_t = eq.Z_t                 # ∂Z/∂θ
        bp2 = (eq.alpha2 * eq.psin_r[:, None])**2 * (r_t**2 + z_t**2) / (eq.J * eq.R)**2
        # 体积平均 <Bp^2> = ∫ Bp^2 J R dθ dρ / V
        bp2_vol_avg = self.get_volume_average_2(eq, bp2)
        li = bp2_vol_avg / b_poloidal_avg**2
        return li

    @staticmethod
    def miller_delta_profile(eq) -> np.ndarray:
        R = np.asarray(eq.R, dtype=np.float64)
        Z = np.asarray(eq.Z, dtype=np.float64)
        nrho = R.shape[0]
        delta = np.zeros(nrho, dtype=np.float64)
        for i in range(nrho):
            Ri = R[i]
            Zi = Z[i]
            a_loc = eq.rho[i] * eq.a
            if a_loc < 1.0e-12:
                continue
            r_geo = eq.Rc[i]
            r_top = float(Ri[np.argmax(Zi)])
            delta[i] = (r_geo - r_top) / a_loc
        return delta