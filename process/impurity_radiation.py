import dataclasses
import logging
import re
from importlib import resources
from pathlib import Path

import numpy as np
from numba import njit
from scipy import integrate

from process import constants
from process.data_structure import impurity_radiation_module
from process.exceptions import ProcessError, ProcessValueError

logger = logging.getLogger(__name__)


def initialise_imprad():
    """Initialises the impurity radiation data structure

    None
    This routine initialises the impurity radiation data.
    """

    errorflag = 0

    table_length = 200  # Number of temperature and Lz values in data file

    f_nd_species_electron = 1.0e0

    #  Hydrogen

    init_imp_element(
        n_species_index=1,
        name_label=impurity_radiation_module.imp_label[0],
        z=1,
        m_species_amu=constants.M_PROTIUM_AMU,  # 1.00782503223 1H
        f_nd_species_electron=f_nd_species_electron,
        len_tab=table_length,
        error=errorflag,
    )

    f_nd_species_electron = 0.0e0

    #  Helium
    init_imp_element(
        n_species_index=2,
        name_label=impurity_radiation_module.imp_label[1],
        z=2,
        m_species_amu=constants.M_HELIUM_AMU,  # 4.002602 (3He,4He) Average mass
        f_nd_species_electron=f_nd_species_electron,
        len_tab=table_length,
        error=errorflag,
    )

    #  Beryllium
    init_imp_element(
        n_species_index=3,
        name_label=impurity_radiation_module.imp_label[2],
        z=4,
        m_species_amu=constants.M_BERYLLIUM_AMU,  # 9.0121831 9Be
        f_nd_species_electron=f_nd_species_electron,
        len_tab=table_length,
        error=errorflag,
    )

    #  Carbon
    init_imp_element(
        n_species_index=4,
        name_label=impurity_radiation_module.imp_label[3],
        z=6,
        m_species_amu=constants.M_CARBON_AMU,  # 12.0096, (12C,13C,14C) Average mass
        f_nd_species_electron=f_nd_species_electron,
        len_tab=table_length,
        error=errorflag,
    )

    #  Nitrogen
    init_imp_element(
        n_species_index=5,
        name_label=impurity_radiation_module.imp_label[4],
        z=7,
        m_species_amu=constants.M_NITROGEN_AMU,  # 14.00643, (14N,15N) Average mass
        f_nd_species_electron=f_nd_species_electron,
        len_tab=table_length,
        error=errorflag,
    )

    #  Oxygen
    init_imp_element(
        n_species_index=6,
        name_label=impurity_radiation_module.imp_label[5],
        z=8,
        m_species_amu=constants.M_OXYGEN_AMU,  # 15.99903, (16O,17O,18O) Average mass
        f_nd_species_electron=f_nd_species_electron,
        len_tab=table_length,
        error=errorflag,
    )

    #  Neon
    init_imp_element(
        n_species_index=7,
        name_label=impurity_radiation_module.imp_label[6],
        z=10,
        m_species_amu=constants.M_NEON_AMU,  # 20.1797 (20Ne,21Ne,22Ne) Average mass
        f_nd_species_electron=f_nd_species_electron,
        len_tab=table_length,
        error=errorflag,
    )

    #  Silicon
    init_imp_element(
        n_species_index=8,
        name_label=impurity_radiation_module.imp_label[7],
        z=14,
        m_species_amu=constants.M_SILICON_AMU,  # 28.084 (28Si,29Si,30Si) Average mass
        f_nd_species_electron=f_nd_species_electron,
        len_tab=table_length,
        error=errorflag,
    )

    #  Argon
    init_imp_element(
        n_species_index=9,
        name_label=impurity_radiation_module.imp_label[8],
        z=18,
        m_species_amu=constants.M_ARGON_AMU,  # 39.948 (40Ar,36Ar,38Ar) Average mass
        f_nd_species_electron=f_nd_species_electron,
        len_tab=table_length,
        error=errorflag,
    )

    #  Iron
    init_imp_element(
        n_species_index=10,
        name_label=impurity_radiation_module.imp_label[9],
        z=26,
        m_species_amu=constants.M_IRON_AMU,  # 55.845 (56Fe,54Fe,57Fe,58Fe) Average mass
        f_nd_species_electron=f_nd_species_electron,
        len_tab=table_length,
        error=errorflag,
    )

    #  Nickel
    init_imp_element(
        n_species_index=11,
        name_label=impurity_radiation_module.imp_label[10],
        z=28,
        m_species_amu=constants.M_NICKEL_AMU,  # 58.6934 (58Ni,60Ni,61Ni,62Ni,64Ni) Average mass
        f_nd_species_electron=f_nd_species_electron,
        len_tab=table_length,
        error=errorflag,
    )

    #  Krypton
    init_imp_element(
        n_species_index=12,
        name_label=impurity_radiation_module.imp_label[11],
        z=36,
        m_species_amu=constants.M_KRYPTON_AMU,  # 83.798 (84Kr,86Kr,82Kr,80Kr,78Kr) Average mass
        f_nd_species_electron=f_nd_species_electron,
        len_tab=table_length,
        error=errorflag,
    )

    #  Xenon
    init_imp_element(
        n_species_index=13,
        name_label=impurity_radiation_module.imp_label[12],
        z=54,
        m_species_amu=constants.M_XENON_AMU,  # 131.293 (132Xe,129Xe,131Xe,134Xe,136Xe) Average mass
        f_nd_species_electron=f_nd_species_electron,
        len_tab=table_length,
        error=errorflag,
    )

    #  Tungsten
    init_imp_element(
        n_species_index=14,
        name_label=impurity_radiation_module.imp_label[13],
        z=74,
        m_species_amu=constants.M_TUNGSTEN_AMU,  # 183.84 (184W,186W,182W,183W,180W) Average mass
        f_nd_species_electron=f_nd_species_electron,
        len_tab=table_length,
        error=errorflag,
    )


@dataclasses.dataclass
class ImpurityDataHeader:
    """Represents a header or metadata section of an impurity data
    file.

    If this is a header for some section of data then the data section
    will be populated with the array of data for which this is a
    header of.
    """

    content: str
    data: list[float] | None = None


def read_impurity_file(impurity_file: Path):
    with open(impurity_file) as f:
        data = f.readlines()

    file_contents: list[ImpurityDataHeader] = []

    for line in data:
        # do not parse comments
        clean_line = line.strip().replace("\n", "")
        if clean_line[0:3].upper() in ["C  ", "C ", "C", "C--"]:
            continue

        if re.fullmatch(r"[0-9\.e+\- ]+", clean_line) is not None:
            header = file_contents[-1]

            new_data = clean_line.split(" ")
            if header.data is None:
                header.data = new_data
            else:
                header.data += new_data
        else:
            file_contents.append(ImpurityDataHeader(clean_line))

    return file_contents


def init_imp_element(
    n_species_index: int,
    name_label: str,
    z: int,
    m_species_amu: float,
    f_nd_species_electron: float,
    len_tab: int,
    error: int,
) -> None:
    """
    Initialise the impurity radiation data for a species.

    :param n_species_index: Position of species in impurity array
    :type n_species_index: int
    :param name_label: Species name
    :type name_label: str
    :param z: Species charge number
    :type z: int
    :param m_species_amu: Species atomic mass (amu)
    :type m_species_amu: float
    :param f_nd_species_electron: Number density / electron density
    :type f_nd_species_electron: float
    :param len_tab: Length of temperature and Lz tables
    :type len_tab: int
    :param error: Error flag; 0 = okay, 1 = missing impurity data
    :type error: int
    :raises ProcessValueError: If illegal impurity number is provided
    :raises FileNotFoundError: If impurity data files are missing
    :raises ProcessError: If required data cannot be located in files

    This routine initialises the impurity radiation data structure
    for a given impurity species. The Lz versus temperature data are
    read in from file.
    """

    if error == 1:
        return

    if n_species_index > len(impurity_radiation_module.impurity_arr_label):
        raise ProcessValueError(
            "Illegal impurity number",
            number=n_species_index,
            max=len(impurity_radiation_module.impurity_arr_label),
        )

    impurity_radiation_module.impurity_arr_label[n_species_index - 1] = name_label
    impurity_radiation_module.impurity_arr_z[n_species_index - 1] = z
    impurity_radiation_module.m_impurity_amu_array[n_species_index - 1] = m_species_amu
    impurity_radiation_module.f_nd_impurity_electron_array[n_species_index - 1] = (
        f_nd_species_electron
    )
    impurity_radiation_module.impurity_arr_len_tab[n_species_index - 1] = len_tab

    if len_tab > 200:
        print(
            f"ERROR: len_tab is {len_tab} but has a maximum value of {impurity_radiation_module.all_array_hotfix_len}"
        )

    impurity_dir = resources.files("process") / "data/lz_non_corona_14_elements/"

    lz_file = impurity_dir / f"{name_label}_lz_tau.dat"
    z_file = impurity_dir / f"{name_label}_z_tau.dat"

    if not lz_file.exists() or not z_file.exists():
        raise FileNotFoundError(
            f"Cannot find one or both of the impurity datafiles: {lz_file}, {z_file}"
        )

    lz_data = read_impurity_file(lz_file)
    z_data = read_impurity_file(z_file)

    Te = None
    lz = None

    for header in lz_data:
        if "Te[eV]" in header.content:
            Te = np.asarray(header.data, dtype=float)

        if "infinite confinement" in header.content:
            lz = np.asarray(header.data, dtype=float)

    if Te is None:
        raise ProcessError(f"Cannot locate Te data in {lz_file}")
    if lz is None:
        raise ProcessError(
            f"Cannot locate Lz for infinite confinement data in {lz_file}"
        )

    zav = None
    for header in z_data:
        if "infinite confinement" in header.content:
            zav = np.asarray(header.data, dtype=float)

    if zav is None:
        raise ProcessError(
            f"Cannot locate Zav for infinite confinement data in {z_file}"
        )

    impurity_radiation_module.temp_impurity_keV_array[n_species_index - 1, :] = Te * 1e-3
    impurity_radiation_module.pden_impurity_lz_nd_temp_array[n_species_index - 1, :] = lz
    impurity_radiation_module.impurity_arr_zav[n_species_index - 1, :] = zav


def z2index(zimp):
    for i in range(len(impurity_radiation_module.impurity_arr_label)):
        if zimp == impurity_radiation_module.impurity_arr_z[i]:
            return i

    # Should only get here if there is a problem
    raise ProcessValueError(
        "Element with the given charge is not in the impurity array", zimp=zimp
    )


def fradcore(rho, radius_plasma_core_norm, f_p_plasma_core_rad_reduction):
    """Finds the fraction of radiation from the core that is subtracted in impurity radiation model.

    Parameters
    ----------
    rho : numpy.array
        normalised minor radius
    radius_plasma_core_norm : float
        normalised radius defining the 'core' region
    f_p_plasma_core_rad_reduction : float
        fraction of radiation from the core region

    Returns
    -------
    numpy.array
        fradcore - array filled with the f_p_plasma_core_rad_reduction
    """
    fradcore = np.zeros(len(rho))
    rho_mask = rho < radius_plasma_core_norm
    fradcore[rho_mask] = f_p_plasma_core_rad_reduction

    return fradcore


def zav_of_te(imp_element_index, teprofile):
    """Calculates electron temperature dependent average atomic number

    Parameters
    ----------
    imp_element_index : int
        Impurity element index
    teprofile : numpy.array
        temperature profile

    Returns
    -------
    numpy.array
        zav_of_te - electron temperature dependent average atomic number
    """
    # less_than_imp_temp_mask = teprofile values less than impurity temperature. greater_than_imp_temp_mask = teprofile values higher than impurity temperature.
    return _zav_of_te_compiled(
        imp_element_index,
        teprofile,
        impurity_radiation_module.temp_impurity_keV_array,
        impurity_radiation_module.impurity_arr_zav,
        impurity_radiation_module.impurity_arr_len_tab,
    )


@njit(cache=True)
def _zav_of_te_compiled(
    imp_element_index: int,
    teprofile: np.array,
    temp_impurity_keV_array: np.array,
    impurity_arr_zav: np.array,
    impurity_arr_len_tab: np.array,
):
    bins = temp_impurity_keV_array[imp_element_index]
    indices = np.digitize(teprofile, bins)
    indices[indices >= bins.shape[0]] = bins.shape[0] - 1
    indices[indices < 0] = 0
    # Use numpy.interp for linear interpolation in log space
    zav_of_te = np.interp(
        np.log(teprofile),
        np.log(temp_impurity_keV_array[imp_element_index, :]),
        impurity_arr_zav[imp_element_index, :],
    )

    less_than_imp_temp_mask = teprofile <= temp_impurity_keV_array[imp_element_index, 0]

    zav_of_te[less_than_imp_temp_mask] = impurity_arr_zav[imp_element_index, 0]
    greater_than_imp_temp_mask = (
        teprofile
        >= temp_impurity_keV_array[
            imp_element_index,
            (impurity_arr_len_tab[imp_element_index]) - 1,
        ]
    )
    zav_of_te[greater_than_imp_temp_mask] = impurity_arr_zav[
        imp_element_index,
        impurity_arr_len_tab[imp_element_index] - 1,
    ]

    return zav_of_te


def pimpden(imp_element_index, neprofile, teprofile):
    """Calculates the impurity radiation density (W/m3)

    Parameters
    ----------
    imp_element_index : int
        Impurity element index
    neprofile : numpy.array
        electron density profile
    teprofile : numpy.array
        electron temperature profile

    Returns
    -------
    numpy.array
        pimpden - total impurity radiation density (W/m3)
    """
    # less_than_imp_temp_mask = teprofile values less than impurity temperature. greater_than_imp_temp_mask = teprofile values higher than impurity temperature.
    bins = impurity_radiation_module.temp_impurity_keV_array[imp_element_index]
    indices = np.digitize(teprofile, bins)
    indices[indices >= bins.shape[0]] = bins.shape[0] - 1
    indices[indices < 0] = 0

    # Use numpy.interp for linear interpolation in log-log space
    pimpden = np.exp(
        np.interp(
            np.log(teprofile),
            np.log(
                impurity_radiation_module.temp_impurity_keV_array[imp_element_index, :]
            ),
            np.log(
                impurity_radiation_module.pden_impurity_lz_nd_temp_array[
                    imp_element_index, :
                ]
            ),
        )
    )

    pimpden = (
        impurity_radiation_module.f_nd_impurity_electron_array[imp_element_index]
        * neprofile
        * neprofile
        * pimpden
    )

    less_than_imp_temp_mask = (
        teprofile
        <= impurity_radiation_module.temp_impurity_keV_array[imp_element_index, 0]
    )
    pimpden[less_than_imp_temp_mask] = (
        impurity_radiation_module.pden_impurity_lz_nd_temp_array[imp_element_index, 0]
    )

    greater_than_imp_temp_mask = (
        teprofile
        >= impurity_radiation_module.temp_impurity_keV_array[
            imp_element_index,
            impurity_radiation_module.impurity_arr_len_tab[imp_element_index] - 1,
        ]
    )
    #  This is okay because Bremsstrahlung will dominate at higher temp.
    pimpden[greater_than_imp_temp_mask] = (
        impurity_radiation_module.pden_impurity_lz_nd_temp_array[
            imp_element_index,
            impurity_radiation_module.impurity_arr_len_tab[imp_element_index] - 1,
        ]
    )

    return pimpden


def element2index(element: str):
    """Returns the index of the `element` in the impurity array with
    a given name
    """
    try:
        return (
            impurity_radiation_module.impurity_arr_label.astype(str)
            .tolist()
            .index(element)
        )
    except ValueError as e:
        raise ProcessValueError(
            f"Element {element} is not found in impurity_arr_label"
        ) from e


class ImpurityRadiation:
    """This class calculates the impurity radiation losses for given temperature and density profiles.
    The considers the  total impurity radiation from the core (pden_impurity_core_rad_total_mw) and total impurity radiation
    (pden_impurity_rad_total_mw) [MW/(m^3)]. The class is used to sum the impurity radiation loss from each impurity
    element to find the total impurity radiation loss.
    """

    def __init__(self, plasma_profile):
        """
        :param plasma_profile: Plasma profile class, parameterises the density and temperature profiles.
        :type plasma_profile: Plasma profile class
        """
        self.plasma_profile = plasma_profile
        self.rho = plasma_profile.neprofile.profile_x
        self.rhodx = plasma_profile.neprofile.profile_dx
        self.imp = np.nonzero(
            impurity_radiation_module.f_nd_impurity_electron_array > 1.0e-30
        )[0]

        self.pimp_profile = np.zeros(self.plasma_profile.profile_size)
        self.pden_impurity_rad_profile = np.zeros(self.plasma_profile.profile_size)
        self.pden_impurity_core_rad_profile = np.zeros(self.plasma_profile.profile_size)

        self.pden_impurity_rad_total_mw = 0.0
        self.pden_impurity_core_rad_total_mw = 0.0

    def map_imprad_profile(self):
        """Map imprad_profile() over each impurity element index."""
        list(map(self.imprad_profile, self.imp))

    def imprad_profile(self, imp_element_index):
        """This routine calculates the impurity radiation losses for given temperature and density profiles.
        References:
            Bremsstrahlung equation from Johner, L(z) data (coronal equilibrium)
            from Marco Sertoli, Asdex-U, ref. Kallenbach et al.
            Johner, Fusion Science and Technology 59 (2011), pp 308-349
            Sertoli, private communication
            Kallenbach et al., Plasma Phys. Control. Fus. 55(2013) 124041
        Authors:
            R Kemp, CCFE, Culham Science Centre
            H Lux, CCFE, Culham Science Centre
            P J Knight, CCFE, Culham Science Centre
            G Turkington, CCFFE, Culham Science Centre

        Parameters
        ----------
        imp_element_index : Int
            Index used to access different impurity radiation elements
        """

        pimp = pimpden(
            imp_element_index,
            self.plasma_profile.neprofile.profile_y,
            self.plasma_profile.teprofile.profile_y,
        )

        self.pimp_profile = np.add(self.pimp_profile, pimp)

    def calculate_radiation_loss_profiles(self):
        """Calculate the Bremsstrahlung (radb), line radiation (radl), total impurity radiation
        from the core (pden_impurity_core_rad_total_mw) and total impurity radiation (pden_impurity_rad_total_mw). Update the stored arrays with
        the values.
        """

        pden_impurity_rad_total = self.pimp_profile * self.rho
        pden_impurity_core_rad_total = self.pimp_profile * (
            self.rho
            * fradcore(
                self.rho,
                impurity_radiation_module.radius_plasma_core_norm,
                impurity_radiation_module.f_p_plasma_core_rad_reduction,
            )
        )

        self.pden_impurity_rad_profile = np.add(
            self.pden_impurity_rad_profile, pden_impurity_rad_total
        )
        self.pden_impurity_core_rad_profile = np.add(
            self.pden_impurity_core_rad_profile, pden_impurity_core_rad_total
        )

    def integrate_radiation_loss_profiles(self):
        """Integrate the radiation loss profiles using the Simpson rule.
        Store the total values for each aspect of impurity radiation loss.
        """

        # 2.0e-6 converts from W/m^3 to MW/m^3 and also accounts for both sides of the plasma
        self.pden_impurity_rad_total_mw = 2.0e-6 * integrate.simpson(
            self.pden_impurity_rad_profile, x=self.rho, dx=self.rhodx
        )
        self.pden_impurity_core_rad_total_mw = 2.0e-6 * integrate.simpson(
            self.pden_impurity_core_rad_profile, x=self.rho, dx=self.rhodx
        )

    def calculate_imprad(self):
        """Call the map function to calculate impurity radiation parameters for each
        impurity element. Calculate the radiation loss profiles, and integrate them to
        find the total values for radiation loss.
        """
        self.map_imprad_profile()
        self.calculate_radiation_loss_profiles()
        self.integrate_radiation_loss_profiles()
