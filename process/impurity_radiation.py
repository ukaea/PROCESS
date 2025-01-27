import dataclasses
import logging
import re
from importlib import resources
from pathlib import Path

import numpy as np
from scipy import integrate

from process.fortran import error_handling, impurity_radiation_module

logger = logging.getLogger(__name__)


def initialise_imprad():
    """
    Initialises the impurity radiation data structure
    author: H Lux, CCFE, Culham Science Centre
    None
    This routine initialises the impurity radiation data.
    """

    errorflag = 0

    table_length = 200  # Number of temperature and Lz values in data file

    frac = 1.0e0

    #  Hydrogen

    init_imp_element(
        no=1,
        label=impurity_radiation_module.imp_label[0],
        z=1,
        amass=1.01e0,
        frac=frac,
        len_tab=table_length,
        error=errorflag,
    )

    frac = 0.0e0

    #  Helium
    init_imp_element(
        2,
        impurity_radiation_module.imp_label[1],
        2,
        4.003e0,
        frac,
        table_length,
        errorflag,
    )

    #  Beryllium
    init_imp_element(
        3,
        impurity_radiation_module.imp_label[2],
        4,
        9.01e0,
        frac,
        table_length,
        errorflag,
    )

    #  Carbon
    init_imp_element(
        4,
        impurity_radiation_module.imp_label[3],
        6,
        12.01e0,
        frac,
        table_length,
        errorflag,
    )

    #  Nitrogen
    init_imp_element(
        5,
        impurity_radiation_module.imp_label[4],
        7,
        14.01e0,
        frac,
        table_length,
        errorflag,
    )

    #  Oxygen
    init_imp_element(
        6,
        impurity_radiation_module.imp_label[5],
        8,
        15.999e0,
        frac,
        table_length,
        errorflag,
    )

    #  Neon
    init_imp_element(
        7,
        impurity_radiation_module.imp_label[6],
        10,
        20.18e0,
        frac,
        table_length,
        errorflag,
    )

    #  Silicon
    init_imp_element(
        8,
        impurity_radiation_module.imp_label[7],
        14,
        28.09e0,
        frac,
        table_length,
        errorflag,
    )

    #  Argon
    init_imp_element(
        9,
        impurity_radiation_module.imp_label[8],
        18,
        39.95e0,
        frac,
        table_length,
        errorflag,
    )

    #  Iron
    init_imp_element(
        10,
        impurity_radiation_module.imp_label[9],
        26,
        55.85e0,
        frac,
        table_length,
        errorflag,
    )

    #  Nickel
    init_imp_element(
        11,
        impurity_radiation_module.imp_label[10],
        28,
        58.70e0,
        frac,
        table_length,
        errorflag,
    )

    #  Krypton
    init_imp_element(
        12,
        impurity_radiation_module.imp_label[11],
        36,
        83.80e0,
        frac,
        table_length,
        errorflag,
    )

    #  Xenon
    init_imp_element(
        13,
        impurity_radiation_module.imp_label[12],
        54,
        131.30e0,
        frac,
        table_length,
        errorflag,
    )

    #  Tungsten
    init_imp_element(
        14,
        impurity_radiation_module.imp_label[13],
        74,
        183.85e0,
        frac,
        table_length,
        errorflag,
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


def init_imp_element(no, label, z, amass, frac, len_tab, error):
    """
    Initialises the impurity radiation data for a species
    author: H Lux, CCFE, Culham Science Centre
    author: P J Knight, CCFE, Culham Science Centre
    no      : input integer  : position of species in impurity array
    label   : input string   : species name
    Z       : input integer  : species charge number
    amass   : input real     : species atomic mass (amu)
    frac    : input real     : number density / electron density
    len_tab : input integer  : length of temperature and Lz tables
    error   : input/output integer : Error flag; 0 = okay, 1 = missing
    impurity data
    This routine initialises the impurity radiation data structure
    for a given impurity species.
    The Lz versus temperature data are read in from file.
    """

    if error == 1:
        return

    if no > len(impurity_radiation_module.impurity_arr_label):
        error_handling.idiags[0] = no
        error_handling.idiags[1] = len(impurity_radiation_module.impurity_arr_label)
        error_handling.report_error(27)

    impurity_radiation_module.impurity_arr_label[no - 1] = label
    impurity_radiation_module.impurity_arr_z[no - 1] = z
    impurity_radiation_module.impurity_arr_amass[no - 1] = amass
    impurity_radiation_module.impurity_arr_frac[no - 1] = frac
    impurity_radiation_module.impurity_arr_len_tab[no - 1] = len_tab

    if len_tab > 200:
        print(
            f"ERROR: len_tab is {len_tab} but has a maximum value of {impurity_radiation_module.all_array_hotfix_len}"
        )

    impurity_label = label.decode("utf-8")
    impurity_dir = resources.files("process") / "data/lz_non_corona_14_elements/"

    lz_file = impurity_dir / f"{impurity_label}_lz_tau.dat"
    z_file = impurity_dir / f"{impurity_label}_z_tau.dat"

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
        raise RuntimeError(f"Cannot locate Te data in {lz_file}")
    if lz is None:
        raise RuntimeError(
            f"Cannot locate Lz for infinite confinement data in {lz_file}"
        )

    zav = None
    for header in z_data:
        if "infinite confinement" in header.content:
            zav = np.asarray(header.data, dtype=float)

    if zav is None:
        raise RuntimeError(
            f"Cannot locate Zav for infinite confinement data in {z_file}"
        )

    impurity_radiation_module.impurity_arr_temp_kev[no - 1, :] = Te * 1e-3
    impurity_radiation_module.impurity_arr_lz_wm3[no - 1, :] = lz
    impurity_radiation_module.impurity_arr_zav[no - 1, :] = zav


def z2index(zimp):
    for i in range(len(impurity_radiation_module.impurity_arr_label)):
        if zimp == impurity_radiation_module.impurity_arr_z[i]:
            return i

    # Should only get here if there is a problem

    error_handling.idiags[0] = zimp
    error_handling.report_error(33)

    # Explicit return
    return None


def fradcore(rho, coreradius, coreradiationfraction):
    """Finds the fraction of radiation from the core that is subtracted in impurity radiation model.

    :param rho: normalised minor radius
    :type rho: numpy.array
    :param coreradius: normalised radius defining the 'core' region
    :type coreradius: float
    :param coreradiationfraction: fraction of radiation from the core region
    :type coreradiationfraction: float
    :return: fradcore - array filled with the coreradiationfraction
    :rtype: numpy.array
    """
    fradcore = np.zeros(len(rho))
    rho_mask = rho < coreradius
    fradcore[rho_mask] = coreradiationfraction

    return fradcore


def zav_of_te(imp_element_index, tprofile):
    """Calculates electron temperature dependent average atomic number

    :param imp_element_index: Impurity element index
    :type imp_element_index: int
    :param tprofile: temperature profile
    :type tprofile: numpy.array
    :return: zav_of_te - electron temperature dependent average atomic number
    :rtype: numpy.array
    """
    # less_than_imp_temp_mask = tprofile values less than impurity temperature. greater_than_imp_temp_mask = tprofile values higher than impurity temperature.
    bins = impurity_radiation_module.impurity_arr_temp_kev[imp_element_index]
    indices = np.digitize(tprofile, bins)
    indices[indices >= bins.shape[0]] = bins.shape[0] - 1
    indices[indices < 0] = 0
    yi = impurity_radiation_module.impurity_arr_zav[imp_element_index, indices - 1]
    xi = np.log(
        impurity_radiation_module.impurity_arr_temp_kev[imp_element_index, indices - 1]
    )
    c = (
        impurity_radiation_module.impurity_arr_zav[imp_element_index, indices] - yi
    ) / (
        np.log(
            impurity_radiation_module.impurity_arr_temp_kev[imp_element_index, indices]
        )
        - xi
    )
    zav_of_te = yi + c * (np.log(tprofile) - xi)
    less_than_imp_temp_mask = (
        tprofile
        <= impurity_radiation_module.impurity_arr_temp_kev[imp_element_index, 0]
    )
    zav_of_te[less_than_imp_temp_mask] = impurity_radiation_module.impurity_arr_zav[
        imp_element_index, 0
    ]
    greater_than_imp_temp_mask = (
        tprofile
        >= impurity_radiation_module.impurity_arr_temp_kev[
            imp_element_index,
            (impurity_radiation_module.impurity_arr_len_tab[imp_element_index]) - 1,
        ]
    )
    zav_of_te[greater_than_imp_temp_mask] = impurity_radiation_module.impurity_arr_zav[
        imp_element_index,
        impurity_radiation_module.impurity_arr_len_tab[imp_element_index] - 1,
    ]

    return zav_of_te


def pimpden(imp_element_index, nprofile, tprofile):
    """Calculates the impurity radiation density (W/m3)

    :param imp_element_index: Impurity element index
    :type imp_element_index: int
    :param nprofile: density profile
    :type nprofile: numpy.array
    :param tprofile: temperature profile
    :type tprofile: numpy.array
    :return: pimpden - total impurity radiation density (W/m3)
    :rtype: numpy.array
    """
    # less_than_imp_temp_mask = tprofile values less than impurity temperature. greater_than_imp_temp_mask = tprofile values higher than impurity temperature.
    bins = impurity_radiation_module.impurity_arr_temp_kev[imp_element_index]
    indices = np.digitize(tprofile, bins)
    indices[indices >= bins.shape[0]] = bins.shape[0] - 1
    indices[indices < 0] = 0

    yi = np.log(
        impurity_radiation_module.impurity_arr_lz_wm3[imp_element_index, indices - 1]
    )
    xi = np.log(
        impurity_radiation_module.impurity_arr_temp_kev[imp_element_index, indices - 1]
    )
    c = (
        np.log(
            impurity_radiation_module.impurity_arr_lz_wm3[imp_element_index, indices]
        )
        - yi
    ) / (
        np.log(
            impurity_radiation_module.impurity_arr_temp_kev[imp_element_index, indices]
        )
        - xi
    )
    pimpden = np.exp(yi + c * (np.log(tprofile) - xi))

    pimpden = (
        impurity_radiation_module.impurity_arr_frac[imp_element_index]
        * nprofile
        * nprofile
        * pimpden
    )

    less_than_imp_temp_mask = (
        tprofile
        <= impurity_radiation_module.impurity_arr_temp_kev[imp_element_index, 0]
    )
    pimpden[less_than_imp_temp_mask] = impurity_radiation_module.impurity_arr_lz_wm3[
        imp_element_index, 0
    ]
    # if not impurity_radiation_module.toolow:  # Only print warning once during a run
    #     impurity_radiation_module.toolow = True
    #     error_handling.fdiags[0] = tprofile
    #     error_handling.report_error(35)

    greater_than_imp_temp_mask = (
        tprofile
        >= impurity_radiation_module.impurity_arr_temp_kev[
            imp_element_index,
            impurity_radiation_module.impurity_arr_len_tab[imp_element_index] - 1,
        ]
    )
    #  This is okay because Bremsstrahlung will dominate at higher temp.
    pimpden[greater_than_imp_temp_mask] = impurity_radiation_module.impurity_arr_lz_wm3[
        imp_element_index,
        impurity_radiation_module.impurity_arr_len_tab[imp_element_index] - 1,
    ]

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
        raise ValueError(f"Element {element} is not found in impurity_arr_label") from e


class ImpurityRadiation:
    """This class calculates the impurity radiation losses for given temperature and density profiles.
    The considers the  total impurity radiation from the core (radcore) and total impurity radiation
    (radtot) [MW/(m^3)]. The class is used to sum the impurity radiation loss from each impurity
    element to find the total impurity radiation loss."""

    def __init__(self, plasma_profile):
        """
        :param plasma_profile: Plasma profile class, parameterises the density and temperature profiles.
        :type plasma_profile: Plasma profile class
        """
        self.plasma_profile = plasma_profile
        self.rho = plasma_profile.neprofile.profile_x
        self.rhodx = plasma_profile.neprofile.profile_dx
        self.imp = np.nonzero(impurity_radiation_module.impurity_arr_frac > 1.0e-30)[0]

        self.pimp_profile = np.zeros(self.plasma_profile.profile_size)
        self.radtot_profile = np.zeros(self.plasma_profile.profile_size)
        self.radcore_profile = np.zeros(self.plasma_profile.profile_size)

        self.radtot = 0.0
        self.radcore = 0.0

    def map_imprad_profile(self):
        """Map imprad_profile() over each impurity element index."""
        list(map(self.imprad_profile, self.imp))

    def imprad_profile(self, imp_element_index):
        """
        This routine calculates the impurity radiation losses for given temperature and density profiles.
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

        :param imp_element_index: Index used to access different impurity radiation elements
        :type imp_element_index: Int
        """

        pimp = pimpden(
            imp_element_index,
            self.plasma_profile.neprofile.profile_y,
            self.plasma_profile.teprofile.profile_y,
        )

        self.pimp_profile = np.add(self.pimp_profile, pimp)

    def calculate_radiation_loss_profiles(self):
        """Calculate the Bremsstrahlung (radb), line radiation (radl), total impurity radiation
        from the core (radcore) and total impurity radiation (radtot). Update the stored arrays with
        the values.
        """

        radtot = self.pimp_profile * self.rho
        radcore = self.pimp_profile * (
            self.rho
            * fradcore(
                self.rho,
                impurity_radiation_module.coreradius,
                impurity_radiation_module.coreradiationfraction,
            )
        )

        self.radtot_profile = np.add(self.radtot_profile, radtot)
        self.radcore_profile = np.add(self.radcore_profile, radcore)

    def integrate_radiation_loss_profiles(self):
        """Integrate the radiation loss profiles using the Simpson rule.
        Store the total values for each aspect of impurity radiation loss."""
        self.radtot = 2.0e-6 * integrate.simpson(
            self.radtot_profile, x=self.rho, dx=self.rhodx
        )
        self.radcore = 2.0e-6 * integrate.simpson(
            self.radcore_profile, x=self.rho, dx=self.rhodx
        )

    def calculate_imprad(self):
        """Call the map function to calculate impurity radiation parameters for each
        impurity element. Calculate the radiation loss profiles, and integrate them to
        find the total values for radiation loss.
        """
        self.map_imprad_profile()
        self.calculate_radiation_loss_profiles()
        self.integrate_radiation_loss_profiles()
