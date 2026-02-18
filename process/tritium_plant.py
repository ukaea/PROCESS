from process import constants


class TritiumPlantMeschini:
    """Class to model the tritium plant inventory and flow rates.

        :reference:

            - Samuele Meschini, S. Ferry, Rémi Delaporte-Mathurin, and D. G. Whyte,
            “Modeling and analysis of the tritium fuel cycle for ARC- and STEP-class D-T fusion power plants,”
            Nuclear Fusion, vol. 63, no. 12, pp. 126005-126005, Sep. 2023,
            doi: https://doi.org/10.1088/1741-4326/acf3fc.
    ‌

    """

    def __init__(self) -> None:
        self.outfile: int = constants.NOUT

    # 1 = Blanket
    # 2 = Tritium Extraction System
    # 3 = First Wall
    # 4 = Divertor
    # 5 = Heat Exchanger
    # 6 = Detritiation System
    # 7 = Vacuum Pump
    # 8 = Fuel cleanup
    # 9 = Isotope Separation System
    # 10 Storage and management
    # 11 = Fuelling system
    # 12 = Tritium sepration membrane


RATE_T_DECAY = 1.78e-9  # Tritium decay rate (1/s)

TBR = 0.5
# Tritium burn efficiency in the plasma
TBE = 0.05

# Direct internal recycling fraction
f_dir = 0.1

# Tritium residence time in the ith component
tau_1 = 5.0
tau_2 = 10.0
tau_3 = 15.0
tau_4 = 20.0
tau_5 = 25.0
tau_6 = 30.0
tau_7 = 35.0
tau_8 = 40.0
tau_9 = 45.0
tau_10 = 50.0
tau_11 = 55.0
tau_12 = 60.0

n_t_burn = 1e-3  # Tritium burn rate in the plasma (kg/s)

i_startup = 5.0

ETA_2 = 0.5

# Flow rate fractions between components
f_5_1 = 0.33
f_5_3 = 0.33
f_5_6 = 1e-4
f_9_6 = 0.1
f_p_3 = 1e-4
f_p_4 = 1e-4
f_5_4 = 0.33

# The non-radioactive loss fraction has been assumed asb in this work.
epsilon_1 = 1e-4
epsilon_2 = 1e-4
epsilon_3 = 0.0
epsilon_4 = 0.0
epsilon_5 = 1e-4
epsilon_6 = 1e-4
epsilon_7 = 1e-4
epsilon_8 = 1e-4
epsilon_9 = 1e-4
epsilon_10 = 0.0
epsilon_11 = 1e-4
epsilon_12 = 1e-4

# Tritium inventory in the ith component (kg)

m_tritium_component_1 = 0.0
m_tritium_component_2 = 0.0
m_tritium_component_3 = 0.0
m_tritium_component_4 = 0.0
m_tritium_component_5 = 0.0
m_tritium_component_6 = 0.0
m_tritium_component_7 = 0.0
m_tritium_component_8 = 0.0
m_tritium_component_9 = 0.0
m_tritium_component_10 = i_startup
m_tritium_component_11 = 0.0
m_tritium_component_12 = 0.0


def calculate_rates(
    m_tritium_component_1,
    m_tritium_component_2,
    m_tritium_component_3,
    m_tritium_component_4,
    m_tritium_component_5,
    m_tritium_component_6,
    m_tritium_component_7,
    m_tritium_component_8,
    m_tritium_component_9,
    m_tritium_component_10,
    m_tritium_component_12,
):
    """Calculate rates of change for all components."""
    d_m_tritium_component_1_dt = (
        TBR * n_t_burn
        + (m_tritium_component_3 / tau_3)
        + (m_tritium_component_4 / tau_4)
        + f_5_1 * (m_tritium_component_5 / tau_5)
        - (m_tritium_component_1 * (((1 + epsilon_1) / tau_1) + RATE_T_DECAY))
    )

    d_m_tritium_component_2_dt = (m_tritium_component_1 / tau_1) - (
        m_tritium_component_2 * (((1 + epsilon_2) / tau_2) + RATE_T_DECAY)
    )

    d_m_tritium_component_3_dt = (
        f_p_3 * (n_t_burn / TBE)
        + f_5_3 * (m_tritium_component_5 / tau_5)
        - (m_tritium_component_3 * (((1 + epsilon_3) / tau_3) + RATE_T_DECAY))
    )

    d_m_tritium_component_4_dt = (
        f_p_4 * (n_t_burn / TBE)
        + f_5_4 * (m_tritium_component_5 / tau_5)
        - (m_tritium_component_4 * (((1 + epsilon_4) / tau_4) + RATE_T_DECAY))
    )

    d_m_tritium_component_5_dt = (1 - ETA_2) * (m_tritium_component_2 / tau_2) - (
        m_tritium_component_5 * (((1 + epsilon_5) / tau_5) + RATE_T_DECAY)
    )

    d_m_tritium_component_6_dt = (
        f_5_6 * (m_tritium_component_5 / tau_5)
        + f_9_6 * (m_tritium_component_9 / tau_9)
        - (m_tritium_component_6 * (((1 + epsilon_6) / tau_6) + RATE_T_DECAY))
    )

    d_m_tritium_component_7_dt = (1 - TBE - f_p_3 - f_p_4) * (n_t_burn / TBE) - (
        m_tritium_component_7 * (((1 + epsilon_7) / tau_7) + RATE_T_DECAY)
    )

    d_m_tritium_component_8_dt = (1 - f_dir) * (m_tritium_component_7 / tau_7) - (
        m_tritium_component_8 * (((1 + epsilon_8) / tau_8) + RATE_T_DECAY)
    )

    d_m_tritium_component_9_dt = (
        (m_tritium_component_6 / tau_6)
        + (m_tritium_component_8 / tau_8)
        - (m_tritium_component_9 * (((1 + epsilon_9) / tau_9) + RATE_T_DECAY))
    )

    d_m_tritium_component_10_dt = (
        (1 - f_9_6) * (m_tritium_component_9 / tau_9)
        + f_dir * (m_tritium_component_7 / tau_7)
        + (m_tritium_component_12 / tau_12)
        - (n_t_burn / TBE)
        - (RATE_T_DECAY * m_tritium_component_10)
    )

    d_m_tritium_component_12_dt = ETA_2 * (m_tritium_component_2 / tau_2) - (
        m_tritium_component_12 * (((1 + epsilon_12) / tau_12) + RATE_T_DECAY)
    )

    return (
        d_m_tritium_component_1_dt,
        d_m_tritium_component_2_dt,
        d_m_tritium_component_3_dt,
        d_m_tritium_component_4_dt,
        d_m_tritium_component_5_dt,
        d_m_tritium_component_6_dt,
        d_m_tritium_component_7_dt,
        d_m_tritium_component_8_dt,
        d_m_tritium_component_9_dt,
        d_m_tritium_component_10_dt,
        d_m_tritium_component_12_dt,
    )


def simulate_tritium_inventory(t_end, dt):
    """
    Simulate tritium inventory evolution over time using Euler method.

    Parameters:
    t_end: End time for simulation (s)
    dt: Time step (s)

    Returns:
    times: Array of time points
    inventories: Array of inventories for each component at each time point
    """
    n_steps = int(t_end / dt) + 1
    times = [i * dt for i in range(n_steps)]

    # Initialize inventory arrays
    inventories = {
        "i_1": [m_tritium_component_1],
        "i_2": [m_tritium_component_2],
        "i_3": [m_tritium_component_3],
        "i_4": [m_tritium_component_4],
        "i_5": [m_tritium_component_5],
        "i_6": [m_tritium_component_6],
        "i_7": [m_tritium_component_7],
        "i_8": [m_tritium_component_8],
        "i_9": [m_tritium_component_9],
        "i_10": [m_tritium_component_10],
        "i_11": [m_tritium_component_11],
        "i_12": [m_tritium_component_12],
    }

    # Current values
    curr = [
        m_tritium_component_1,
        m_tritium_component_2,
        m_tritium_component_3,
        m_tritium_component_4,
        m_tritium_component_5,
        m_tritium_component_6,
        m_tritium_component_7,
        m_tritium_component_8,
        m_tritium_component_9,
        m_tritium_component_10,
        m_tritium_component_11,
        m_tritium_component_12,
    ]

    for step in range(1, n_steps):
        # Calculate rates
        rates = calculate_rates(
            curr[0],
            curr[1],
            curr[2],
            curr[3],
            curr[4],
            curr[5],
            curr[6],
            curr[7],
            curr[8],
            curr[9],
            curr[11],
        )

        # Update inventories using Euler method
        curr[0] += rates[0] * dt
        curr[1] += rates[1] * dt
        curr[2] += rates[2] * dt
        curr[3] += rates[3] * dt
        curr[4] += rates[4] * dt
        curr[5] += rates[5] * dt
        curr[6] += rates[6] * dt
        curr[7] += rates[7] * dt
        curr[8] += rates[8] * dt
        curr[9] += rates[9] * dt
        curr[11] += rates[10] * dt
        # i_11 remains 0

        # Store values
        for idx, key in enumerate([
            "i_1",
            "i_2",
            "i_3",
            "i_4",
            "i_5",
            "i_6",
            "i_7",
            "i_8",
            "i_9",
            "i_10",
            "i_11",
            "i_12",
        ]):
            inventories[key].append(curr[idx])

    return times, inventories
