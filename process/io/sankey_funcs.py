"""
Library of Sankey plotting routine

Author: Hanni Lux (Hanni.Lux@ukaea.uk)
        Matti Coleman (Matti.Coleman@ukaea.uk)

Updated 13/09/2019: Adam Brown (adam.brown@ukaea.uk)
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.sankey import Sankey
from numpy import sqrt

from process.io.mfile import MFile


def plot_full_sankey(
    mfilename="MFILE.DAT",
):  # Plots the power flow from PROCESS as a Sankey Diagram
    # ------------------------------- Pulling values from the MFILE -------------------------------

    m_file = MFile(mfilename)

    # Used in [PLASMA]
    fusion_power = m_file.data["fusion_power"].get_scan(-1)  # Fusion power (MW)
    p_hcd_injected_total_mw = m_file.data["p_hcd_injected_total_mw"].get_scan(
        -1
    )  # Total auxiliary injected power (MW)
    p_plasma_ohmic_mw = m_file.data["p_plasma_ohmic_mw"].get_scan(
        -1
    )  # Ohmic heating power (MW)
    totalplasma = (
        fusion_power + p_hcd_injected_total_mw + p_plasma_ohmic_mw
    )  # Total Power in plasma (MW)
    neutron_power_total = m_file.data["neutron_power_total"].get_scan(
        -1
    )  # Neutron fusion power (MW)
    non_alpha_charged_power = m_file.data["non_alpha_charged_power"].get_scan(
        -1
    )  # Non-alpha charged particle power (MW)
    pcharohmmw = (
        non_alpha_charged_power + p_plasma_ohmic_mw
    )  # The ohmic and charged particle power (MW)
    alpha_power_total = m_file.data["alpha_power_total"].get_scan(
        -1
    )  # Alpha power (MW)
    palpinjmw = (
        alpha_power_total + p_hcd_injected_total_mw
    )  # Alpha particle and HC&D power (MW)

    # Used in [NEUTRONICS]
    emultmw = m_file.data["emultmw"].get_scan(
        -1
    )  # Energy multiplication in blanket (MW)
    p_blkt_nuclear_heat_total_mw = m_file.data["p_blkt_nuclear_heat_total_mw"].get_scan(
        -1
    )  # Total Nuclear heating in the blanket (MW)
    pnucemblkt = (
        p_blkt_nuclear_heat_total_mw - emultmw
    )  # External nuclear heating in blanket (MW)
    p_div_nuclear_heat_total_mw = m_file.data["p_div_nuclear_heat_total_mw"].get_scan(
        -1
    )  # Nuclear heating in the divertor (MW)
    p_fw_nuclear_heat_total_mw = m_file.data["p_fw_nuclear_heat_total_mw"].get_scan(
        -1
    )  # Nuclear heating in the first wall (MW)
    pnucshld = m_file.data["pnucshld"].get_scan(
        -1
    )  # Nuclear heating in the shield (MW)
    ptfnuc = m_file.data["ptfnuc"].get_scan(-1)  # Nuclear heating in the TF coil (MW)

    # Used in [CHARGEP]
    p_plasma_separatrix_mw = m_file.data["p_plasma_separatrix_mw"].get_scan(
        -1
    )  # Charged particle power deposited on divertor (MW)
    f_alpha_plasma = m_file.data["f_alpha_plasma"].get_scan(
        -1
    )  # Fraction of alpha power deposited in plasma
    p_fw_alpha_mw = alpha_power_total * (
        1 - f_alpha_plasma
    )  # Alpha particles hitting first wall (MW)
    p_plasma_rad_mw = m_file.data["p_plasma_rad_mw"].get_scan(
        -1
    )  # Total radiation Power (MW)

    # Used in [RADIATION]
    p_div_rad_total_mw = p_plasma_rad_mw * m_file.data["f_ster_div_single"].get_scan(
        -1
    )  # Radiation deposited on the divertor (MW)
    p_fw_hcd_rad_total_mw = p_plasma_rad_mw * m_file.data["f_a_fw_hcd"].get_scan(
        -1
    )  # Radiation deposited on HCD (MW)
    p_fw_rad_total_mw = (
        p_plasma_rad_mw - p_div_rad_total_mw - p_fw_hcd_rad_total_mw
    )  # Radiation deposited in the FW (MW)

    # Used in [DIVERTOR]
    htpmw_div = m_file.data["htpmw_div"].get_scan(-1)  # Divertor coolant pumping power
    pthermdiv = m_file.data["pthermdiv"].get_scan(
        -1
    )  # Total power extracted from divertor (MW)

    # Used in [FIRST_WALL]
    pthermfw_blkt = m_file.data["pthermfw_blkt"].get_scan(
        -1
    )  # Power extracted blanket & FW (MW)
    htpmw_fw_blkt = m_file.data["htpmw_fw_blkt"].get_scan(
        -1
    )  # Pump Power in FW and blanket (MW)
    htpmwblkt = htpmw_fw_blkt / 2  # Pump power in blanket (MW)
    htpmwfw = htpmw_fw_blkt / 2  # Pump power in FW (MW)
    pthermfw = (
        pthermfw_blkt - htpmwblkt - p_blkt_nuclear_heat_total_mw
    )  # Power extracted 1st wall (MW)
    # porbitloss = m_file.data['porbitloss'].get_scan(-1) # Charged P. on FW before thermalising
    # p_beam_shine_through_mw = m_file.data['p_beam_shine_through_mw'].get_scan(-1) # Injection shine-through to 1st wall

    # Initialising x and y variables for adjusting 'Plasma Heating' branch tip location
    y_adj_1 = 0
    y_adj_2 = 0

    # Loop 1 to get 'Plasma Heating' branch tip coords; loop 2 to match 'PLASMA' branch
    for _ in range(2):
        # The visual settings of the Sankey Plot
        plt.rcParams.update({"font.size": 9})
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, xticks=[], yticks=[], frameon=False)
        sankey = Sankey(
            ax=ax, unit="MW", margin=0.5, format="%1.0f", scale=1.0 / (totalplasma)
        )

        # --------------------------------------- PLASMA - 0 --------------------------------------

        # Fusion, Injected, Ohmic, -Charged P.-Ohmic, -Alphas-Injected, -Neutrons
        plasma = [
            fusion_power,
            p_hcd_injected_total_mw,
            p_plasma_ohmic_mw,
            -pcharohmmw,
            -palpinjmw,
            -neutron_power_total,
        ]
        sankey.add(
            flows=plasma,
            # [left(in), down(in), down(in), up(out), up(out), right(out)]
            orientations=[0, -1, -1, 1, 1, 0],
            trunklength=0.5,
            pathlengths=[0.5, 0.25, 0.25, 0.75, 0.25 + 0.5 * y_adj_1, 0.0],
            # labels=["Fusion","H&CD", "Ohmic", "Charged P.", "Alphas", "Neutrons"])
            labels=[None, None, None, None, None, None],
        )

        # Check to see if the fusion components balance
        if _ == 0 and sqrt(sum(plasma) ** 2) > 0.1:
            print("FUSION power balance =", sum(plasma), "\n")
            exit()

        if _ == 1:
            print(sankey.finish()[0])
        if _ == 1:
            print(sankey.finish()[0].patch)
        if _ == 1:
            print(type(sankey.finish()[0].patch))

        # ------------------------------------- NEUTRONICS - 1 ------------------------------------

        # Neutrons, -Divertor, -1st wall, -Shield, -TF coils, -Blanket+Energy Mult.
        neutrons = [
            neutron_power_total,
            -p_div_nuclear_heat_total_mw,
            -p_fw_nuclear_heat_total_mw,
            -pnucshld,
            -ptfnuc,
            -pnucemblkt,
        ]
        sankey.add(
            flows=neutrons,
            # left(in), up(out), up(out), up(out), up(out), right(out)
            orientations=[0, 1, 1, 1, 1, 0],
            trunklength=0.5,
            pathlengths=[0.3, 0.25, 0.25, 0.25, 0.25, 0.15],
            prior=0,  # PLASMA
            connect=(5, 0),  # Neutrons
            # labels=["Neutrons", "Divertor", "1st Wall", "Shield", "TF coils", "Blanket"])
            labels=[None, None, None, None, None, None],
        )

        # Checking to see if the neutronics components balance
        if _ == 0 and sqrt(sum(neutrons) ** 2) > 0.1:
            print("NEUTRONS power balance =", sum(neutrons), "\n")
            exit()

        # Check to see if connections balance
        if _ == 0:
            check = sankey.finish()
            diff1_1 = check[0].flows[5] + check[1].flows[0]
            plt.close()
            if diff1_1 > 0.1:
                print("Neutrons [0][5] and [1][0] difference =", diff1_1)
                exit()

        # --------------------------------- CHARGED PARTICLES - 2 ---------------------------------

        # Charge P.+Ohmic, Alpha+Injected, -Divertor, -1st Wall, -Photons
        chargedp = [
            pcharohmmw,
            palpinjmw,
            -p_plasma_separatrix_mw,
            -p_fw_alpha_mw,
            -p_plasma_rad_mw,
        ]
        sankey.add(
            flows=chargedp,
            # down(in), down(in), up(out), up(out), right(out)
            orientations=[-1, -1, 1, -1, 0],
            trunklength=0.5,
            pathlengths=[0.75, 0.25 + 0.5 * y_adj_1, 0.25, 0.25, 0.25],
            prior=0,  # PLASMA
            connect=(3, 0),  # Charged P.+Ohmic
            # labels=["Charged P.", "Alphas", "Divertor", "1st Wall", "Photons"])
            labels=[None, None, None, None, None],
        )

        if _ == 0 and sqrt(sum(chargedp) ** 2) > 0.1:
            print("CHARGEDP power balance =", sum(chargedp))
            exit()

        # Check to see if connections balance
        if _ == 0:
            check = sankey.finish()
            diff2_1 = check[0].flows[3] + check[2].flows[0]
            diff2_2 = check[0].flows[4] + check[2].flows[1]
            plt.close()
            if diff2_1 > 0.1:
                print("Charged P.+Ohmic [0][3] and [2][0] difference =", diff2_1)
                exit()
            if diff2_2 > 0.1:
                print("Alphas+Injected [0][4] and [2][1] difference =", diff2_2)
                exit()

        # ------------------------------------- RADIATION - 3 -------------------------------------

        # Photons, -1st Wall, -Divertor, -H&CD
        radiation = [
            p_plasma_rad_mw,
            -p_fw_rad_total_mw,
            -p_div_rad_total_mw,
            -p_fw_hcd_rad_total_mw,
        ]
        sankey.add(
            flows=radiation,
            # right(in), up(out), up(out), up(out)
            orientations=[
                0,
                -1,
                1,
                1,
            ],
            trunklength=0.5,
            pathlengths=[0.25, 0.25, 0.25, 0.25],
            prior=2,  # CHARGED PARTICLES
            connect=(4, 0),  # Charged P.
            # labels=["Photons", "1st Wall", "Divertor", "H&CD"])
            labels=[None, None, None, None],
        )

        if _ == 0 and sqrt(sum(radiation) ** 2) > 0.1:
            print("RADIATION power balance =", sum(radiation))
            exit()

        if _ == 0:
            check = sankey.finish()
            diff3_1 = check[2].flows[4] + check[3].flows[0]
            plt.close()
            if diff3_1 > 0.1:
                print("Photons [2][4] and [3][0] difference =", diff3_1)
                exit()

        # -------------------------------------- DIVERTOR - 4 -------------------------------------

        # Charged P., Neutrons, Photons, Coolant Pumping, Total Divertor
        divertor = [
            p_plasma_separatrix_mw,
            p_div_nuclear_heat_total_mw,
            p_div_rad_total_mw,
            htpmw_div,
            -pthermdiv,
        ]
        sankey.add(
            flows=divertor,
            # down(in), up(in), down(in), up(in), right(out)
            orientations=[-1, -1, -1, -1, 0],
            trunklength=0.5,
            pathlengths=[0.25, 0.25, 0.25, 0.25 - 0.5 * y_adj_2, 0.25],
            prior=2,  # CHARGED PARTICLES
            connect=(2, 0),  # Charged P. --> None
            # labels=["Charged P.", "Neutrons", "Photons", "Coolant Pumping", "Divertor Power"])
            labels=[None, None, None, None, None],
        )

        if _ == 0 and sqrt(sum(divertor) ** 2) > 0.1:
            print("DIVERTOR power balance =", sum(divertor))
            exit()

        if _ == 0:
            check = sankey.finish()
            diff4_1 = check[1].flows[1] + check[4].flows[0]
            diff4_2 = check[2].flows[3] + check[4].flows[3]
            plt.close()
            if diff4_1 > 0.1:
                print("Neutrons [1][1] and [4][0] difference =", diff4_1)
                exit()
            if diff4_2 > 0.1:
                print("Charged P. [2][3] and [4][3] difference =", diff4_2)
                exit()

        # ---------------------------------------- 1ST WALL - 5 ---------------------------------------

        # Alphas, Neutrons, Photons, Coolant Pumping, Total 1st Wall
        first_wall = [
            p_fw_alpha_mw,
            p_fw_nuclear_heat_total_mw,
            p_fw_rad_total_mw,
            htpmwfw,
            -pthermfw,
        ]
        sankey.add(
            flows=first_wall,
            orientations=[0, -1, 1, -1, 0],
            trunklength=0.5,
            pathlengths=[0.25, 0.25, 0.25, 0.25, 0.25],
            prior=1,
            connect=(2, 1),
            # labels=["Alphas", "Neutrons", "Radiation", "Coolant Pumping", "FW Power"])
            labels=[None, None, None, None, None],
        )

        if _ == 0 and sqrt(sum(first_wall) ** 2) > 0.1:
            print("FIRST_WALL power balance =", sum(first_wall))
            exit()

        """# -------------------------------------- BLANKET - 6 --------------------------------------

        # Blanket - Energy mult., Energy Mult., pumping power, Blanket
        BLANKET = [pnucemblkt, emultmw, htpmwblkt, -pthermblkt]
        sankey.add(flows=BLANKET,
                   # left(in), down(in), down(in), right(out)
                   orientations=[0, -1, -1, 0],
                   trunklength=0.5,
                   pathlengths=[0.25, 0.25, 0.25, 0.25],
                   #prior=1, # NEUTRONICS
                   #connect=(1, 0), # Blanket --> None
                   labels=[None, "Energy Mult.", "Coolant Pumping", "Blanket"])

        # Checking to see if the blanket components balance
        if _ == 0:
            if sqrt(sum(BLANKET)**2) > 0.1:
                print("BLANKET power balance =", sum(BLANKET), "\n")
                exit()

        # Check to see if connections balance
        if _ == 0:
            check = sankey.finish()
            diff = check[1].flows[1]+check[3].flows[0]
            if diff > 0.1:
                print("The difference between [1][1] and [3][0] =", diff)
                exit()"""

        """# --------------------------------------- SHIELD - 7 --------------------------------------

        # Neutrons, Coolant pumping, Total power
        SHIELD = [pnucshld, htpmw_shld, -pthermshld]
        sankey.add(flows=SHIELD,
                   orientations=[-1, -1, 1],
                   trunklength=0.5,
                   pathlengths=[0.25, 0.25 ,0.25],
                   #prior=2,
                   #connect=(5, 0),
                   labels=["Neutrons", "Coolant Pumping", "Shield Power"])

        if _ == 0:
            if sqrt(sum(SHIELD)**2) > 0.1:
                print("SHIELD power balance =", sum(SHIELD))
                exit()"""

        """# ------------------------------------ PRIMARY HEAT - 7 -----------------------------------

        # 1st wall, Blanket, Shield, Divertor, Total thermal power
        HEAT = [pthermfw, pthermblkt, pthermshld, pthermdiv, -pthermmw]
        sankey.add(flows=HEAT,
                   orientations=[1, 0, -1, 1, 0],
                   trunklength=0.5,
                   pathlengths=[0.25, 0.25 ,0.25, 0.25, 0.25],
                   #prior=2,
                   #connect=(5, 0),
                   labels=["1st Wall", "Blanket", "Shield", "Divertor", "Total Power"])

        if _ == 0:
            if sqrt(sum(HEAT)**2) > 0.1:
                print("PRIMARY power balance =", sum(HEAT))
                exit()"""

        """# ------------------------------- ELECTRICITY CONVERSION - 8 ------------------------------

        # Total thermal, Elctricty conversion loss, Gross Electricity
        GROSS = [pthermmw, -pelectloss, -pgrossmw]
        sankey.add(flows=GROSS,
                   orientations=[0, -1, 0],
                   trunklength=0.5,
                   pathlengths=[0.25, 0.25 ,0.25],
                   #prior=2,
                   #connect=(5, 0),
                   labels=["Thermal Power", "Conversion loss", "Gross Electricity"])

        if _ == 0:
            if sqrt(sum(GROSS)**2) > 0.1:
                print("GROSS power balance =", sum(GROSS))
                exit()"""

        # ------------------------------ RECIRCULATED ELECTRICITY - 9 -----------------------------

        """# ---------------------------------------- HCD - 11 ----------------------------------------

        # HCD loss + injected, -injected, -HCD loss
        HCD = [pinjht+p_hcd_injected_total_mw, -p_hcd_injected_total_mw, -pinjht]
        assert(sum(HCD)**2 < 0.5)
        sankey.add(flows=HCD,
                   # [down(in), up(out), down(out)]
                   orientations=[-1, 1, -1],
                   #prior=0, # PLASMA
                   #connect=(1, 1), # H&CD --> None
                   trunklength=0.5,
                   pathlengths=[0.25, 0.25, 0.25],
                   labels=['H&CD power', None, 'H&CD loss'])"""

        fig.tight_layout()

        if _ == 0:
            plt.close()

        # Matching PLASMA and CHARGED PARTICLES 'Alphas' branches
        # x_adj_1, y_adj_1 = diagrams[2].tips[1] - diagrams[0].tips[4]
        # Matching CHARGED PARTICLES and DIVERTOR 'Charged P.' branches
        # x_adj_2, y_adj_2 = diagrams[4].tips[3] - diagrams[2].tips[3]
        # x_adj_3, y_adj_3 = diagrams[3].tips[3] - diagrams[4].tips[0]

    # --------------------------------------- Label Positioning ---------------------------------------

    # Munipulating the positioning of the branch labels
    # -ve to left and down; +ve to right and up
    # pos[0] = x-axis; pos[1] = y-axis

    """for d in diagrams:
        y = 0
        for t in d.texts:
            pos = tuple(np.ndarray.tolist(d.tips[y]))
            t.set_position(pos)
            if t == diagrams[0].texts[0]: # Fusion Power
                t.set_horizontalalignment('right')
                t.set_position((pos[0]-0.2,pos[1]))
            if t == diagrams[0].texts[1]: # H&CD
                t.set_horizontalalignment('right')
                t.set_position((pos[0]-0.5*(p_hcd_injected_total_mw/totalplasma)-0.05,pos[1]))
            if t == diagrams[0].texts[2]: # Ohmic
                t.set_horizontalalignment('left')
                t.set_position((pos[0]+0.5*(p_plasma_ohmic_mw/totalplasma)+0.05,pos[1]))
            if t == diagrams[0].texts[3]: # Neutrons
                t.set_horizontalalignment('right')
                t.set_position((pos[0]-0.2,pos[1]))
            if t == diagrams[0].texts[4]: # Charged Particles
                t.set_horizontalalignment('right')
                t.set_position((pos[0]-0.5*(non_alpha_charged_power/totalplasma)-0.05,pos[1]))
            if t == diagrams[0].texts[5]: # Alphas
                t.set_horizontalalignment('left')
                t.set_position((pos[0]+0.5*(alpha_power_total/totalplasma)+0.05,pos[1]-0.1))
            if t == diagrams[1].texts[0]: # H&CD power
                t.set_horizontalalignment('right')
                t.set_position((pos[0]-0.5*((pinjht+p_hcd_injected_total_mw)/totalplasma)-0.05,pos[1]))
            if t == diagrams[1].texts[2]: # H&CD losses
                t.set_horizontalalignment('left')
                t.set_position((pos[0]+(pinjht/totalplasma)+0.05,pos[1]))
            if t == diagrams[2].texts[1]: # Energy Multiplication
                t.set_horizontalalignment('center')
                t.set_position((pos[0],pos[1]-0.2))
            if t == diagrams[2].texts[2]: # Blanket
                t.set_horizontalalignment('right')
                t.set_position((pos[0]-0.2,pos[1]))
            if t == diagrams[2].texts[3]: # Divertor
                t.set_horizontalalignment('right')
                t.set_position((pos[0]-0.5*(p_div_nuclear_heat_total_mw/totalplasma)-0.1,pos[1]))
            if t == diagrams[3].texts[2]: # Rad.FW
                t.set_horizontalalignment('right')
                t.set_position((pos[0],pos[1]+0.5*(p_fw_rad_total_mw/totalplasma)+0.15))
            if t == diagrams[3].texts[3]: # Charged P.
                t.set_horizontalalignment('left')
                t.set_position((pos[0]+0.5*((p_plasma_separatrix_mw+p_fw_alpha_mw)/totalplasma)+0.1,pos[1]+0.05))
            if t == diagrams[3].texts[4]: # Rad. Div.
                t.set_horizontalalignment('right')
                t.set_position((pos[0]-0.5*(p_div_rad_total_mw/totalplasma)-0.1,pos[1]))
            y += 1"""


####################################################################################################


def plot_sankey(mfilename="MFILE.DAT"):  # Plot simplified power flow Sankey Diagram
    # ------------------------------- Pulling values from the MFILE -------------------------------

    m_file = MFile(mfilename)

    # Used in [PLASMA]
    fusion_power = m_file.data["fusion_power"].get_scan(-1)  # Fusion Power (MW)
    p_hcd_injected_total_mw = m_file.data["p_hcd_injected_total_mw"].get_scan(
        -1
    )  # Total auxiliary injected Power (MW)
    p_plasma_ohmic_mw = m_file.data["p_plasma_ohmic_mw"].get_scan(
        -1
    )  # Ohmic heating Power (MW)
    totalplasma = (
        fusion_power + p_hcd_injected_total_mw + p_plasma_ohmic_mw
    )  # Total Power in plasma (MW)

    # Used in [DEPOSITION]
    p_plasma_rad_mw = m_file.data["p_plasma_rad_mw"].get_scan(
        -1
    )  # Total radiation Power (MW)
    f_ster_div_single = m_file.data["f_ster_div_single"].get_scan(
        -1
    )  # Area fraction taken up by divertor
    fdiv_2 = m_file.data["2*f_ster_div_single"].get_scan(
        -1
    )  # Area fraction taken up by double null divertor
    if (
        fdiv_2 > 0
    ):  # Takes into account old MFILE representation of double null divertor
        f_ster_div_single = fdiv_2
    p_div_rad_total_mw = (
        p_plasma_rad_mw * f_ster_div_single
    )  # Radiation deposited on the divertor (MW)
    f_a_fw_hcd = m_file.data["f_a_fw_hcd"].get_scan(
        -1
    )  # Area fraction covered by HCD and diagnostics
    p_fw_hcd_rad_total_mw = (
        p_plasma_rad_mw * f_a_fw_hcd
    )  # Radiation deposited on HCD and diagnostics (MW)
    p_fw_rad_total_mw = (
        p_plasma_rad_mw - p_div_rad_total_mw - p_fw_hcd_rad_total_mw
    )  # Radiation deposited in the blanket (MW)
    p_plasma_separatrix_mw = m_file.data["p_plasma_separatrix_mw"].get_scan(
        -1
    )  # power to conducted to the divertor region (MW)
    p_div_nuclear_heat_total_mw = m_file.data["p_div_nuclear_heat_total_mw"].get_scan(
        -1
    )  # nuclear heating in the divertor (MW)
    p_fw_nuclear_heat_total_mw = m_file.data["p_fw_nuclear_heat_total_mw"].get_scan(
        -1
    )  # nuclear heating in the first wall (MW)
    p_blkt_nuclear_heat_total_mw = m_file.data["p_blkt_nuclear_heat_total_mw"].get_scan(
        -1
    )  # nuclear heating in the blanket (MW)
    pnucshld = m_file.data["pnucshld"].get_scan(
        -1
    )  # nuclear heating in the shield (MW)
    pnuc_cp_sh = m_file.data["pnuc_cp_sh"].get_scan(
        -1
    )  # nuclear heating in the CP shield (MW)
    emultmw = m_file.data["emultmw"].get_scan(-1)  # Blanket energy multiplication (MW)
    alpha_power_total = m_file.data["alpha_power_total"].get_scan(
        -1
    )  # Alpha power (MW)
    f_alpha_plasma = m_file.data["f_alpha_plasma"].get_scan(
        -1
    )  # Fraction of alpha power deposited in plasma
    p_fw_alpha_mw = alpha_power_total * (
        1 - f_alpha_plasma
    )  # Alpha power hitting 1st wall (MW)
    itart = m_file.data["itart"].get_scan(
        -1
    )  # switch for spherical tokamak (ST) models

    # Power deposited on divertor (MW)
    totaldivetc = (
        p_plasma_separatrix_mw + p_div_nuclear_heat_total_mw + p_div_rad_total_mw
    )
    # Power deposited on Blanket (MW)
    totalblktetc = (
        p_fw_nuclear_heat_total_mw
        + p_blkt_nuclear_heat_total_mw
        + pnucshld
        + p_fw_rad_total_mw
        + p_fw_alpha_mw
        - emultmw
    )

    if itart == 0:
        # Power deposited in CP (MW) (None here)
        totalcpetc = 0.0
    elif itart == 1:
        # Power deposited in CP (MW)
        totalcpetc = pnuc_cp_sh

    # Used in [BLANKETSETC]
    pthermfw_blkt = m_file.data["pthermfw_blkt"].get_scan(
        -1
    )  # Heat for electricity (MW)
    htpmw_fw_blkt = m_file.data["htpmw_fw_blkt"].get_scan(
        -1
    )  # 1st wall & blanket pumping (MW)
    pthermmw_p = pthermfw_blkt - htpmw_fw_blkt  # Heat - pumping power (MW)

    # Used in [PRIMARY]
    pgrossmw = m_file.data["pgrossmw"].get_scan(-1)  # gross electric power (MW)

    # Used in [NET]
    pnetelmw = m_file.data["pnetelmw"].get_scan(-1)  # net electric power (MW)
    precircmw = pgrossmw - pnetelmw  # Recirculating power (MW)

    # Used in [RECIRC]
    crypmw = m_file.data["crypmw"].get_scan(-1)  # cryogenic plant power (MW)
    fachtmw = m_file.data["fachtmw"].get_scan(-1)  # facility heat removal (MW)
    tfacpd = m_file.data["tfacpd"].get_scan(
        -1
    )  # total steady state TF coil AC power demand (MW)
    trithtmw = m_file.data["trithtmw"].get_scan(
        -1
    )  # power required for tritium processing (MW)
    vachtmw = m_file.data["vachtmw"].get_scan(-1)  # vacuum pump power (MW)
    pfwpmw = m_file.data["pfwpmw"].get_scan(
        -1
    )  # Total mean wall plug power for PFC & CS (MW)
    ppumpmw = (
        m_file.data["ppump"].get_scan(-1) / 1e6
    )  # Set pumping power to MW by dividing by 1e6

    # Energy required for rest of power plant (MW)
    pcoresystems = crypmw + fachtmw + tfacpd + trithtmw + vachtmw + pfwpmw + ppumpmw
    pinjwp = m_file.data["pinjwp"].get_scan(-1)  # injector wall plug power (MW)
    htpmw = m_file.data["htpmw"].get_scan(
        -1
    )  # heat transport system electrical pump power (MW)

    # Initialising x and y variables for adjusting 'Plasma Heating' branch tip location
    x_adj, y_adj = 0, 0

    # Loop 1 to get 'Plasma Heating' branch tip coords; loop 2 to match 'PLASMA' branch
    for _ in range(2):
        # ------------------------------------ Visual Settings ------------------------------------

        plt.rcParams.update({"font.size": 9})  # Setting font size to 9
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, xticks=[], yticks=[], frameon=False)
        sankey = Sankey(
            ax=ax, unit="MW", margin=0.0, format="%1.0f", scale=1.0 / (totalplasma)
        )

        # --------------------------------------- PLASMA - 0 --------------------------------------

        # Fusion power, Injected power + ohmic power, - total plasma power
        plasma = [
            fusion_power,
            p_hcd_injected_total_mw + p_plasma_ohmic_mw,
            -totalplasma,
        ]
        sankey.add(
            flows=plasma,
            orientations=[0, -1, 0],  # [right(in), down(in), right(out)]
            pathlengths=[
                0.5,
                0.8 + 0.5 * y_adj,
                -0.1 + 0.5 * x_adj,
            ],  # 'Plasma Heating' adjust
            labels=["Fusion Power", None, "Plasma"],
        )

        # --------------------------------- ENERGY DEPOSITION - 1 ---------------------------------

        # Plasma power, - divertor deposited power, - blanket deposited power
        deposition = [totalplasma, -totalblktetc - totaldivetc - totalcpetc]
        # Check if difference >2 between plasma and divertor + blanket
        if _ == 1 and sqrt(sum(deposition) ** 2) > 2:
            print(
                "\ncomponents power balance difference =",
                totalplasma - totaldivetc - totalblktetc - totalcpetc,
            )
        sankey.add(
            flows=deposition,
            orientations=[0, 0],  # [right(in), up(in), right(out)]
            prior=0,  # PLASMA
            connect=(2, 0),  # Plasma --> None
            pathlengths=[0.2, 0.2 + 0.5 * x_adj],  # 'Plasma Heating' adjust
            labels=[None, "Blanket/etc."],
        )

        # -------------------------------------- BLANKET - 2 --------------------------------------

        # Blanket deposited power, blanket energy multiplication, - primary heat
        blanketsetc = [
            totalblktetc + totaldivetc + totalcpetc,
            emultmw,
            -pthermmw_p - totaldivetc - totalcpetc - pnucshld,
        ]
        # Check if difference >2 between primary heat and blanket + blanket multiplication
        if _ == 1 and sqrt(sum(blanketsetc) ** 2) > 2:
            print(
                "blankets etc. power balance",
                totalblktetc + emultmw,
                -pthermmw_p - pnucshld,
            )
        sankey.add(
            flows=blanketsetc,
            orientations=[0, -1, 0],  # [right(in), down(in), right(out)]
            prior=1,  # DEPOSITION
            connect=(1, 0),  # Blanket/etc. --> None
            pathlengths=[0.5, 0.25, 0.0],
            labels=[None, "Energy Mult.", "Primary Heat"],
        )

        # ------------------------------------- HEAT LOSS - 3 -------------------------------------

        # Primary heat, -Gross electric power, -difference (loss)
        primary = [
            pthermmw_p + totaldivetc + totalcpetc + pnucshld,
            -pgrossmw,
            -pthermmw_p + pgrossmw - totaldivetc - totalcpetc - pnucshld,
        ]
        sankey.add(
            flows=primary,
            orientations=[0, -1, 0],  # [right(in), down(out), right(out)]
            prior=2,  # BLANKETSETC
            connect=(2, 0),  # Primary Heat --> None
            pathlengths=[0.2, 0.7, 0.4],
            labels=[None, "Gross electric", "Losses"],
        )

        # ------------------------------------ ELECTRICITY - 4 ------------------------------------

        # If net electric is +ve or -ve changes the flow organisation
        if pnetelmw >= 0:  # net electric is +ve
            # Gross electric power, -net electric power, -recirculated power
            net = [pgrossmw, -pnetelmw, -precircmw]
            sankey.add(
                flows=net,
                orientations=[0, 0, -1],  # [down(in), down(out), left(out)]
                prior=3,  # PRIMARY
                connect=(1, 0),  # Gross electric --> None
                pathlengths=[0.1, 0.25, 1.5],
                labels=[None, "Net elec.", "Recirc. Power"],
            )
        elif pnetelmw < 0:  # net electric is -ve
            # Gross electric power, -net electric power, -recirculated power
            net = [-pnetelmw, pgrossmw, -precircmw]
            sankey.add(
                flows=net,
                orientations=[0, -1, 0],  # [left(in), down(in), left(out)]
                prior=3,  # PRIMARY
                connect=(1, 1),  # Gross electric --> None
                pathlengths=[0.25, 1.0, 0.5],
                labels=["Net elec.", None, "Recirc. Power"],
            )

        # -------------------------------- RECIRCULATING POWER - 5 --------------------------------

        # Recirculated power, -Core Systems, -Heating System
        recirc = [precircmw, -pcoresystems - htpmw, -pinjwp + ppumpmw]
        # Check if difference >2 between recirculated power and the output sum
        if sum(recirc) ** 2 > 2:
            print(
                "Recirc. Power Balance",
                precircmw,
                -pcoresystems + ppumpmw - pinjwp - htpmw,
            )
        sankey.add(
            flows=recirc,
            orientations=[0, 1, 0],  # [left(in), down(out), left(out)]
            prior=4,  # NET
            connect=(2, 0),  # Recirc. Power --> None
            pathlengths=[0.1, 0.25, 0.8],
            labels=[None, "Core Systems", "Heating System"],
        )

        # --------------------------------------- LOSSES - 6 --------------------------------------

        # HCD: Heating system, -Plasma heating, -losses
        hcd = [
            pinjwp - ppumpmw,
            -p_hcd_injected_total_mw,
            -pinjwp + p_hcd_injected_total_mw + ppumpmw,
        ]
        sankey.add(
            flows=hcd,
            orientations=[0, -1, 0],  # [left(in), up(out), left(out)]
            prior=5,  # RECIRC
            connect=(2, 0),  # Heating System --> None
            pathlengths=[0.5, 0.8 + 0.5 * y_adj, 0.4],  # 'Plasma Heating' adjust
            labels=[None, "Plasma Heating", "Losses"],
        )

        # Colelcting Sankey diagram and applying a condensed layout
        diagrams = sankey.finish()
        fig.tight_layout()

        # Difference in branch tip locations for 'Plasma Heating'
        x_adj, y_adj = diagrams[0].tips[1] - diagrams[6].tips[1]

    # --------------------------------------- Label Positioning ---------------------------------------

    # Munipulating the positioning of the branch labels
    # -ve to left and down; +ve to right and up
    # pos[0] = x-axis; pos[1] = y-axis
    for d in diagrams:
        for y, t in enumerate(d.texts):
            pos = tuple(np.ndarray.tolist(d.tips[y]))
            t.set_position(pos)
            if t == diagrams[0].texts[0]:  # Fusion Power
                t.set_horizontalalignment("left")
                t.set_position((
                    pos[0] - 0.35,
                    pos[1] + 0.5 * (fusion_power / totalplasma) + 0.2,
                ))
            if t == diagrams[0].texts[2]:  # Plasma
                t.set_horizontalalignment("right")
                t.set_position((pos[0] - 0.25, pos[1]))
            if t == diagrams[1].texts[1]:  # Blanket/etc.
                t.set_horizontalalignment("right")
                t.set_position((pos[0] - 0.2, pos[1]))
            if t == diagrams[2].texts[1]:  # Energy Mult.
                t.set_position((pos[0], pos[1] - 0.3))
            if t == diagrams[2].texts[2]:  # Primary Heat
                t.set_horizontalalignment("right")
                t.set_position((pos[0] - 0.25, pos[1]))
            if t == diagrams[3].texts[1]:  # Gross Electric
                t.set_horizontalalignment("right")
                t.set_position((
                    pos[0] - 0.5 * (pgrossmw / totalplasma) - 0.1,
                    pos[1] + 0.1,
                ))
            if t == diagrams[3].texts[2]:  # Losses
                t.set_horizontalalignment("right")
                t.set_position((pos[0] - 0.2, pos[1]))
            if pnetelmw >= 1:
                if t == diagrams[4].texts[1]:  # Net electric
                    t.set_horizontalalignment("center")
                    t.set_position((pos[0], pos[1] - 0.2))
            elif pnetelmw < 1 and t == diagrams[4].texts[0]:  # Net electric
                t.set_horizontalalignment("left")
                t.set_position((pos[0] + 0.2, pos[1]))
            if t == diagrams[4].texts[2]:  # Recirc. Power
                if pnetelmw >= 1:
                    t.set_position((
                        pos[0] + 0.15,
                        pos[1] + 0.5 * (precircmw / totalplasma) + 0.2,
                    ))
                elif pnetelmw < 1:
                    t.set_horizontalalignment("left")
                    t.set_position((pos[0] + 0.2, pos[1]))
            if t == diagrams[5].texts[1]:  # Core Systems
                t.set_position((pos[0], pos[1] - 0.2))
            if t == diagrams[5].texts[2]:  # Heating System
                if pnetelmw >= 1:
                    t.set_position((
                        pos[0] + 0.15,
                        pos[1] + 0.5 * (pinjwp / totalplasma) + 0.2,
                    ))
                if pnetelmw < 1:
                    t.set_position((
                        pos[0] + 0.15,
                        pos[1] + 0.5 * (pinjwp / totalplasma) + 0.2,
                    ))
            if t == diagrams[6].texts[1]:  # Plasma Heating
                t.set_horizontalalignment("left")
                t.set_position((
                    pos[0] + 0.5 * (p_hcd_injected_total_mw / totalplasma) + 0.1,
                    pos[1] - 0.05,
                ))
            if t == diagrams[6].texts[2]:  # Losses
                t.set_horizontalalignment("left")
                t.set_position((
                    pos[0] + 0.15,
                    pos[1]
                    - 0.5 * ((pinjwp - p_hcd_injected_total_mw) / totalplasma)
                    - 0.2,
                ))
