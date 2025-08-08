import xml.etree.ElementTree as ET
import pandas as pd
import math
from pathlib import Path
from buiseigenschappen import BuisEigenschappen
from staal_beton_diagram import StaalBetonKolomCalculator
from buckling_calculator import BucklingCalculator


def parse_xml_file(xml_file_path: Path) -> tuple:

    # Ensure the file exists
    if not xml_file_path.exists():
        raise FileNotFoundError(f"The file {xml_file_path} does not exist.")

    # Parse the XML file
    tree = ET.parse(xml_file_path)
    root = tree.getroot()

    # Namespace used in the XML
    namespace = {"ns": "http://www.scia.cz"}

    return root, namespace


# Fixed values definition for calculations
def get_fixed_values(t) -> dict:
    if t == 8:
        return {
            "EI": 147337,  # inclusief 8r20
            "N_plRd": 8271,  # beton is 14Mpa aangehouden ipv 16.67 MPa
            "V_plRd": 1257,
            "Alpha": 0.21,  # curve a
            "EA_ratio": 2908 / 5842,  # ratio stalen buis tov totale EA inclusief 8r20
            "N_plrdbuis": 4154,
            "beta": 0.66,
            "M_plRd": 674,  # bij fy = 300 alleen buispaal
        }
    elif t == 12.5:
        return {
            "EI": 203552,  # inclusief 8r20
            "N_plRd": 11629,
            "V_plRd": 2800,
            "Alpha": 0.21,  # curve a
            "EA_ratio": 4507 / 7361,  # ratio stalen buis tov totale EA inclusief 8r20
            "N_plrdbuis": 7619,  # was gelijk aan 6438, maar dat is berekend met fy=300 (incorrect)
            "beta": 0.66,
            "M_plRd": 1227,  # bij fy = 355 alleen buispaal
        }
    else:
        raise ValueError("Unsupported value of t. Only t=8 or t=12.5 are allowed.")


# Extract "Slankheid staal" table data
def extract_slankheid_data_for_steelslenderness(root, namespace) -> dict:
    slankheid_data = {}
    for row in root.findall(
        ".//ns:table[@name='Steel slenderness']//ns:row", namespace
    ):
        member = row.find("ns:p0", namespace).get("v")  # Member name
        l_y = row.find("ns:p6", namespace).get("v")  # l_{y} value
        l_z = row.find("ns:p12", namespace).get("v")  # l_{z} value
        slankheid_data[member] = {"l_y": float(l_y), "l_z": float(l_z)}
    return slankheid_data


# Extract "Stability" table data
def extract_slankheid_data_for_stability(root, namespace) -> dict:
    slankheid_data = {}
    for row in root.findall(".//ns:table[@name='Stability']//ns:obj", namespace):
        member = row.find(".//ns:row/ns:p3", namespace).get("v")  # Member name
        l_y = row.find("ns:p9[@x='1']", namespace).get("v")  # l_{y} value
        l_z = row.find("ns:p14[@x='1']", namespace).get("v")  # l_{z} value
        slankheid_data[member] = {"l_y": float(l_y), "l_z": float(l_z)}
    return slankheid_data


# Extract "Results on 1D members" table data
def extract_results_1d(
    file_name: str,
    folder: Path,
    root,
    namespace,
    slankheid_data: dict,
    fixed_values: dict,
    buis: BuisEigenschappen,
    plot: bool = False,
) -> pd.DataFrame:
    EI = fixed_values["EI"]
    # N_plRd = fixed_values["N_plRd"]
    V_plRd = fixed_values["V_plRd"]
    Alpha = fixed_values["Alpha"]
    # EA_ratio = fixed_values["EA_ratio"]
    # N_plrdbuis = fixed_values["N_plrdbuis"]
    beta = fixed_values["beta"]
    # M_plRd = fixed_values["M_plRd"]

    data = []
    headers = [
        "Name",
        "dx",
        "Case",
        "Analysistype",
        "N",
        "V_y",
        "V_z",
        "M_x",
        "M_y",
        "M_z",
        "l_{y}",
        "l_{z}",
        "Max_l",
        "abs(N)",
        "Root_Vy_Vz",
        "Root_My_Mz",
        "Unity_Check_Vy_Vz",
        "Ncr",
        "lamda",
        "Xsi",
        "Unity_Check_N_Buckling",
        "k",
        "MN_in_polygon",
        "mu",
        "UC_M_ec4",
    ]

    for row in root.findall(
        ".//ns:table[@name='1D internal forces']//ns:row", namespace
    ):
        row_data = [row.find(f"ns:p{i}", namespace).get("v") for i in range(11)]
        member_name = row_data[0]  # Name column
        # Add the l_{y} and l_{z} values from the Slankheid table
        l_y_value = slankheid_data.get(member_name, {}).get("l_y", None)
        l_z_value = slankheid_data.get(member_name, {}).get("l_z", None)
        if l_y_value is not None:
            l_y_value = round(l_y_value * 1000)  # Convert to mm and round to 0 decimals
        if l_z_value is not None:
            l_z_value = round(l_z_value * 1000)  # Convert to mm and round to 0 decimals
        max_l_value = max(l_y_value, l_z_value) if l_y_value and l_z_value else None

        dx = round(float(row_data[1]), 2) if row_data[1] else 0.0

        # Divide N, V_y, V_z, M_x, M_y, M_z by 1000
        N = round(float(row_data[3]) / 1000, 2) if row_data[3] else 0.0
        V_y = round(float(row_data[4]) / 1000, 2) if row_data[4] else 0.0
        V_z = round(float(row_data[5]) / 1000, 2) if row_data[5] else 0.0
        M_x = round(float(row_data[6]) / 1000, 2) if row_data[6] else 0.0
        M_y = round(float(row_data[7]) / 1000, 2) if row_data[7] else 0.0
        M_z = round(float(row_data[8]) / 1000, 2) if row_data[8] else 0.0

        Ned = abs(N)

        # Calculate the square root of V_y^2 + V_z^2
        root_vy_vz = round(math.sqrt(V_y**2 + V_z**2), 2)

        # Calculate the square root of M_y^2 + M_z^2
        root_my_mz = round(math.sqrt(M_y**2 + M_z**2), 2)

        # Unity check: root_vy_vz / V_plRd
        unity_check_vy_vz = root_vy_vz / V_plRd if V_plRd else None

        # Unity check: N_buckling
        N_plRd = buis.NplRd()
        N_cr = math.pi**2 * EI / (max_l_value**2) * 10**6 if max_l_value else None
        Lambda = math.sqrt(buis.NplRk() / N_cr) if N_cr else None
        Psi = 0.5 * (1 + Alpha * (Lambda - 0.2) + Lambda**2) if Lambda else None
        Xhi = 1 / (Psi + math.sqrt(Psi**2 - Lambda**2)) if Psi and Lambda else None
        if Xhi and Xhi > 1:
            Xhi = 1
        unity_check_n_buckling = Ned / (Xhi * N_plRd) if Xhi and N_plRd else None

        # Unity checks: interactie buiging en normaalkracht
        # volgens art 6.2
        k = beta / (1 - Ned / N_cr)  # k-factor voor tweede orde effecten
        if k < 1:
            k = 1

        if row_data[2].startswith("NC"):
            analysis = "NonLin"
            k = 1
        else:
            analysis = "Lin"
        M_Ed = k * root_my_mz  # 2de orde moment

        r = 0  # verhoudingsgetal kopmomenten

        staalbeton = StaalBetonKolomCalculator(
            M_Ed, Ned, r, Xhi, buis.Nc(), buis.NplRd(), buis.MplRd(), buis.MmaxRd()
        )
        bool_in_polygon = staalbeton.check_point_within_polygon()
        mu, UC_M_ec4 = staalbeton.bending_moment_check()
        if plot:
            staalbeton.plot_diagrams(
                folder
                / f"{file_name} {analysis} {member_name} dx_{dx}_interactie_diagram.png"
            )
            BucklingCalculator(
                10**9,
                EI,
                max_l_value,
                1000,
                round(N_plRd, 2),
                "curve_a",
                True,
                1000 * buis.NplRk(),
            ).plot(
                folder
                / f"{file_name} {member_name} lbuc_{max_l_value}mm buckling_diagram.png"
            )

        # check the connection between the column and the concrete pile cap
        project_name = file_name.split(" ")[0]
        if project_name in [
            "N03",
            "N07",
            "N09",
            "N10",
            "N12",
            "Z01",
            "Z03",
            "Z07",
            "Z09",
            "Z10",
            "Z12",
        ]:
            paalkop_check = False
            if paalkop_check:
                if "Buispaal_fase_1" in member_name or "Buispaal_fase_2" in member_name:
                    # if "Buispaal_fase_1" in member_name:
                    if (
                        dx == 25.6
                        and analysis == "NonLin"
                        or dx == 26
                        and analysis == "NonLin"
                        or dx == 27.27
                        and analysis == "NonLin"
                    ):
                        print(
                            f"Member {member_name} with dx {dx} Ned: {Ned} kN, root_my_mz: {root_my_mz} kNm"
                        )
                        # Check if value exceedes the capacity of the concrete section
                        if root_my_mz > 200:
                            raise ValueError(
                                f"Member {member_name} with dx {dx} and analysis {analysis} has a root_my_mz value of {root_my_mz} kNm, which exceeds the limit."
                            )

        row_data[1] = dx
        row_data[3] = analysis
        row_data[4] = N
        row_data[5] = V_y
        row_data[6] = V_z
        row_data[7] = M_x
        row_data[8] = M_y
        row_data[9] = M_z
        del row_data[10]
        row_data.append(l_y_value)
        row_data.append(l_z_value)
        row_data.append(max_l_value)
        row_data.append(Ned)
        row_data.append(root_vy_vz)
        row_data.append(root_my_mz)
        row_data.append(
            round(unity_check_vy_vz, 3) if unity_check_vy_vz is not None else None
        )
        row_data.append(round(N_cr, 3) if N_cr is not None else None)
        row_data.append(round(Lambda, 3) if Lambda is not None else None)
        row_data.append(round(Xhi, 3) if Xhi is not None else None)
        row_data.append(
            round(unity_check_n_buckling, 3)
            if unity_check_n_buckling is not None
            else None
        )
        row_data.append(round(k, 3))
        row_data.append(bool_in_polygon)
        row_data.append(round(mu, 3))
        row_data.append(round(UC_M_ec4, 3))
        data.append(row_data)

    return pd.DataFrame(data, columns=headers)


def get_highest_unity_checks(results_1d_df):
    unity_check_columns = [
        "Unity_Check_Vy_Vz",
        "Unity_Check_N_Buckling",
        "MN_in_polygon",
        "UC_M_ec4",
    ]

    highest_values = {}
    for analysis_type in ["Lin", "NonLin"]:
        analysis_df = results_1d_df[results_1d_df["Analysistype"] == analysis_type]
        highest_values[analysis_type] = {
            col: (
                analysis_df[col].max()
                if col != "MN_in_polygon"
                else not analysis_df[col].eq(False).any()
            )
            for col in unity_check_columns
        }

    highest_values_df = pd.DataFrame(highest_values).T
    return highest_values_df


def check_unity_exceedance(highest_unity_checks):
    for col in highest_unity_checks.columns:
        if (highest_unity_checks[col] > 1.00).any():
            return "Conclusion: EXCEEDANCE of unity check", "EXCEEDS"
    return "Conclusion: All unity checks < 1.0", "Satisfies"


if __name__ == "__main__":

    # Define project_names based on your selection criteria
    # Example: Use a variable or condition to select which set to use
    selection = "single"  # Change this to "N", "Z", or "single" as needed

    if selection == "all":
        version = "v3.2"
        project_names = [
            "N02",
            "N03",
            "N04",
            "N05",
            "N06",
            "N08",
            "N09",
            "N11",
            "N12",
            "N13",
            "N14",
            "N15",
            "Z02",
            "Z03",
            "Z05",
            "Z06",
            "Z07",
            "Z08",
            "Z09",
            "Z10",
            "Z11",
            "Z12",
            "Z13",
            "Z14",
            "Z15",
        ]
    elif selection == "N":
        version = "v3.5"
        project_names = ["N02", "N03", "N04", "N06", "N07", "N08", "N09"]
    elif selection == "Z":
        version = "v3.5"
        project_names = ["Z06", "Z10", "Z15"]
    elif selection == "single":
        version = "v3.5"
        project_names = ["N01"]
    else:
        project_names = []
    folder_path = Path(
        "P:/1329/132975/WIP/PR-0314-Voorafgaande_Werken/PR-0488 - Voorafgaande werken TD16/31_paaljukken/Herbeschouwing_paaljukken_2025_01_01"
    )
    folder_path = folder_path / f"{version} - rekenuitvoer"
    for project_name in project_names:
        file_name = f"{project_name} - {version}_AJOR"
        file_name_extension = file_name + ".xml"
        xml_file_path = folder_path / file_name_extension

        # Parse the XML file
        root, namespace = parse_xml_file(xml_file_path)

        if project_name in [
            "N02",
            "N05",
            "N08",
            "N11",
            "N14",
            "Z02",
            "Z05",
            "Z08",
            "Z11",
            "Z14",
        ]:
            t = 12.5
            buis = BuisEigenschappen(t=12.5, fyd=355)

        else:
            t = 8
            buis = BuisEigenschappen(t=8, fyd=355)

        print(f"project: {project_name}    t = {t}mm")

        fixed_values = get_fixed_values(t)

        # Main execution
        # slankheid_data = extract_slankheid_data_for_stability(root, namespace)
        slankheid_data = extract_slankheid_data_for_steelslenderness(root, namespace)
        # Create folder if it does not exist
        output_folder = folder_path / file_name
        # output_folder.mkdir(parents=True, exist_ok=True)
        results_1d_df = extract_results_1d(
            file_name,
            output_folder,
            root,
            namespace,
            slankheid_data,
            fixed_values,
            buis,
            False,
        )
        # Display or save the table
        # print(results_1d_df)
        csv_name = file_name + "_results.csv"
        csv_path = folder_path / csv_name
        results_1d_df.to_csv(csv_path, index=False)  # Optional: Save to CSV

        highest_unity_checks = get_highest_unity_checks(results_1d_df)
        print(highest_unity_checks)
        result, tag = check_unity_exceedance(highest_unity_checks)
        print(result)
        csv_output_path = folder_path / f"{file_name}_summary_{tag}.csv"
        highest_unity_checks.to_csv(csv_output_path, index=True)
