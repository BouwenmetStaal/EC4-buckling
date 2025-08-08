from blueprints.structural_sections.steel.steel_cross_sections.standard_profiles.chs import (CHS,
)
from blueprints.codes.eurocode.en_1993_1_1_2005.chapter_3_materials.table_3_1 import (
    SteelStrengthClass,
)
from blueprints.materials.steel import SteelMaterial

from blueprints.structural_sections.steel.steel_cross_sections.i_profile import (
    ISteelProfile,
)
from buckling_calculator import BucklingCalculator

from buiseigenschappen import BuisEigenschappen
from staal_beton_diagram import StaalBetonKolomCalculator
from pathlib import Path
import math
steel_material = SteelMaterial(steel_class=SteelStrengthClass.S355)

#Dubbele HEB snede, toetsing om de sterke as

for chs_section in CHS:
    l_buc = 5600  # Buckling length in mm
    selected_curve = "curve_c"

    diameter = chs_section.diameter
    thickness = chs_section.thickness - 0.5  # Adjusted thickness for corrosion
    fyd = 355
   
    buis = BuisEigenschappen(t=thickness, diameter=diameter, fyd=fyd)

    r = 1
        
    M_Ed = 650
    Ned = 500
    Xhi = 0.8
    staalbeton = StaalBetonKolomCalculator(
        M_Ed, Ned, r, Xhi, buis.Nc(), buis.NplRd(), buis.MplRd(), buis.MmaxRd()
    )
    
    m = buis.MplRd()
    
    if math.isnan(m):
        continue

    mprd = math.floor(m)
    if mprd > M_Ed and mprd < M_Ed + 100:

        dt = diameter/thickness
        criterion = 90*(235/fyd)
        if dt < criterion:
            plastic = "Plastic analysis allowed"
        else:
            plastic = "Does not fulfill d/t criterion for plastic analysis"

        print(f"CHS Section: {chs_section.name} MplRd: {mprd} kNm, {plastic}")

        folder = Path(r"C:\Users\AJOR\OneDrive - Witteveen+Bos\LEKA_palen")

        staalbeton.plot_diagrams(
            folder
            / f"{chs_section.name}_{mprd}_kNm_interaction_and_stability_diagram.png"
        )
    # BucklingCalculator(
    #     10**9,
    #     EI,
    #     row.Max_l,
    #     1000,
    #     round(buis.NplRd(), 2),
    #     "curve_a",
    #     True,
    #     round(1000*buis.NplRk(), 2),
    # ).plot(
    #     folder
    #     / f"{file_name}_{member_name}_lbuc_{row.Max_l}mm_buckling_diagram.png"
    # )
