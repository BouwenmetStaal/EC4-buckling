
import matplotlib.pyplot as plt
from typing import Tuple
from matplotlib.figure import Figure
from scipy.optimize import fsolve
from Afgekapte_cirkel import TruncatedCircle, Circle
from staal_beton_diagram import StaalBetonKolomCalculator
from buckling_calculator import BucklingCalculator
from dataclasses import dataclass
import math

@dataclass
class BuisEigenschappen:
    t: float
    diameter: float = 559
    fyd: float = 300 #for t=8 mm, zie tabel 6.3 EC4
    fcd: float = 16.67 #C25/35
    fsd: float = 400/1.15 #MPa rekenwaarde B400B wapeningsstaal
    wapening_diameter: float = 20 #mm
    staven_wapening: int = 8 #aantal staven wapening
    @property
    def staal_r(self):
        return 0.5*self.diameter

    @property
    def beton_r(self):
        return self.staal_r - self.t

    

    # Define the function to solve
    def equation(self,h_n):
        area_beton_strip = TruncatedCircle(self.beton_r, h_n).area_of_strip()
        area_steel_strip = TruncatedCircle(self.staal_r, h_n).area_of_strip() - TruncatedCircle(self.beton_r, h_n).area_of_strip()
        return area_beton_strip*self.fcd + area_steel_strip*2*self.fyd - Circle(self.beton_r).area()*self.fcd

    def solve_hn(self):
    # Solve the equation
        solutions = fsolve(self.equation, [100])  # Initial guess [1]
        h_n = solutions[0]
        return h_n


    def MmaxRd(self):
        Wpl_staal = Circle(self.staal_r).plastic_section_modulus() - Circle(self.beton_r).plastic_section_modulus()
        Wpl_beton = Circle(self.beton_r).plastic_section_modulus()

        MmaxRd = (0.5 * Wpl_beton * self.fcd + Wpl_staal * self.fyd) / 1e6  # wapening niet meegenomen
        return MmaxRd

    def MplRd(self):
        h_n = self.solve_hn()

        area_beton_strip = TruncatedCircle(self.beton_r, h_n).area_of_strip()
        area_steel_strip = TruncatedCircle(self.staal_r, h_n).area_of_strip() - TruncatedCircle(self.beton_r, h_n).area_of_strip()
        a = area_beton_strip*self.fcd
        b = area_steel_strip*2*self.fyd 
        c = Circle(self.beton_r).area()*self.fcd
        result = a + b - c # should be zero

        Wpl_staal_n = TruncatedCircle(self.staal_r, h_n).plastic_section_modulus() - TruncatedCircle(self.beton_r, h_n).plastic_section_modulus()
        Wpl_beton_n = TruncatedCircle(self.beton_r, h_n).plastic_section_modulus()

        MnRd = (0.5 * Wpl_beton_n * self.fcd + Wpl_staal_n * self.fyd) / 1e6  # wapening niet meegenomen
        MplRd = self.MmaxRd() - MnRd
        return MplRd
    
    def NplRd(self):
        #NEN-EN 1994-1-1, art. 6.7.3.2
        """
        Calculate the design axial load capacity (Npl,Rd) of a composite column section.
        Returns:
            float: The design axial load capacity (Npl,Rd) in kN.
        """
        alpha = 1.0
        A_staal = Circle(self.staal_r).area() - Circle(self.beton_r).area()
        A_wapening = self.staven_wapening *Circle(0.5*self.wapening_diameter).area()#8 staven van 20 mm
        A_beton = Circle(self.beton_r).area() - A_wapening
        NplRd = (A_staal*self.fyd + (A_beton-A_wapening)*alpha*self.fcd+ A_wapening*alpha*self.fsd)/1000 #kN
        return NplRd #kN
    
    def NplRk(self):
        """	
        Calculate the characteristic axial load capacity (Npl,Rk) of a composite column section.
        Returns:
            float: The characteristic axial load capacity (Npl,Rk) in kN.
        """
        #NEN-EN 1994-1-1, art.
        alpha = 1.0
        A_staal = Circle(self.staal_r).area() - Circle(self.beton_r).area()
        A_wapening = self.staven_wapening *Circle(0.5*self.wapening_diameter).area()#8 staven van 20 mm
        A_beton = Circle(self.beton_r).area() - A_wapening
        NplRd = (A_staal*(355) + (A_beton-A_wapening)*alpha*(25)+ A_wapening*alpha*500)/1000 #kN
        return NplRd #kN
    
    def Nc (self):
        #NEN-EN 1994-1-1, art.
        alpha = 1.0
        A_wapening = self.staven_wapening *Circle(0.5*self.wapening_diameter).area()#8 staven van 20 mm
        A_beton = Circle(self.beton_r).area() - A_wapening #oppervlakte beton in mm2

        Nc = (A_beton * alpha *self.fcd)/1000 #kN
        return Nc
    
    def staalbijdragefactor(self):
        staal_A = Circle(self.staal_r).area() - Circle(self.beton_r).area()
        staalbijdragefactor = (staal_A*self.fyd)/ (self.NplRd()*1000)
        criterion_text = "Can the column be considered as a steel-concrete column?"
        latex_formula = r"$\delta = \frac{A_{staal} \cdot f_{yd}}{N_{pl,Rd} \cdot 1000}$"
        if staalbijdragefactor < 0.2:
            text = f"$\delta$:{staalbijdragefactor:.2f} < 0.2\n No, the column must be considered as a reinforced concrete column"
            uitkomst = False
        elif staalbijdragefactor > 0.9:
            text = f"$\delta$:{staalbijdragefactor:.2f} > 0.9\n No, the column must be considered as a fully steel column"
            uitkomst = False
        else:
            text = f"$\delta$: 0.2 < {staalbijdragefactor:.2f} < 0.9\n Yes, the column must be considered as a steel-concrete column"
            uitkomst = True

        combined_text = f"{criterion_text}\n{latex_formula}\n{text}"
        return combined_text, uitkomst
    
    def plastic_analysis_criterion(self):
        dt = self.diameter / self.t
        criterion = 90 * (235 / self.fyd)
        criterion_text = "Can the plastic capacity be considered?"
        latex_formula = r"$\frac{d}{t} < 90 \cdot \frac{235}{f_{yd}}$"
        if dt < criterion:
            plastic = r"$\frac{d}{t} = "+f"{dt:.2f}"+r" < 90 \cdot \frac{235}{"+f"{self.fyd}"+r"}$ = "+f"{criterion:.2f}"+f" \n Yes, Plastic analysis allowed"
            uitkomst = True
        else:
            plastic = r"$\frac{d}{t} = "+f"{dt:.2f}"+r" \geq 90 \cdot \frac{235}{"+f"{self.fyd}"+r"}$ = "+f"{criterion:.2f}"+f" \n No, Does not fulfill d/t criterion for plastic analysis"
            uitkomst = False
        combined_text = f"{criterion_text}\n{latex_formula}\n{plastic}"
        return combined_text, uitkomst
    
    def criteria_info(self):
        text_a, bool_a = self.staalbijdragefactor()
        text_b, bool_b = self.plastic_analysis_criterion()
        combined_text = f"The cross-section must meet the following two requirements:\n\n\n{text_a}\n\n\n{text_b}\n"
        if bool_a and bool_b:
            return combined_text, True
        else:
            return combined_text, False

    

    def figure_criteria_info(self) -> Tuple[Figure, bool]:
        """
        Generates a matplotlib figure displaying criteria information as centered text.

        Calls the `criteria_info` method to retrieve the criteria text and a boolean value.
        Creates a figure with the criteria text displayed in the center, with axes turned off.

        Returns:
            Tuple[Figure, bool]: The generated figure containing the criteria information and the boolean value.
        """
        text_a, bool_a = self.criteria_info()
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.axis('off')
        ax.text(0.5, 0.5, text_a, ha='center', va='center', wrap=True, fontsize=12)
        return fig, bool_a

    def EIeff(self):
        """
        Calculate the effective flexural stiffness (EIeff) of a composite column section.
        This method calculates the effective flexural stiffness according to NEN-EN 1994-1-1, art.6.7.3.3 (3).
        Returns:
            float: The effective flexural stiffness (EIeff) in *10^9 Nmm^2.
        """
        
        Ke = 0.6
        E_staal = 210000   
        E_beton = 10500
        I_staal = Circle(self.staal_r).moment_of_inertia() - Circle(self.beton_r).moment_of_inertia()
        I_beton = Circle(self.beton_r).moment_of_inertia()
        hoh = 200
        arm_1 = hoh * math.sin(0.125 * math.pi)
        arm_2 = hoh * math.sin(0.375 * math.pi)
        wap = Circle(0.5*self.wapening_diameter)
        I_wap_1 = 4*(wap.moment_of_inertia() + wap.area() * arm_1**2) #steiners rule
        I_wap_2 = 4*(wap.moment_of_inertia() + wap.area() * arm_2**2) #steiners rule
        I_wapening = I_wap_1 + I_wap_2
        
        EIeff = E_staal*I_staal  + E_staal*I_wapening+ Ke* E_beton*I_beton
        return EIeff / 1e9
    
    def create_MN_curve(self):
        list_Mrd = []
        list_Nrd = []
        import numpy as np
        for x in np.linspace(0, 2*self.staal_r, 100):
            print(f'{x:.0f} mm')
            min_y = -self.staal_r
            neutral_axis = x-self.staal_r
            max_y = self.staal_r
            beton = TruncatedCircle(self.beton_r, 0)
            staal = TruncatedCircle(self.staal_r, 0)
            Mrd_beton = self.fcd * beton.plastic_section_modulus_assymetric(neutral_axis, max_y)
            Mrd_staal_druk = self.fyd * (staal.plastic_section_modulus_assymetric(neutral_axis, max_y)- beton.plastic_section_modulus_assymetric(neutral_axis, max_y))
            Mrd_staal_trek = self.fyd * (staal.plastic_section_modulus_assymetric(min_y, neutral_axis) - beton.plastic_section_modulus_assymetric(min_y, neutral_axis))
            Mrd = Mrd_beton + Mrd_staal_druk - Mrd_staal_trek
            Nrd_beton = self.fcd * beton.area_of_assymmetric_strip(neutral_axis, max_y)
            Nrd_staal_druk = self.fyd * (staal.area_of_assymmetric_strip(neutral_axis, max_y) - beton.area_of_assymmetric_strip(neutral_axis, max_y))
            Nrd_staal_trek = self.fyd * (staal.area_of_assymmetric_strip(min_y, neutral_axis) - beton.area_of_assymmetric_strip(min_y, neutral_axis))
            Nrd = Nrd_beton + Nrd_staal_druk - Nrd_staal_trek
            list_Mrd.append(Mrd*1e-6)
            list_Nrd.append(Nrd*1e-3)
        return list_Mrd, list_Nrd
    
    def plot_MN_curve(self, filename):
        list_Mrd, list_Nrd = self.create_MN_curve()
        plt.plot(list_Mrd, list_Nrd)
        plt.xlabel('M [kNm]')
        plt.ylabel('N [kN]')
        plt.axhline(0, color='black',linewidth=0.5)
        plt.axvline(0, color='black',linewidth=0.5)
        plt.grid(color = 'gray', linestyle = '--', linewidth = 0.5)
        plt.savefig(filename)
        plt.show()
        return plt
    
    def plot_interactie_diagram(self, MEd, Ned, r, w_buc, path:str):
        list_Mrd, list_Nrd = self.create_MN_curve()
        staalbeton = StaalBetonKolomCalculator(MEd, Ned, r, w_buc, self.Nc(), self.NplRd(), MplRd, MmaxRd, list_Nrd, list_Mrd)
        fig = staalbeton.plot_diagrams(path)

    
    
    def calculate_I_pipe(self, t):
        I = math.pi / 64 * (self.diameter**4 - (self.diameter - 2 * t)**4)
        return I

            
    
    
    



if __name__ == "__main__":
    Ned = 2000 #kN
    MEd = 2 #kNm
    r = 1
    l_buc = 7000
    fy = 300
    diameter = 559
    fcd = 25/1.5
 
    w_buc = 0.86

    buis = BuisEigenschappen(8, diameter, fy, fcd, 0,0,0)
    buis = BuisEigenschappen(8, diameter, fy, fcd)

    fig, _ = buis.figure_criteria_info()
    fig.show()

    #buis.plot_MN_curve('testplots/buispaal_MN_curve.png')

    MmaxRd = buis.MmaxRd()
    MplRd = buis.MplRd()
    EIeff = buis.EIeff()

    buckling = BucklingCalculator(10**9, round(EIeff), l_buc, 1000, round(buis.NplRd(),2), "curve_a")
    buckling.plot('testplots/buispaal_buckling.png')
    w_buc = buckling.reduction_factor

    print(f"{MmaxRd:.2f} kNm")
    print(f"{MplRd:.2f} kNm")
    print(f"EIeff: {EIeff:.0f} ")
    EIeff = BuisEigenschappen(8,300).EIeff()
    print(f"EIeff: {EIeff:.0f} ")
    list_Mrd, list_Nrd = buis.create_MN_curve()
    buis.plot_interactie_diagram(MEd, Ned,r,w_buc, 'testplots/buispaal_interactie_diagram2.png')

    #staalbeton = StaalBetonKolomCalculator(MEd, Ned, r, w_buc, buis.Nc(), buis.NplRd(), MplRd, MmaxRd, list_Nrd, list_Mrd)

    #fig = staalbeton.plot_diagrams('testplots/buispaal_interactie_diagram.png')
