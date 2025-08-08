import math
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from dataclasses import dataclass, field
from pathlib import Path

@dataclass
class BucklingCalculator:
    """
    BucklingCalculator is a dataclass for calculating and visualizing the buckling resistance of columns according to Eurocode 3 (EC3) and Eurocode 4 (EC4).
    Attributes:
        E (float): Young's modulus of the material [MPa or N/mm²].
        I (float): Second moment of area (moment of inertia) [mm⁴].
        l_buc (float): Buckling length of the column [mm].
        A (float): Cross-sectional area [mm²].
        fy (float): Yield strength of the material [MPa or N/mm²].
        selected_curve (str): The selected buckling curve key (e.g., "curve_a", "curve_b", etc.).
        alpha_values (dict): Dictionary mapping buckling curve names to their imperfection factors (α).
        selected_alpha (float): Imperfection factor (α) for the selected curve (set in __post_init__).
        staalbetonkolom (bool): If True, use EC4 (composite column) formulas; otherwise, use EC3.
        N_pl_rk (float): Plastic resistance [N] in  of the composite column (used if staalbetonkolom is True).

    Usage:
        Create an instance of BucklingCalculator with the required parameters, then call main() or plot() to perform calculations and visualize results.
    """
    E: float
    I: float
    l_buc: float
    A: float
    fy: float
    selected_curve: str
    staalbetonkolom: bool = False
    N_pl_rk: float = None
    alpha_values: dict = field(default_factory=lambda: {
        "curve_a": 0.21,
        "curve_b": 0.34,
        "curve_c": 0.49,
        "curve_d": 0.76,
        "curve_e": 1.13
    })
    selected_alpha: float = field(init=False)
    
    
    def __post_init__(self):
        self.selected_alpha = self.alpha_values[self.selected_curve]

    def calculate_reduction_factor_x(self, alpha, lambda_rel):
        phi = 0.5 * (1 + alpha * (lambda_rel - 0.2) + lambda_rel**2)
        x = 1 / (phi + math.sqrt(phi**2 - lambda_rel**2))
        return min(x, 1.0)

    def calculate_euler_buckling_load(self):
        return (math.pi**2 * self.E * self.I) / (self.l_buc**2)

    def calculate_relative_slenderness(self, N_cr):
        if self.staalbetonkolom:
            return math.sqrt(self.N_pl_rk / N_cr)
        else:
            return math.sqrt((self.A * self.fy) / N_cr)

    def calculate_axial_buckling_resistance(self, chi):
        return chi * self.A * self.fy

    def plot_results(self, lambda_rel, reduction_factor, N_cr, N_b)->plt.Figure:
        if lambda_rel < 3.0:            
            lambda_rel_values = [i * 0.01 for i in range(301)]
        else:
            lambda_rel_values = [i * 0.01 for i in range(int(lambda_rel * 100) + 1)]
        fig = plt.figure(figsize=(10, 10))
        gs = GridSpec(2, 1, height_ratios=[1, 2], figure=fig)

        ax1 = fig.add_subplot(gs[0])
        ax1.axis("off")
        if self.staalbetonkolom:
            ax1.set_title("EC4 Buckling Formulas")
            formula_text = (
                r"$N_{cr} = \frac{\pi^2 \cdot EI_{eff}}{l_{buc}^2}, \quad N_{cr} = \frac{\pi^2 \cdot "+f"{self.E}"+r" \cdot "+f"{self.I}"+r"}{"f"{self.l_buc}^2"+r"}"+" = " + f"{N_cr:.0f} \, N = {N_cr/1000:.1f} \, kN$" + "\n"
                r"$\lambda_{rel} = \sqrt{\frac{N_{pl,Rk}}{N_{cr}}}, \quad \lambda_{rel} = \sqrt{\frac{"+f"{self.N_pl_rk}"+r"}{N_{cr}}} = " + f"{lambda_rel:.3f}$" + "\n"
                r"$\varphi = 0.5 \cdot (1 + \alpha \cdot (\lambda_{rel} - 0.2) + \lambda_{rel}^2) $" + "\n"
                r"$\chi = \frac{1}{\varphi + \sqrt{\varphi^2 - \lambda_{rel}^2}} = " + f"{reduction_factor:.3f}$" + "\n"
                r"$N_b = \frac{\chi \cdot N_{pl,Rd}}{\gamma_{M1}}, \quad N_b = \frac{"+f"{reduction_factor:.3f}"+" \cdot "+f"{self.A}"+" \cdot "+f"{self.fy}"+"}{1.0} = " + f"{N_b:.0f} \, N = {N_b/1000:.1f} \, kN$"
            )
        else:
            ax1.set_title("EC3 Buckling Formulas")
            formula_text = (
                r"$N_{cr} = \frac{\pi^2 \cdot EI}{l_{buc}^2}, \quad N_{cr} = \frac{\pi^2 \cdot "+f"{self.E}"+r" \cdot "+f"{self.I}"+r"}{"f"{self.l_buc}^2"+r"}"+" = " + f"{N_cr:.0f} \, N = {N_cr/1000:.1f} \, kN$" + "\n"
                r"$\lambda_{rel} = \sqrt{\frac{A \cdot f_y}{N_{cr}}}, \quad \lambda_{rel} = \sqrt{\frac{"+f"{self.A}"+" \cdot "+f"{self.fy}"+r"}{N_{cr}}} = " + f"{lambda_rel:.3f}$" + "\n"
                r"$\varphi = 0.5 \cdot (1 + \alpha \cdot (\lambda_{rel} - 0.2) + \lambda_{rel}^2) $" + "\n"
                r"$\chi = \frac{1}{\varphi + \sqrt{\varphi^2 - \lambda_{rel}^2}} = " + f"{reduction_factor:.3f}$" + "\n"
                r"$N_b = \frac{\chi \cdot A \cdot f_y}{\gamma_{M1}}, \quad N_b = \frac{"+f"{reduction_factor:.3f}"+" \cdot "+f"{self.A}"+" \cdot "+f"{self.fy}"+"}{1.0} = " + f"{N_b:.0f} \, N = {N_b/1000:.1f} \, kN$"
            )

        ax1.text(0.5, 0.5, formula_text, fontsize=12, ha="center", va="center", multialignment="center")
        

        ax2 = fig.add_subplot(gs[1])
        for curve, alpha in self.alpha_values.items():
            reduction_factors = [self.calculate_reduction_factor_x(alpha, lambda_rel) for lambda_rel in lambda_rel_values]
            ax2.plot(lambda_rel_values, reduction_factors, label=f"{curve} (α={alpha})")

        ax2.scatter([lambda_rel], [reduction_factor], color="red", zorder=5, label=f"$({self.selected_curve} "+r", \quad \lambda_{rel}="+f"{lambda_rel:.3f}"+r", \quad \chi="+f"{reduction_factor:.3f})$")
        ax2.plot([lambda_rel, lambda_rel], [0, reduction_factor], color="gray", linestyle="--", linewidth=1)
        ax2.plot([0, lambda_rel], [reduction_factor, reduction_factor], color="gray", linestyle="--", linewidth=1)
        ax2.text(lambda_rel, 0, "$\lambda_{rel} = "+f"{lambda_rel:.3f}$", ha="center", va="bottom", fontsize=10, color="blue")
        ax2.text(0, reduction_factor, "$\chi = "+f"{reduction_factor:.3f}$", ha="left", va="center", fontsize=10, color="blue")
        ax2.set_title("Reduction Factor $\chi$ vs Relative Slenderness $\lambda_{rel}$")
        ax2.set_xlabel("Relative Slenderness $\lambda_{rel}$")
        ax2.set_ylabel("Reduction Factor $\chi$")
        ax2.legend()
        ax2.grid(True)
        plt.tight_layout()
        return fig
    
    @property
    def reduction_factor(self):
        N_cr = self.calculate_euler_buckling_load()
        lambda_rel = self.calculate_relative_slenderness(N_cr)
        reduction_factor = self.calculate_reduction_factor_x(self.selected_alpha, lambda_rel)
        return reduction_factor
    
    @property
    def N_b(self):
        reduction_factor = self.reduction_factor
        return self.calculate_axial_buckling_resistance(reduction_factor)

    def main(self):
        N_cr = self.calculate_euler_buckling_load()

        lambda_rel = self.calculate_relative_slenderness(N_cr)

        reduction_factor = self.calculate_reduction_factor_x(self.selected_alpha, lambda_rel)

        N_b = self.calculate_axial_buckling_resistance(reduction_factor)

        return self.plot_results(lambda_rel, reduction_factor, N_cr, N_b)
    
    def plot(self, path:Path):
        figure = self.main()
        figure.savefig(path)
        plt.close(figure)
        

if __name__ == "__main__":
    E = 210000  # Flexural rigidity in N/mm^2
    I = 3692 * 10**4  # mm^4
    l_buc = 10000  # Buckling length in mm
    A = 5383  # Cross-sectional area in mm^2
    fy = 355  # Yield strength in MPa
    selected_curve = "curve_b"

    calculator = BucklingCalculator(E, I, l_buc, A, fy, selected_curve, staalbetonkolom=True, N_pl_rk=1000)
    figure = calculator.main()
    figure.show()
    pass
