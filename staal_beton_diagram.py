from matplotlib.path import Path
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

@dataclass
class StaalBetonKolomCalculator:
    Med: float
    Ned: float
    r: float
    w_buc: float
    Nc: float
    NplRd: float
    M_plRd: float
    M_maxRd: float
    N_list_detailed: list = None
    M_list_detailed: list = None

    @property
    def N_list(self):
        return [self.NplRd, self.Nc, 0.5*self.Nc , 0]
    
    @property
    def M_list(self):
        return [0, self.M_plRd, self.M_maxRd, self.M_plRd]


    def calculate_parameters(self):
        Nc = self.Nc
        NplRd = self.NplRd
        Ned = self.Ned
        r = self.r
        w_buc = self.w_buc
        k_c = Nc/NplRd
        k_d = Ned/NplRd
        k_n = 0.25*(1-r)

        if k_d < k_c:
            mu_d = 1
            mu_k = ((1-w_buc)/(1-k_c)) * mu_d
            if k_d < k_n:
                mu = 1
            else:
                mu = 1 - (k_d - k_n)/(w_buc - k_n) * mu_k    
        else:   
            mu_d = 1 - (k_d - k_c)/(1 - k_c)
            mu_k = ((1-w_buc)/(1-k_d)) * mu_d
            x_coord = (k_d - k_n)/(w_buc - k_n) * mu_k
            mu = mu_d - x_coord

        return k_c, k_d, k_n, mu_k, mu_d, mu
    
    def bending_moment_check(self):
        #NEN-EN 1994-1-1, art. 6.7.3.6
        _, _, _, _, _, mu = self.calculate_parameters()
        if mu < 0:
            mu = 0.0000001

        alpha_m = 0.9 #voor de staalsoorten S235 t/m S355;
        M_plRd = self.M_plRd
        M_Rd = mu * M_plRd
        UC = (1/alpha_m) * self.Med / M_Rd
        return mu, UC




    def check_point_within_polygon(self):
        Med = self.Med
        Ned = self.Ned
        # Create a polygon from the points in M_list and N_list
        polygon_points = np.column_stack((self.M_list + [0], self.N_list + [0]))
        polygon = Path(polygon_points)

        # Check if the point is within the polygon or on the border
        if polygon.contains_point((Med, Ned)) or polygon.contains_point((Med, Ned), radius=-1e-9):
            return True
        else:
            return False


    def plot_diagrams(self, path:str):
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
        ax1 = self.plot_interaction_diagram(ax1)
        ax2 = self.plot_stability_diagram(ax2)
        ax3 = self.plot_formula(ax3)
        plt.tight_layout()
        fig.savefig(path)
        plt.close(fig)

    def plot_interaction_diagram(self, ax1)->None:
        M_list = self.M_list
        N_list = self.N_list
        MEd = self.Med
        Ned = self.Ned
        if self.check_point_within_polygon():
            point_color = 'green'
        else:
            point_color = 'red'
        
        # First graph
        ax1.plot(M_list, N_list, marker='o')
        ax1.set_xlabel('Moment (kNm)')
        ax1.set_ylabel('$N_{rd}$', rotation=0, labelpad=20)
        ax1.set_title('Interaction Diagram')
        Nmax = N_list[0]
        if self.Ned < Nmax:
            ax1.set_ylim(0, Nmax)
        else:
            ax1.set_ylim(0, self.Ned+100)
        Mmax = M_list[2]
        if self.Med < Mmax:
            ax1.set_xlim(0, Mmax)
        else:
            ax1.set_xlim(0, self.Med+20)

        # Plot the point (MEd, Ned) and check if it is within the boundaries
        ax1.plot(MEd, Ned, 'ro')  # Plot the point in red
        ax1.annotate(f'({MEd}, {Ned})', (MEd, Ned), textcoords="offset points", xytext=(0,10), ha='center', color='red')

        

        # Plot the point with the determined color
        ax1.plot(MEd, Ned, marker='o', color=point_color)
        ax1.annotate(f'({MEd}, {Ned})', (MEd, Ned), textcoords="offset points", xytext=(0,10), ha='center', color=point_color)
        # Add vertical and horizontal lines
        ax1.axhline(y=Ned, xmin=0, xmax=MEd/ax1.get_xlim()[1], color=point_color, linestyle='--', linewidth=0.7)
        ax1.axvline(x=MEd, ymin=0, ymax=Ned/ax1.get_ylim()[1], color=point_color, linestyle='--', linewidth=0.7)
        # Annotate the coordinates of the nodes
        for i, (M, N) in enumerate(zip(M_list, N_list)):
            if i == 0:
                ax1.annotate(f'$N_{{pl,Rd}} = {N:.0f} kN$', (M, N), textcoords="offset points", xytext=(0,10), ha='center')
            else:
                ax1.annotate(f'({M:.0f} kNm, {N:.0f} kN)', (M, N), textcoords="offset points", xytext=(0,10), ha='center')

        # Add horizontal dashed lines
        for y, x in zip(N_list, M_list):
            ax1.plot([0, x], [y, y], color='black', linestyle='--', linewidth=0.7)

        if self.N_list_detailed is not None and self.M_list_detailed is not None:
            ax1.plot(self.M_list_detailed, self.N_list_detailed, color='blue', linestyle='-', linewidth=0.7)

        return ax1


    def plot_stability_diagram(self, ax2):
        k_c, k_d, k_n, mu_k, mu_d, mu = self.calculate_parameters()
        w_buc = self.w_buc
        chi_list = [1, k_c, 0]
        mu_list = [0,1,1]
        # Second graph (example: same data but with different style)
        ax2.plot(mu_list, chi_list, marker='x', linestyle='-', color='r')
        ax2.set_xlabel(r'$\frac{M_{rd}}{M_{pl,Rd}}$')
        ax2.set_ylabel(r'$\frac{N_{rd}}{N_{pl,Rd}}$', rotation=0, labelpad=20)
        ax2.set_title('Stability Diagram')
        ax2.set_ylim(0, 1)  # Set the y-axis limits to start at 0 and end at NplRd
        ax2.set_xlim(0, 1)  # Set the x-axis limits to start at 0 and end at M_d

        # Add horizontal dashed line
        ax2.axhline(y=k_c, color='black', linestyle='--', linewidth=0.7)
        ax2.annotate(f'$k_c = {k_c:.2f}$', xy=(0, k_c), xytext=(-30, 0), textcoords='offset points', ha='center')
        ax2.annotate(f'$k_n = {k_n:.2f}$', xy=(0, k_n), xytext=(-30, 0), textcoords='offset points', ha='center')
        ax2.annotate(f'$k_d = {k_d:.2f}$', xy=(0, k_d), xytext=(-30, 0), textcoords='offset points', ha='center')
        ax2.plot([0, mu_k], [w_buc, w_buc], color='black', linestyle='--', linewidth=0.7)
        # Add vertical dashed line
        ax2.plot([mu_k, mu_k], [0, w_buc], color='black', linestyle='--', linewidth=0.7)
        ax2.annotate(f'$w_{{buc}} = {w_buc:.2f}$', xy=(0, w_buc), xytext=(-30, 0), textcoords='offset points', ha='center')
        ax2.annotate(f'$\mu_{{k}} = {mu_k:.2f}$', xy=(mu_k, 0.05), xytext=(0, -10), textcoords='offset points', ha='center')
        ax2.annotate(f'$\mu_{{d}} = {mu_d:.2f}$', xy=(mu_d, 0.05), xytext=(0, -10), textcoords='offset points', ha='center')
        # Add vertical dashed line from (mu_d, k_d) to (mu_d, 0)
        ax2.plot([mu_d, mu_d], [k_d, 0], color='blue', linestyle='--', linewidth=0.7)
        # Add diagonal dashed line from (0, k_n) to (mu_k, w_buc)
        ax2.plot([0, mu_k], [k_n, w_buc], color='black', linestyle='--', linewidth=0.7)

        x_coord = (k_d - k_n)/(w_buc - k_n) * mu_k
        if x_coord < 0:
            x_coord = 0
        ax2.annotate(f'({x_coord:.2f}, {k_d:.2f})', xy=(x_coord, k_d), xytext=(10, 10), textcoords='offset points', ha='center', color='blue')
        # Add diagonal dashed line from previous coordinate to (1.0, k_c)
        ax2.annotate('', xy=(x_coord, k_d), xytext=(mu_d, k_d), arrowprops=dict(arrowstyle='<->', color='blue'))
        distance = mu_d - x_coord
        if distance != mu:
            raise ValueError(f'The distance between the point ({x_coord:.2f}, {k_d:.2f}) and the point ({mu_d:.2f}, {k_d:.2f}) is not equal to mu ({mu:.2f})')
        ax2.annotate(f'$\mu = {distance:.2f}$', xy=((x_coord + 1.0) / 2, k_d), xytext=(0, 10), textcoords='offset points', ha='center', color='blue')
        return ax2

    def plot_formula(self, ax):
        mu, uc = self.bending_moment_check()
        ax.axis('off')
        formula1 = r"$\frac{M_{IImax}}{\mu_{d} * M_{Pl,Rd}} < \alpha_m$"
        formula2 = r"$\frac{" + f"{self.Med:.0f} kNm" + "}{" + f"{mu:.2f}" + " * " + f"{self.M_plRd:.0f} kNm" + "} < 0.9$"
        formula3 = f"UC = {uc:.2f}"
        if uc > 1:
            color = 'red'
        else:
            color = 'green'
        ax.text(0.5, 0.6, formula1, fontsize=12, ha='center', va='center')
        ax.text(0.5, 0.4, formula2, fontsize=12, ha='center', va='center', color=color)
        ax.text(0.5, 0.2, formula3, fontsize=12, ha='center', va='center', color=color)
        return ax


if __name__ == "__main__":
    Ned = 1000 #kN
    MEd = 150 #kNm
    r = 1

    fcd = 20/1.5 #MPa
    alpha = 0.85
    A_c = 300**2 - 5380 #oppervlakte beton in mm2
    
    Nc = (A_c * alpha *fcd)/1000 #kN
    NplRd = 3034 #kN

    w_buc = 0.86
    M_d = 206 #kNm
    M_n = 14
    M_plRd = M_d - M_n

    staalbeton = StaalBetonKolomCalculator(MEd, Ned, r, w_buc, Nc, NplRd, M_plRd, M_d)
    k_c, k_d, k_n, mu_k, mu_d, mu = staalbeton.calculate_parameters()
    bool_point_in_polygon = staalbeton.check_point_within_polygon()
    print(bool_point_in_polygon)

    fig = staalbeton.plot_diagrams('testplots/interaction_stability_diagrams.png')
    

    print(mu)
    print(staalbeton.bending_moment_check())

    
