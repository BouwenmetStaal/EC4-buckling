import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import quad
from dataclasses import dataclass

#https://engcourses-uofa.ca/books/statics/moments-of-inertia-of-area/rectangular-moment-of-inertia/
#https://en.wikipedia.org/wiki/Section_modulus

@dataclass
class TruncatedCircle:
    r: float
    y_flat: float

    # Functie om de cirkel te plotten
    def circle_with_flat_top(self):        
        r = self.r
        y_flat = self.y_flat
        # Hoeken voor de cirkel (boven- en onderzijde)
        theta = np.linspace(0, 2 * np.pi, 1000)
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        
        # Selecteer alleen punten binnen het vlakke gebied
        mask = (y <= y_flat) & (y >= -y_flat)
        x_circle = x[mask]
        y_circle = y[mask]

        x_combined = x_circle
        y_combined = y_circle

        return x_combined, y_combined


    def moment_of_inertia(self):
        """
        Calculate the moment of inertia for a strip on a semicircle.

        Parameters:
            R (float): Radius of the semicircle.
            y_max (float): Upper limit for the strip.

        Returns:
            float: Moment of inertia about the x-axis.
        """
        R = self.r
        y_max = self.y_flat
        # Define the integrand
        def integrand(theta):
            return (R**4 / 4) * (1 - np.cos(4 * theta))

        # Calculate the bounds in terms of theta
        theta_max = np.arcsin(y_max / R)

        # Perform the integration
        I_x, _ = quad(integrand, -theta_max, theta_max)
        
        return I_x
    
    def oppervlakte_strook(self, debug=False):
        """
        Bereken de oppervlakte van een horizontale strook binnen een cirkel.

        Parameters:
        R (float): Straal van de cirkel.
        h (float): Hoogte van de strook vanaf de x-as.
        debug (bool): Indien True, worden tussenresultaten geprint.

        Returns:
        float: Oppervlakte van de strook.
        """
        R = self.r
        h = self.y_flat
        if not (0 <= h <= R):
            A_half = 0.5 * math.pi * R**2
            return A_half
        if R <= 0:
            raise ValueError("De straal moet groter zijn dan 0.")

        # Bereken halve oppervlakte van de cirkel
        A_half = 0.5 * math.pi * R**2
        
        # Bereken de hoek in radialen
        theta = math.pi - 2 * math.asin(h / R)

        # Oppervlakte van de sector
        A_sector = 0.5 * R**2 * theta

        # Oppervlakte van de driehoek
        A_triangle = 0.5 * R**2 * math.sin(theta)

        # Oppervlakte van het cirkelsegment
        A_segment = A_sector - A_triangle

        # Oppervlakte van de strook
        A_strook = A_half - A_segment
        A_stroken = 2*A_strook

        if debug:
            print(f"Halve oppervlakte van de cirkel (A_half): {A_half:.2f} cm²")
            print(f"Hoek in radialen (theta): {theta:.4f} rad")
            print(f"Oppervlakte van de sector (A_sector): {A_sector:.2f} cm²")
            print(f"Oppervlakte van de driehoek (A_triangle): {A_triangle:.2f} cm²")
            print(f"Oppervlakte van het cirkelsegment (A_segment): {A_segment:.2f} cm²")
            print(f"Oppervlakte van de strook (A_strook): {A_strook:.2f} cm²")

        return A_stroken

    def area_of_strip(self):
        """
        Calculate the area of the strip on a semicircle.

        Parameters:
            R (float): Radius of the semicircle.
            y_max (float): Upper limit for the strip.

        Returns:
            float: Area of the strip.
        """
        R = self.r
        y_max = self.y_flat
        # Define the integrand for area
        def area_integrand(y):
            return 2 * np.sqrt(R**2 - y**2)

        # Perform the integration
        A, _ = quad(area_integrand, -y_max, y_max)
        
        return A
    
    def area_of_assymmetric_strip(self, min_y, max_y):
        """
        Calculate the area of the strip on a semicircle.

        Parameters:
            R (float): Radius of the semicircle.
            y_max (float): Upper limit for the strip.

        Returns:
            float: Area of the strip.
        """
        R = self.r
        y_max = self.y_flat
        # Define the integrand for area
        def area_integrand(y):
            return 2 * np.sqrt(R**2 - y**2)

        # Perform the integration
        if min_y < -R:
            min_y = -R
        if max_y > R:
            max_y = R

        A, _ = quad(area_integrand, min_y, max_y)
        if np.isnan(A):
            A = 0
        
        return A


    def elastic_section_modulus(self):
        """
        Calculate the section modulus for the strip on a semicircle.

        Parameters:
            R (float): Radius of the semicircle.
            y_max (float): Upper limit for the strip.

        Returns:
            float: Section modulus about the x-axis.
        """
        y_max = self.y_flat
        # Moment of inertia
        I_x = self.moment_of_inertia()
        
        # Maximum distance from the x-axis (y_max)
        Z = I_x / y_max

        return Z
    
    def centroid_of_strip(self):
        """
        Calculate the centroid of the strip on a semicircle.

        Parameters:
            R (float): Radius of the semicircle.
            y_max (float): Upper limit for the strip.

        Returns:
            float: Distance from the x-axis to the centroid of the strip.
        """
        R = self.r
        y_max = self.y_flat
        # Define the integrand for the first moment of area
        def centroid_integrand(y):
            return 2 * np.sqrt(R**2 - y**2) * y

        # Perform the integration for the first moment of area
        first_moment, _ = quad(centroid_integrand, 0, y_max)

        # Calculate the area of the strip
        A = 0.5*self.area_of_strip()

        # Centroid is the first moment divided by the area
        if first_moment != 0:
            y_centroid = first_moment / A
            return y_centroid
        else:
            return 0
        
    def centroid_of_assym_strip(self, min_y, max_y):
        """
        Calculate the centroid of the strip on a semicircle.

        Parameters:
            R (float): Radius of the semicircle.
            y_max (float): Upper limit for the strip.

        Returns:
            float: Distance from the x-axis to the centroid of the strip.
        """
        R = self.r
        # Define the integrand for the first moment of area
        if min_y < -R:
            min_y = -R
        if max_y > R:
            max_y = R
        def centroid_integrand(y):
            return 2 * np.sqrt(R**2 - y**2) * y

        # Perform the integration for the first moment of area
        first_moment, _ = quad(centroid_integrand, min_y, max_y)

        # Calculate the area of the strip
        A = self.area_of_assymmetric_strip(min_y, max_y)

        # Centroid is the first moment divided by the area
        if first_moment != 0 and A != 0:
            y_centroid = first_moment / A
            return y_centroid
        else:
            return 0
    
    def plastic_section_modulus(self):
        """
        Calculate the plastic section modulus for the strip on a semicircle.

        Parameters:
            R (float): Radius of the semicircle.
            y_max (float): Upper limit for the strip.

        Returns:
            float: Plastic section modulus about the x-axis.
        """
        R = self.r
        y_max = self.y_flat
        # Calculate the area of the strip
        A = self.area_of_strip()

        # Calculate the distance from the x-axis to the centroid
        y_centroid = self.centroid_of_strip()

        # Plastic section modulus is area times the centroid distance
        Z_plastic = A * y_centroid

        return Z_plastic
    
    def plastic_section_modulus_assymetric(self, min_y, max_y):
        """
        Calculate the plastic section modulus for the strip between min_y and max_y on the circle.

        Returns:
            float: Plastic section modulus about the x-axis.
        """
        # Calculate the area of the strip
        A = self.area_of_assymmetric_strip(min_y, max_y)
        y_centroid = self.centroid_of_assym_strip(min_y, max_y)
        Z_plastic = A * y_centroid
        print(f"Z_plastic: {Z_plastic:.2f} mm³ centroid: {y_centroid:.2f} mm A: {A:.2f} mm²")
        return Z_plastic

@dataclass   
class Circle:
    radius: float
        
    def area(self):
        """
        Calculate the area of a full circle.

        Returns:
            float: Area of the full circle.
        """
        R = self.radius
        A_full = np.pi * R**2
        return A_full
    
    def moment_of_inertia(self):
        """
        Calculate the moment of inertia for a full circle.

        Returns:
            float: Moment of inertia about the x-axis.
        """
        R = self.radius
        I_x_full = (np.pi * R**4) / 4
        return I_x_full

    def elastic_section_modulus(self):
        """
        Calculate the section modulus for a full circle.

        Returns:
            float: Section modulus about the x-axis.
        """
        R = self.radius
        # Moment of inertia for a full circle
        I_x_full = (np.pi * R**4) / 4
            
        # Maximum distance from the x-axis (R)
        Z_full = I_x_full / R

        return Z_full

    def plastic_section_modulus(self):
        """
        Calculate the plastic section modulus for a full circle.

        Parameters:

        Returns:
            float: Plastic section modulus about the x-axis.
        """
        R = self.radius
        # The plastic section modulus for a full circle is given by:
        Z_plastic = (4 * R**3) / 3

        return Z_plastic

if __name__ == "__main__":
    # Parameters voor de cirkel
    r = 0.5*559-8  # Straal
    y_flat = 94.33   # Hoogte van de vlakke bovenzijde
    circle = TruncatedCircle(r, y_flat)

    # Genereer de punten
    x, y = circle.circle_with_flat_top()

    # Plot the area of the strip for different values of y_flat
    y_flat_values = np.linspace(0, r, 100)
    areas = [TruncatedCircle(r, y_flat).area_of_strip() for y_flat in y_flat_values]

    # Add the area of the full circle
    full_circle_area = Circle(r).area()

    # Plot the moment of inertia for different values of y_flat
    moments_of_inertia = [TruncatedCircle(r, y_flat).moment_of_inertia() for y_flat in y_flat_values]

    # Add the moment of inertia of the full circle
    full_circle_moment_of_inertia = Circle(r).moment_of_inertia()

    # Plot the section modulus for different values of y_flat
    elastic_section_moduli = [TruncatedCircle(r, y_flat).elastic_section_modulus() for y_flat in y_flat_values]
    # Plot the plastic section modulus for different values of y_flat
    plastic_section_moduli = [TruncatedCircle(r, y_flat).plastic_section_modulus() for y_flat in y_flat_values]


    # Add the section modulus of the full circle
    full_circle_elastic_section_modulus = Circle(r).elastic_section_modulus()
    full_circle_plastic_section_modulus = Circle(r).plastic_section_modulus()

    # Create a figure with a single subplot without axis and data
    fig, ax = plt.subplots(figsize=(20, 20))
    ax.axis('off')

    # Create a grid for the plots
    grid = plt.GridSpec(3, 2, wspace=0.4, hspace=0.4)

    # Plot the truncated circle on the left
    ax0 = plt.subplot(grid[:, 0])
    ax0.plot(x, y, label='Afgekapte Cirkel')
    ax0.fill(x, y, 'b', alpha=0.3)
    ax0.set_aspect('equal', adjustable='box')

    # Plot the full circle in red
    theta_full = np.linspace(0, 2 * np.pi, 1000)
    x_full = r * np.cos(theta_full)
    y_full = r * np.sin(theta_full)
    ax0.plot(x_full, y_full, color='red', linestyle='--', label='Volledige Cirkel')

    ax0.set_title('Afgekapte Cirkel met Vlakke Bovenzijde')
    ax0.set_xlabel('x')
    ax0.set_ylabel('y')
    ax0.legend()
    ax0.grid(True)

    # Plot the area of the strip on the top right
    ax1 = plt.subplot(grid[0, 1])
    ax1.plot(y_flat_values, areas, label='Area of the strip')
    ax1.scatter([r], [full_circle_area], color='red', zorder=5, label='Full Circle Area')
    ax1.scatter([y_flat], [circle.area_of_strip()], color='blue', zorder=5, label='Current Circle Area')
    ax1.set_title(f'Area of the strip for different y_flat values [mm2]')
    ax1.text(y_flat, circle.area_of_strip(), f'{circle.area_of_strip():.2f}', fontsize=9, verticalalignment='bottom', horizontalalignment='right')
    ax1.text(r, full_circle_area, f'{full_circle_area:.2f}', fontsize=9, verticalalignment='bottom', horizontalalignment='right', color='red')
    ax1.set_xlabel('y_flat')
    ax1.set_ylabel('Area')
    ax1.legend()
    ax1.grid(True)

    # Plot the section modulus in the middle right
    ax2 = plt.subplot(grid[1, 1])    
    ax2.plot(y_flat_values, plastic_section_moduli, label='Plastic Section Modulus')
    ax2.plot(y_flat_values, elastic_section_moduli, label='Elastic Section Modulus')
    ax2.scatter([r], [full_circle_elastic_section_modulus], color='red', zorder=5, label='Full Circle Elastic Section Modulus')
    ax2.scatter([r], [full_circle_plastic_section_modulus], color='red', zorder=5, label='Full Circle Plastic Section Modulus')
    ax2.scatter([y_flat], [circle.elastic_section_modulus()], color='blue', zorder=5, label='Current Circle Elastic Section Modulus')
    ax2.scatter([y_flat], [circle.plastic_section_modulus()], color='blue', zorder=5, label='Current Circle Plastic Section Modulus')
    ax2.set_title(f'Section Modulus for different y_flat values [mm3]')
    ax2.text(y_flat, circle.elastic_section_modulus(), f'{circle.elastic_section_modulus():.2f}', fontsize=9, verticalalignment='bottom', horizontalalignment='right')
    ax2.text(y_flat, circle.plastic_section_modulus(), f'{circle.plastic_section_modulus():.2f}', fontsize=9, verticalalignment='bottom', horizontalalignment='right')
    ax2.text(r, full_circle_elastic_section_modulus, f'{full_circle_elastic_section_modulus:.2f}', fontsize=9, verticalalignment='bottom', horizontalalignment='right', color='red')
    ax2.text(r, full_circle_plastic_section_modulus, f'{full_circle_plastic_section_modulus:.2f}', fontsize=9, verticalalignment='bottom', horizontalalignment='right', color='red')
    ax2.set_xlabel('y_flat')
    ax2.set_ylabel('Section Modulus')
    ax2.legend()
    ax2.grid(True)

    # Plot the moment of inertia on the bottom right
    ax3 = plt.subplot(grid[2, 1])
    ax3.plot(y_flat_values, moments_of_inertia, label='Moment of Inertia')
    ax3.scatter([r], [full_circle_moment_of_inertia], color='red', zorder=5, label='Full Circle Moment of Inertia')
    ax3.scatter([y_flat], [circle.moment_of_inertia()], color='blue', zorder=5, label='Current Circle Moment of Inertia')
    ax3.set_title(f'Moment of Inertia for different y_flat values [mm4]')
    ax3.text(y_flat, circle.moment_of_inertia(), f'{circle.moment_of_inertia():.2f}', fontsize=9, verticalalignment='bottom', horizontalalignment='right')
    ax3.text(r, full_circle_moment_of_inertia, f'{full_circle_moment_of_inertia:.2f}', fontsize=9, verticalalignment='bottom', horizontalalignment='right', color='red')
    ax3.set_xlabel('y_flat')
    ax3.set_ylabel('Moment of Inertia')
    ax3.legend()
    ax3.grid(True)

    # Adjust layout
    plt.show()
    # Save the figure as an image
    fig.savefig(f"truncated_circle_analysis_2r_{2*r}_yflat_{y_flat}.png", dpi=300)
