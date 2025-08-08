import numpy as np
import matplotlib.pyplot as plt

def bereken_interactie_diagram(b, h, f_cd, A_s, f_yd, d_s, aantal_punten=100):
    """
    Bereken het interactiediagram voor een rechthoekige doorsnede van een staal-beton profiel.

    Parameters:
        b (float): Breedte van de betondoorsnede (mm).
        h (float): Hoogte van de betondoorsnede (mm).
        f_cd (float): Betondruksterkte (N/mm^2).
        A_s (float): Doorsnede van de wapening (mm^2).
        f_yd (float): Rekengrenswaarde staal (N/mm^2).
        d_s (float): Effectieve hoogte van de wapening (mm).
        aantal_punten (int): Aantal punten in het diagram.

    Returns:
        momenten (list): Lijst van momenten (kNm).
        krachten (list): Lijst van krachten (kN).
    """
    momenten = []
    krachten = []

    # Itereer over verschillende neutrale lijnhoogtes (h_n)
    for h_n in np.linspace(0.01, h, aantal_punten):
        # Drukkracht in beton
        beta = 0.42  # Locatie van de resulterende kracht voor een parabolisch-rechthoekige verdeling
        if h_n <= h:
            C_c = 0.85 * f_cd * b * h_n  # Resultante drukkracht in beton (N)
        else:
            C_c = 0  # Buiten de neutrale zone draagt beton geen kracht bij
        z_c = beta * h_n  # Hefboomarm beton (mm)

        # Trekkracht in staal
        if h_n < d_s:
            T_s = A_s * f_yd  # Trekkracht in staal (N)
        else:
            T_s = 0  # Als staal volledig in druk is, geen trekbijdrage
        z_s = d_s - h_n / 2  # Hefboomarm staal (mm)

        # Controle op evenwicht van krachten
        if T_s <= C_c:
            # Reductie in betonkracht nodig
            C_c = T_s

        # Momentberekening
        M = (C_c * z_c + T_s * z_s) / 1e6  # Moment (kNm)

        # Opslaan krachten en momenten
        momenten.append(M)
        krachten.append((C_c - T_s) / 1e3)  # Netto kracht (kN)

    return momenten, krachten

# Parameters
b = 300  # mm
h = 500  # mm
f_cd = 25  # N/mm^2
A_s = 2000  # mm^2
f_yd = 500  # N/mm^2
d_s = 450  # mm

# Bereken interactiediagram
momenten, krachten = bereken_interactie_diagram(b, h, f_cd, A_s, f_yd, d_s)

# Plotten
plt.figure(figsize=(8, 6))
plt.plot(momenten, krachten, label="Interactie Diagram")
plt.xlabel("Moment (kNm)")
plt.ylabel("Normaalkracht (kN)")
plt.title("Interactie Diagram voor Staal-Beton Profiel")
plt.grid(True)
plt.legend()
plt.show()
