import viktor as vkt
from viktor.parametrization import ViktorParametrization
from viktor.views import ImageResult, ImageView, WebResult, WebView

from io import StringIO
from buckling_calculator import BucklingCalculator
from pathlib import Path
from buiseigenschappen import BuisEigenschappen


class Parametrization(vkt.Parametrization):
    welcome = vkt.Text("# EC4 Buckling calculation \n## Flexural buckling analysis according to EN-1993-1-1 art. 6.3.1.2\n Define your input here:")
    diameter = vkt.NumberField('Diameter [mm]', default=559)
    webthickness = vkt.NumberField('Moment of Inertia [mm^4]', default=8)
    l_buc = vkt.NumberField('Buckling length [mm]', default=10000)
    fcd = vkt.NumberField('f_cd [N/mm^2]', default=17.7)
    fy = vkt.NumberField('Yield strength [MPa]', default=355)
    selected_curve = vkt.OptionField('Buckling curve', options=["curve_a", "curve_b", "curve_c", "curve_d", "curve_e"], default="curve_b")
    pass # Welcome to VIKTOR! You can add your input fields here. Happy Coding!


class Controller(vkt.Controller):
    parametrization = Parametrization
    @ImageView("Buckling calculation", duration_guess=1)

    def createPlot(self, params, **kwargs):
        buis = BuisEigenschappen(t=params.webthickness, diameter=params.diameter, fyd=params.fy)
                
        l_buc = params.l_buc
        fy = params.fy
        selected_curve = params.selected_curve
        N_pl_rk = 0
        calculator = BucklingCalculator(buis.EIeff(), 1, l_buc, buis.NplRd, fy, selected_curve, True, N_pl_rk=buis.NplRd)
        figure = calculator.main()
        redfactor = calculator.reduction_factor
        figure.show()
        svg_data = StringIO()
        figure.savefig(svg_data, format='svg')
        return ImageResult(svg_data)
    @vkt.ImageView("Table 6.2")
    def create_img_result(self, params, **kwargs):
        image_path = Path(__file__).parent / 'Table_6_2.png'
        return vkt.ImageResult.from_path(image_path)

#viktor-cli publish --registered-name ec4-buckling --tag v0.1.0