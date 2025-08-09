import viktor as vkt
from viktor.parametrization import ViktorParametrization
from viktor.views import ImageResult, ImageView, WebResult, WebView

from io import StringIO
from buckling_calculator import BucklingCalculator
from pathlib import Path
from buiseigenschappen import BuisEigenschappen
from staal_beton_diagram import StaalBetonKolomCalculator


class Parametrization(vkt.Parametrization):
    welcome = vkt.Text("# EC4 Buckling calculation \n## Flexural buckling analysis according to EN-1993-1-1 art. 6.3.1.2\n Define your input here:")
    diameter = vkt.NumberField('Diameter [mm]', default=559)
    webthickness = vkt.NumberField('Moment of Inertia [mm^4]', default=8)
    l_buc = vkt.NumberField('Buckling length [mm]', default=10000)
    fcd = vkt.NumberField('f_cd [N/mm^2]', default=17.7)
    fy = vkt.NumberField('Yield strength [MPa]', default=355)
    My_ed = vkt.NumberField('M,Ed [kNm]', default=500)
    N_ed = vkt.NumberField('N,Ed [kN]', default=1000)

class Controller(vkt.Controller):
    parametrization = Parametrization
    @ImageView("NM and buckling diagram", duration_guess=1)
    def createPlot_1(self, params, **kwargs):
        dt = params.diameter / params.webthickness
        criterion = 90 * (235 / fy)
        if dt < criterion:
            plastic = "Plastic analysis allowed"
        else:
            plastic = "Does not fulfill d/t criterion for plastic analysis"


        buis = BuisEigenschappen(t=params.webthickness, diameter=params.diameter, fyd=params.fy)

        text_bijdrage = buis.staalbijdragefactor

                
        l_buc = params.l_buc
        fy = params.fy
        selected_curve = "curve_b"
        EIeff = buis.EIeff()
        calculator = BucklingCalculator(10**9,
                round(buis.EIeff()),
                l_buc,
                1000,
                round(buis.NplRd(), 2),
                selected_curve,
                True,
                round(1000*buis.NplRk(), 2),)
        Xhi = calculator.reduction_factor
        M_ed = params.My_ed 
        N_ed = params.N_ed
        r = 0
        staalbeton = StaalBetonKolomCalculator(
                M_ed, N_ed, r, Xhi, buis.Nc(), buis.NplRd(), buis.MplRd(), buis.MmaxRd()
            )
        
        fig, _ = staalbeton.return_diagrams()
        fig.show()
        svg_data = StringIO()
        fig.savefig(svg_data, format='svg')
        
        return ImageResult(svg_data)
    @ImageView("Buckling calculation", duration_guess=1)
    def createPlot_2(self, params, **kwargs):
        buis = BuisEigenschappen(t=params.webthickness, diameter=params.diameter, fyd=params.fy)
                
        l_buc = params.l_buc
        fy = params.fy
        selected_curve = "curve_b"
        EIeff = buis.EIeff()
        calculator = BucklingCalculator(10**9,
                round(buis.EIeff()),
                l_buc,
                1000,
                round(buis.NplRd(), 2),
                selected_curve,
                True,
                round(1000*buis.NplRk(), 2),)
        figure = calculator.main()
        redfactor = calculator.reduction_factor
        figure.show()
        svg_data = StringIO()
        figure.savefig(svg_data, format='svg')
        return ImageResult(svg_data)
    @vkt.ImageView("Table 6.3")
    def create_img_result_1(self, params, **kwargs):
        image_path = Path(__file__).parent / 'EC4_tabel6_3.png'
        return vkt.ImageResult.from_path(image_path)
    
    @vkt.ImageView("Table 6.5")
    def create_img_result_2(self, params, **kwargs):
        image_path = Path(__file__).parent / 'EC4_tabel6_5.png'
        return vkt.ImageResult.from_path(image_path)

#viktor-cli publish --registered-name ec4-buckling --tag v0.0.3