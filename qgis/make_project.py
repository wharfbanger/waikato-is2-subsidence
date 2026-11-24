"""Builds outputs/waikato_is2.qgz: the pipeline's maps in the styles of the paper's figures (rate classes and colours of Figures 3-4; Figure 5a-b classes for the InSAR zones,
grey hillshade, S-map Organic Soils outlines). Run with the QGIS Python, e.g. C:\\OSGeo4W\\bin\\python-qgis-ltr.bat qgis/make_project.py <outputs dir> <organic.shp>"""
import os, sys
from qgis.core import (QgsApplication, QgsProject, QgsVectorLayer, QgsRasterLayer, QgsFillSymbol, QgsSingleBandGrayRenderer, QgsContrastEnhancement,
                       QgsColorRampShader, QgsRasterShader, QgsSingleBandPseudoColorRenderer, QgsGraduatedSymbolRenderer, QgsRendererRange,
                       QgsCoordinateReferenceSystem, QgsCategorizedSymbolRenderer, QgsRendererCategory, QgsFeatureRequest)
from qgis.PyQt.QtGui import QColor

OUT, ORG = sys.argv[1], sys.argv[2]; M = f"{OUT}/maps"
BREAKS = [(-1e9, -40, "#08306b", "≤ −40"), (-40, -30, "#1f5aa6", "−40 to −30"), (-30, -20, "#2a78d6", "−30 to −20"), (-20, -10, "#6fa6e8", "−20 to −10"),
          (-10, -5, "#b5d2f4", "−10 to −5"), (-5, 5, "#efeee9", "−5 to +5"), (5, 10, "#f5c4c3", "+5 to +10"), (10, 20, "#e8807f", "+10 to +20"), (20, 1e9, "#b22b2b", "> +20")]
BRK5 = [(-1e9, -60.05, "#0a2654", "< −60"), (-60.05, -50, "#164380", "−60 to −50"), (-50, -40, "#2560aa", "−50 to −40"), (-40, -30, "#3c7fd0", "−40 to −30"),
        (-30, -20, "#6aa2e3", "−30 to −20"), (-20, -10, "#9cc3ef", "−20 to −10"), (-10, -5, "#c6dbf5", "−10 to −5"),
        (-5, 5, "#efeee9", "−5 to +5"), (5, 10, "#f5c4c3", "+5 to +10"), (10, 20, "#e8807f", "+10 to +20"), (20, 1e9, "#b22b2b", "> +20")]   # Figure 5a-b
app = QgsApplication([], False); app.initQgis(); P = QgsProject.instance(); P.setCrs(QgsCoordinateReferenceSystem("EPSG:2193")); root = P.layerTreeRoot()


def add(L, group):
    assert L.isValid(), L.source(); P.addMapLayer(L, False); g = root.findGroup(group) or root.addGroup(group); g.addLayer(L); return L


def rate_raster(name, path, group):
    L = QgsRasterLayer(path, name); fcn = QgsColorRampShader(); fcn.setColorRampType(QgsColorRampShader.Discrete)
    fcn.setColorRampItemList([QgsColorRampShader.ColorRampItem(hi if hi < 1e8 else float("inf"), QColor(c), lab) for lo, hi, c, lab in BREAKS])
    sh = QgsRasterShader(); sh.setRasterShaderFunction(fcn); L.setRenderer(QgsSingleBandPseudoColorRenderer(L.dataProvider(), 1, sh)); return add(L, group)


def rate_vector(name, path, field, group, breaks=BREAKS, subset=None, outline="white"):
    L = QgsVectorLayer(path, name, "ogr"); L.setSubsetString(f'"{field}" IS NOT NULL' + (f" AND {subset}" if subset else ""))
    L.setRenderer(QgsGraduatedSymbolRenderer(field, [QgsRendererRange(lo, hi, QgsFillSymbol.createSimple({"color": c, "outline_color": outline, "outline_width": "0.15"}), lab) for lo, hi, c, lab in breaks]))
    return add(L, group)


o = add(QgsVectorLayer(ORG, "Organic Soils (S-map)", "ogr"), "Base")
o.renderer().setSymbol(QgsFillSymbol.createSimple({"style": "no", "outline_color": "#000000", "outline_width": "0.3"}))
for nm, fld in (("a ICESat-2 Oct 2018 – Jun 2022 (mm/yr)", "rate_2018_2022"), ("b ICESat-2 Jul 2022 – Oct 2026 (mm/yr)", "rate_2022_2026")):
    # zones at -10 to -60 in steps of 10, strictest drawn last (on top): each area takes the rate of the strictest zone it belongs to
    L = rate_vector(nm, f"{M}/insar_zones.gpkg", fld, "InSAR zones, ICESat-2 rate (Figure 5a-b)", BRK5, '"threshold" % 10 = 0', "0,0,0,0")
    r = L.renderer(); r.setOrderByEnabled(True); r.setOrderBy(QgsFeatureRequest.OrderBy([QgsFeatureRequest.OrderByClause('"threshold"', False)]))
for nm, f in (("a LiDAR − Copernicus, rel. mineral-soil farmland (mm/yr)", "lidcop_rel_240m.tif"), ("b InSAR, rel. mineral-soil farmland (mm/yr)", "insar_rel_240m.tif")):
    rate_raster(nm, f"{M}/{f}", "Rate maps, 240 m (Figure 3)")
for nm, fld in (("a LiDAR − Copernicus (mm/yr)", "dem_rate"), ("b ICESat-2 (mm/yr; also Figure 3c)", "rate"), ("c InSAR (mm/yr)", "insar_rate")):
    rate_vector(nm, f"{M}/hex_3km.gpkg", fld, "Hexagons, 3 km (Figure 4)")
hs = QgsRasterLayer(f"{M}/hillshade_cop30.tif", "Hillshade (Copernicus GLO-30)"); r = QgsSingleBandGrayRenderer(hs.dataProvider(), 1)
ce = QgsContrastEnhancement(hs.dataProvider().dataType(1)); ce.setContrastEnhancementAlgorithm(QgsContrastEnhancement.StretchToMinimumMaximum); ce.setMinimumValue(0); ce.setMaximumValue(255)
r.setContrastEnhancement(ce); hs.setRenderer(r); hs.setOpacity(0.55); add(hs, "Base")
assert P.write(f"{OUT}/waikato_is2.qgz"); print("wrote", f"{OUT}/waikato_is2.qgz")
