"""Builds outputs/waikato_is2.qgz: the pipeline's maps in the styles of the paper's figures (rate classes and colours of Figures 3-5,
grey hillshade, S-map Organic Soils outlines). Run with the QGIS Python, e.g. C:\\OSGeo4W\\bin\\python-qgis-ltr.bat qgis/make_project.py <outputs dir> <organic.shp>"""
import os, sys
from qgis.core import (QgsApplication, QgsProject, QgsVectorLayer, QgsRasterLayer, QgsFillSymbol, QgsSingleBandGrayRenderer, QgsContrastEnhancement,
                       QgsColorRampShader, QgsRasterShader, QgsSingleBandPseudoColorRenderer, QgsGraduatedSymbolRenderer, QgsRendererRange,
                       QgsCoordinateReferenceSystem, QgsCategorizedSymbolRenderer, QgsRendererCategory)
from qgis.PyQt.QtGui import QColor

OUT, ORG = sys.argv[1], sys.argv[2]; M = f"{OUT}/maps"
BREAKS = [(-1e9, -40, "#08306b", "≤ −40"), (-40, -30, "#1f5aa6", "−40 to −30"), (-30, -20, "#2a78d6", "−30 to −20"), (-20, -10, "#6fa6e8", "−20 to −10"),
          (-10, -5, "#b5d2f4", "−10 to −5"), (-5, 5, "#efeee9", "−5 to +5"), (5, 10, "#f5c4c3", "+5 to +10"), (10, 20, "#e8807f", "+10 to +20"), (20, 1e9, "#b22b2b", "> +20")]
app = QgsApplication([], False); app.initQgis(); P = QgsProject.instance(); P.setCrs(QgsCoordinateReferenceSystem("EPSG:2193")); root = P.layerTreeRoot()


def add(L, group):
    assert L.isValid(), L.source(); P.addMapLayer(L, False); g = root.findGroup(group) or root.addGroup(group); g.addLayer(L); return L


def rate_raster(name, path, group):
    L = QgsRasterLayer(path, name); fcn = QgsColorRampShader(); fcn.setColorRampType(QgsColorRampShader.Discrete)
    fcn.setColorRampItemList([QgsColorRampShader.ColorRampItem(hi if hi < 1e8 else float("inf"), QColor(c), lab) for lo, hi, c, lab in BREAKS])
    sh = QgsRasterShader(); sh.setRasterShaderFunction(fcn); L.setRenderer(QgsSingleBandPseudoColorRenderer(L.dataProvider(), 1, sh)); return add(L, group)


def rate_vector(name, path, field, group):
    L = QgsVectorLayer(path, name, "ogr"); L.setSubsetString(f'"{field}" IS NOT NULL')
    L.setRenderer(QgsGraduatedSymbolRenderer(field, [QgsRendererRange(lo, hi, QgsFillSymbol.createSimple({"color": c, "outline_color": "white", "outline_width": "0.15"}), lab) for lo, hi, c, lab in BREAKS]))
    return add(L, group)


o = add(QgsVectorLayer(ORG, "Organic Soils (S-map)", "ogr"), "Base")
o.renderer().setSymbol(QgsFillSymbol.createSimple({"style": "no", "outline_color": "#000000", "outline_width": "0.3"}))
for nm, f in (("Distance bands inside Organic Soils", "distance_bands.gpkg"), ("InSAR zones", "insar_zones.gpkg")):
    L = add(QgsVectorLayer(f"{M}/{f}", nm, "ogr"), "Maximum-rate subsets (Figure 5)"); L.renderer().setSymbol(QgsFillSymbol.createSimple({"style": "no", "outline_color": "#e0b800", "outline_width": "0.4"}))
for nm, f in (("a LiDAR − Copernicus, rel. mineral-soil farmland (mm/yr)", "lidcop_rel_240m.tif"), ("b Repeat LiDAR, rel. mineral-soil farmland (mm/yr)", "repeat_rel_240m.tif"),
              ("c InSAR, rel. mineral-soil farmland (mm/yr)", "insar_rel_240m.tif")):
    rate_raster(nm, f"{M}/{f}", "Rate maps, 240 m (Figure 3)")
for nm, fld in (("ICESat-2 (mm/yr)", "rate"), ("LiDAR − Copernicus (mm/yr)", "dem_rate"), ("InSAR (mm/yr)", "insar_rate")):
    rate_vector(nm, f"{M}/hex_3km.gpkg", fld, "Hexagons, 3 km (Figure 4)")
hs = QgsRasterLayer(f"{M}/hillshade_cop30.tif", "Hillshade (Copernicus GLO-30)"); r = QgsSingleBandGrayRenderer(hs.dataProvider(), 1)
ce = QgsContrastEnhancement(hs.dataProvider().dataType(1)); ce.setContrastEnhancementAlgorithm(QgsContrastEnhancement.StretchToMinimumMaximum); ce.setMinimumValue(0); ce.setMaximumValue(255)
r.setContrastEnhancement(ce); hs.setRenderer(r); hs.setOpacity(0.55); add(hs, "Base")
assert P.write(f"{OUT}/waikato_is2.qgz"); print("wrote", f"{OUT}/waikato_is2.qgz")
