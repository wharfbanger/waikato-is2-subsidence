"""Step 9 (Table S1): final outputs.
- outputs/tables/Table3.xlsx, Table4.xlsx and paper_values.json (every number the paper quotes, keyed by name);
- outputs/plots/*.xlsx: the data of each plot (Figure 4d hexagon scatter, Figure 5 maximum rates; Figure 6 written by steps 3 and 7);
- outputs/maps/*.tif|gpkg: 240 m rate maps (steps 4-5), hexagons (step 6), hillshade, distance bands and InSAR zones;
- outputs/waikato_is2.qgz: QGIS project with these layers in the styles of the paper's figures (built with PyQGIS, qgis/make_project.py)."""
import json, os, subprocess, numpy as np, pandas as pd, geopandas as gpd, rasterio
from lib.maps import read

J = lambda out, f: json.load(open(f"{out}/tables/{f}"))


def hillshade(cfg, out):
    dem, tf = read(f"{cfg['paths']['work']}/grid30/cop30_nzvd2016.tif"); gy, gx = np.gradient(np.nan_to_num(dem), 30.0)
    slope = np.arctan(np.hypot(gx, gy)); aspect = np.arctan2(-gx, gy); az, alt = np.radians(315), np.radians(45)
    hs = np.clip(255 * (np.sin(alt) * np.cos(slope) + np.cos(alt) * np.sin(slope) * np.cos(az - aspect)), 0, 255).astype(np.uint8)
    with rasterio.open(f"{out}/maps/hillshade_cop30.tif", "w", driver="GTiff", height=hs.shape[0], width=hs.shape[1], count=1, dtype="uint8", crs="EPSG:2193", transform=tf, compress="deflate") as d: d.write(hs, 1)


def run(cfg):
    out, work = cfg["paths"]["out"], cfg["paths"]["work"]; os.makedirs(f"{out}/plots", exist_ok=True)
    I3, LC, RL, IN, HX, MX, WE, SI, RV, SE = (J(out, f) for f in ("icesat2_rates.json", "lidcop.json", "repeat_lidar.json", "insar.json", "hexagons.json", "max_rates.json", "weather.json", "south_island.json", "repeat_vs_insar.json", "sensitivity.json"))
    rel, ab = I3["relative to mineral-soil farmland"], I3["absolute"]; lc = next(iter(LC.values())); tc = RL["tie"]["classes"]; ic = IN["classes (240 m map at segment locations)"]
    se = lambda d: f"{d['rate']:+.1f} ± {d['se']:.1f}".replace("-", "−")
    T3 = pd.DataFrame([["Organic-Soil farmland", se(rel["Organic-Soil farmland"]), ic["Organic-Soil farmland"], lc["median_rel"], tc["Organic-Soil farmland"]["median"]],
                       ["Mineral-soil farmland", "reference", ic["mineral-soil farmland"], "reference", tc["mineral-soil farmland"]["median"]],
                       ["Hauraki Plains mineral-soil farmland", "–", ic["Hauraki Plains mineral-soil farmland"], "–", tc["Hauraki Plains mineral-soil farmland"]["median"]]],
                      columns=["Land class", "ICESat-2", "InSAR", "LiDAR − Copernicus", "Repeat LiDAR, tied"])
    rows = []
    for k, v in MX.items():
        if not isinstance(v, dict) or "km2" not in v: continue
        f_ = lambda d: f"{d['rate']:+.1f} ± {d['se']:.1f} ({d['passes']})".replace("-", "−") if d and d.get("rate") is not None else "–"
        rows.append([("Distance inside Organic Soils (m)" if k.startswith("band") else "Farmland zoned by InSAR (mm/yr)"), k.split("_", 1)[1].replace("_", "–").replace("le–", "≤ −"),
                     v["km2"], f_(v["Oct 2018 - Jun 2022"]), f_(v["Jul 2022 - Oct 2026"]), f_(v["2018-2026"]), v.get("InSAR map median", "–"), v.get("LiDAR - Copernicus median", "–")])
    T4 = pd.DataFrame(rows, columns=["Subset", "Selection", "Area (km²)", "ICESat-2 Oct 2018 – Jun 2022", "ICESat-2 Jul 2022 – Oct 2026", "ICESat-2 2018–2026", "InSAR map median", "LiDAR − Copernicus median"])
    with pd.ExcelWriter(f"{out}/tables/Table3.xlsx") as xw: T3.to_excel(xw, index=False)
    with pd.ExcelWriter(f"{out}/tables/Table4.xlsx") as xw: T4.to_excel(xw, index=False)
    V = {"icesat2": I3, "lidar_copernicus": LC, "repeat_lidar": RL, "insar": IN, "hexagons": HX, "max_rates": MX, "weather": WE, "south_island": SI, "repeat_vs_insar": RV, "sensitivity": SE}
    json.dump(V, open(f"{out}/tables/paper_values.json", "w"), indent=1)
    hx = gpd.read_file(f"{out}/maps/hex_3km.gpkg")
    with pd.ExcelWriter(f"{out}/plots/fig4d_hexagons.xlsx") as xw: pd.DataFrame(hx.drop(columns="geometry")).to_excel(xw, index=False)
    with pd.ExcelWriter(f"{out}/plots/fig5_max_rates.xlsx") as xw: T4.to_excel(xw, index=False)
    hillshade(cfg, out)
    for nm in ("band", "insar_le"):   # target polygons of step 6 as one layer each
        fs = sorted(f for f in os.listdir(f"{work}/is2rate/max/polys") if f.startswith(nm))
        gpd.GeoDataFrame(pd.concat([gpd.read_file(f"{work}/is2rate/max/polys/{f}").assign(name=f[:-5]) for f in fs], ignore_index=True), crs=2193).to_file(f"{out}/maps/{'distance_bands' if nm == 'band' else 'insar_zones'}.gpkg")
    q = cfg.get("qgis", {})
    if q.get("python"):
        r = subprocess.run([q["python"], os.path.abspath("qgis/make_project.py"), os.path.abspath(out), os.path.abspath(f"{cfg['paths']['raw']}/soils/smap_organic_waikato.shp")], capture_output=True, text=True)
        print(r.stdout[-1500:], r.stderr[-1500:], flush=True)
    print("outputs written", flush=True)
