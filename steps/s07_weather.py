"""Step 7 (Table S1): rainfall and soil moisture (Open-Meteo historical archive; downloaded once to data/raw/weather/, retrieval date in config).
- rain_default: daily rainfall, archive default (ECMWF IFS analysis, 9 km, from 2017) at the rain point;
- rain_era5: daily rainfall, ERA5 model, same point (Section 5.1);
- soil_era5land: hourly 28-100 cm soil moisture, ERA5-Land, at the soil point.
Outputs: outputs/tables/weather.json (period means, water years, ERA5 periods), outputs/plots/fig6a_weather.xlsx (monthly anomalies, 30-day soil moisture)."""
import json, os, urllib.request, numpy as np, pandas as pd

URL = "https://archive-api.open-meteo.com/v1/archive"


def get(path, lat, lon, start, end, var, kind="daily", model=None):
    if not os.path.exists(path):
        q = f"{URL}?latitude={lat}&longitude={lon}&start_date={start}&end_date={end}&{kind}={var}&timezone=Pacific%2FAuckland" + (f"&models={model}" if model else "")
        urllib.request.urlretrieve(q, path)
    d = json.load(open(path))[kind]; return pd.DataFrame(d).assign(time=lambda x: pd.to_datetime(x.time)).set_index("time")


def run(cfg):
    w = cfg["weather"]; raw, out = f"{cfg['paths']['raw']}/weather", cfg["paths"]["out"]; os.makedirs(raw, exist_ok=True); os.makedirs(f"{out}/plots", exist_ok=True)
    rain = get(f"{raw}/rain_default.json", *w["rain_point"], "2017-01-01", "2026-09-30", "precipitation_sum")["precipitation_sum"]
    era5 = get(f"{raw}/rain_era5.json", *w["rain_point"], "2010-01-01", "2026-10-31", "precipitation_sum", model="era5")["precipitation_sum"]
    sm = get(f"{raw}/soil_era5land.json", *w["soil_point"], "2017-01-01", "2026-09-30", w["soil_layer"], kind="hourly", model="era5_land")[w["soil_layer"]].resample("D").mean()
    per = lambda s, a, b: s[(s.index >= a) & (s.index < b)]
    rate = lambda s, a, b: round(float(per(s, a, b).sum() / (len(per(s, a, b)) / 365.25)))
    R = {"periods (IFS rain, mm/yr; soil moisture, m3/m3)": {}}
    for k, (a, b) in {"Jul 2018 - Jun 2022": ("2018-07-01", "2022-07-01"), "Jul 2022 - Jun 2026": ("2022-07-01", "2026-07-01")}.items():
        R["periods (IFS rain, mm/yr; soil moisture, m3/m3)"][k] = dict(rain=rate(rain, a, b), soil_moisture=round(float(per(sm, a, b).mean()), 3))
    y = rain[(rain.index >= "2018-07-01") & (rain.index < "2026-07-01")]; wy = np.where(y.index.month >= 7, y.index.year, y.index.year - 1)
    R["water years (Jul-Jun), IFS rain mm"] = {f"{k}-{str(k + 1)[2:]}": round(float(v)) for k, v in y.groupby(wy).sum().items()}
    dec = lambda t: pd.Timestamp(int(t), 1, 1) + pd.Timedelta(days=(t - int(t)) * 365.25); c0 = cfg["copernicus"]["date"]
    R["ERA5 rain, mm/yr"] = {"LiDAR - Copernicus interval (2013.0 - 2024.24)": rate(era5, dec(c0), dec(cfg["lidar_dates"]["waikato_2024"])),
                             "ICESat-2 record (Oct 2018 - Oct 2026)": rate(era5, "2018-10-01", "2026-11-01"),
                             "Oct 2018 - Jun 2022": rate(era5, "2018-10-01", "2022-07-01"), "Jul 2022 - Oct 2026": rate(era5, "2022-07-01", "2026-11-01")}
    m = rain.resample("MS").sum(); m = m[m.index >= "2017-01-01"]; b0, b1 = w["anomaly_base"]; base = m[(m.index.year >= b0) & (m.index.year <= b1)]
    anom = m - m.index.month.map(base.groupby(base.index.month).mean()).values
    with pd.ExcelWriter(f"{out}/plots/fig6a_weather.xlsx") as xw:
        pd.DataFrame(dict(rain_mm=m, anomaly_mm=anom)).to_excel(xw, sheet_name="monthly rainfall")
        sm.rolling(w["smooth_days"], center=True, min_periods=15).mean().rename("soil_moisture_30d").to_frame().to_excel(xw, sheet_name="soil moisture")
    json.dump(R, open(f"{out}/tables/weather.json", "w"), indent=1); print(json.dumps(R, indent=1), flush=True)
