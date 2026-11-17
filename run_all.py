"""Master script: runs the processing steps of Harvey (2026), 'ICESat-2 confirms peat subsidence in the Waikato, New Zealand', in order.
Each step reads config.yaml and the outputs of earlier steps, and writes to the output folder. Usage (as LiCSBAS batch_LiCSBAS.sh):
    python run_all.py                  # all steps
    python run_all.py --start 3 --end 4
"""
import argparse, importlib, os, time, yaml

STEPS = {1: "s01_prepare", 2: "s02_icesat2", 3: "s03_icesat2_rates", 4: "s04_dem_rates", 5: "s05_insar",
         6: "s06_compare", 7: "s07_weather", 8: "s08_south_island", 9: "s09_outputs"}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--start", type=int, default=1); ap.add_argument("--end", type=int, default=9)
    ap.add_argument("--config", default="config.yaml"); a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config, encoding="utf-8"))
    def merge(d, u):
        for k, v in u.items(): d[k] = merge(d.get(k, {}), v) if isinstance(v, dict) else v
        return d
    if os.path.exists("config.local.yaml"): merge(cfg, yaml.safe_load(open("config.local.yaml", encoding="utf-8")))
    for n in range(a.start, a.end + 1):
        t = time.time(); print(f"--- step {n}: {STEPS[n]}", flush=True)
        importlib.import_module(f"steps.{STEPS[n]}").run(cfg)
        print(f"--- step {n} done in {time.time() - t:.0f} s", flush=True)
