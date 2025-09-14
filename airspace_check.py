#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Airspace incursion checker for IGC files using OpenAIP GeoJSON.

Highlights
- Lädt standardmäßig ALLE *.geojson aus ./ressources, IGCs aus ./input, CSVs nach ./output
- Filter: controlled (A–D; E unterdrückt), restricted (ED-R/RESTRICTED), TMZ, RMZ, both/none
- Robuste Datumserkennung aus HFDTE-Headern
- Höhenkorrektur: fester Offset 50 m (H = h − 50), abschaltbar mit --no-geoid-correct
- Schnelle Geometriesuche via STRtree (Shapely)
- Pro Ereignis: CSV-Row + 2 Plots
  • Hauptplot: Airspace-Grenze + Track, OSM-Basemap (sofern verfügbar)
  • Detailplot: Segment vor/während/nach Verletzung (Fokus auf Kernsegment), OSM-Basemap
"""

import argparse
import csv
import datetime as dt
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from shapely.geometry import shape, Point, LineString, Polygon, MultiPolygon
from shapely.geometry.base import BaseGeometry
from shapely.strtree import STRtree
from shapely.ops import transform as shp_transform

import re

# ----------------------------
# Fixed offset altitude "correction"
# ----------------------------

def _apply_fixed_offset(track, offset_m: float = 50.0):
    if not track:
        return track
    new = []
    for (t, lon, lat, alt) in track:
        if alt is None:
            new.append((t, lon, lat, alt))
        else:
            new.append((t, lon, lat, float(alt) - float(offset_m)))
    return new

# ----------------------------
# IGC parsing
# ----------------------------

def parse_igc_datetime(lines: List[str]) -> Optional[dt.datetime]:
    """
    Robustly parse the IGC flight date from H-sentences.
    Supports e.g.:
      HFDTE300423
      HFDTEDATE:300423
      HFDTE:30/04/23
      HFDTE 30-04-2023
    Returns UTC date at 00:00:00.
    """
    # pass 1: canonical 'HFDTE'
    for line in lines:
        u = line.strip().upper()
        if u.startswith("HFDTE"):
            m = re.search(r"HFDTE[^0-9]*([0-9]{6,8})", u)
            if m:
                s = m.group(1)
                try:
                    if len(s) == 6:  # DDMMYY
                        day = int(s[0:2]); month = int(s[2:4]); yy = int(s[4:6])
                        year = 2000 + yy if yy < 80 else 1900 + yy
                    elif len(s) == 8:  # DDMMYYYY
                        day = int(s[0:2]); month = int(s[2:4]); year = int(s[4:8])
                    else:
                        continue
                    return dt.datetime(year, month, day, tzinfo=dt.timezone.utc)
                except Exception:
                    continue
    # pass 2: any H-line with DTE/DATE and DD[sep]MM[sep]YY(YY)
    date_pat = re.compile(r"(\d{1,2})[./\-\s]?(\d{1,2})[./\-\s]?(\d{2,4})")
    for line in lines:
        u = line.strip().upper()
        if u.startswith("H") and ("DTE" in u or "DATE" in u):
            m = date_pat.search(u)
            if m:
                d, mth, y = m.groups()
                try:
                    day = int(d); month = int(mth); year = int(y)
                    if year < 100:
                        year = 2000 + year if year < 80 else 1900 + year
                    return dt.datetime(year, month, day, tzinfo=dt.timezone.utc)
                except Exception:
                    continue
    return None

def dm_to_deg(dm: str, is_lat: bool) -> float:
    if is_lat:
        deg = int(dm[0:2]); mins = float(dm[2:4] + "." + dm[4:7])
    else:
        deg = int(dm[0:3]); mins = float(dm[3:5] + "." + dm[5:8])
    return deg + mins / 60.0

def parse_b_record(line: str, base_date: Optional[dt.datetime]) -> Optional[Tuple[dt.datetime, float, float, Optional[float]]]:
    try:
        if not line.startswith("B") or len(line) < 35:
            return None
        hh = int(line[1:3]); mm = int(line[3:5]); ss = int(line[5:7])
        lat_dm = line[7:14]; lat_hem = line[14]
        lon_dm = line[15:23]; lon_hem = line[23]
        gps_alt = None
        try:
            gps_alt = int(line[30:35])
        except Exception:
            gps_alt = None

        lat = dm_to_deg(lat_dm, True); lat = -lat if lat_hem == 'S' else lat
        lon = dm_to_deg(lon_dm, False); lon = -lon if lon_hem == 'W' else lon

        if base_date is not None:
            t = dt.datetime(base_date.year, base_date.month, base_date.day, hh, mm, ss, tzinfo=dt.timezone.utc)
        else:
            t = dt.datetime(1970, 1, 1, hh, mm, ss, tzinfo=dt.timezone.utc)
        return (t, lon, lat, float(gps_alt) if gps_alt is not None else None)
    except Exception:
        return None

def read_igc_track_with_meta(path: Path):
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        lines = [l.strip() for l in f if l.strip()]
    base_date = parse_igc_datetime(lines)
    track = []
    for line in lines:
        if line.startswith("B"):
            rec = parse_b_record(line, base_date)
            if rec: track.append(rec)
    track.sort(key=lambda x: x[0])
    meta = {}
    return track, meta

# ----------------------------
# Altitude parsing
# ----------------------------

def ft_to_m(x: float) -> float:
    return x * 0.3048

def parse_openaip_alt(limit_obj: Any) -> Tuple[float, bool]:
    if limit_obj is None:
        return (float("inf"), True)
    if isinstance(limit_obj, (int, float)):
        return (float(limit_obj), False)
    if isinstance(limit_obj, str):
        s = limit_obj.strip().upper().replace(" ", "")
        if s in ("GND", "SFC"): return (0.0, False)
        if s in ("UNL", "UNLIM", "UNLIMITED"): return (float("inf"), True)
        if s.startswith("FL"):
            try: return (ft_to_m(float(s[2:]) * 100.0), False)
            except Exception: return (float("inf"), True)
        try:
            if s.endswith("FT"): return (ft_to_m(float(s[:-2])), False)
            if s.endswith("M"): return (float(s[:-1]), False)
            return (float(s), False)
        except Exception:
            return (float("inf"), True)
    if isinstance(limit_obj, dict):
        val = limit_obj.get("value")
        unit = limit_obj.get("unit")
        ref = limit_obj.get("referenceDatum")
        if ref == 0:  # SFC/GND
            return (0.0, False)
        if unit == 6:  # FL
            try: return (ft_to_m(float(val) * 100.0), False)
            except Exception: return (float("inf"), True)
        if unit == 1:  # FT
            try: return (ft_to_m(float(val)), False)
            except Exception: return (float("inf"), True)
        try: return (float(val), False)
        except Exception: return (float("inf"), True)
    return (float("inf"), True)

def extract_alt_bounds(props: Dict[str, Any]) -> Tuple[float, float]:
    lower_obj = props.get("lowerLimit") or props.get("lower") or props.get("floor")
    upper_obj = props.get("upperLimit") or props.get("upper") or props.get("ceiling")
    lower_m, _ = parse_openaip_alt(lower_obj)
    upper_m, _ = parse_openaip_alt(upper_obj)
    return (lower_m, upper_m)

# ----------------------------
# Airspace model
# ----------------------------

class Airspace:
    def __init__(self, name: str, as_type: str, icao_class: Optional[Union[str,int]], geom: BaseGeometry, lower_m: float, upper_m: float):
        self.name = name
        self.as_type = as_type
        self.icao_class = icao_class
        self.geom = geom
        self.lower_m = lower_m
        self.upper_m = upper_m

    def is_inside(self, lon: float, lat: float, alt_m: Optional[float]) -> bool:
        pt = Point(lon, lat)
        inside_2d = self.geom.contains(pt) or self.geom.touches(pt)
        if not inside_2d:
            return False
        if alt_m is None:
            return True
        if math.isfinite(self.lower_m) and alt_m < self.lower_m: return False
        if math.isfinite(self.upper_m) and alt_m > self.upper_m: return False
        return True

def load_airspaces(geojson_path: Path) -> List['Airspace']:
    data = json.loads(geojson_path.read_text(encoding="utf-8"))
    feats = data.get("features", [])
    airspaces: List[Airspace] = []
    for f in feats:
        try:
            geom = shape(f["geometry"])
            if not isinstance(geom, (Polygon, MultiPolygon)):
                continue
            props = f.get("properties", {})
            name = str(props.get("name") or props.get("_id") or "Unbenannt")
            as_type = str(props.get("type", "Unknown"))
            icao_class = props.get("icaoClass")
            lower_m, upper_m = extract_alt_bounds(props)
            airspaces.append(Airspace(name, as_type, icao_class, geom, lower_m, upper_m))
        except Exception:
            continue
    return airspaces

def load_airspaces_many(paths: list) -> List['Airspace']:
    all_spaces: List[Airspace] = []
    for p in paths:
        try:
            all_spaces.extend(load_airspaces(Path(p)))
        except Exception:
            continue
    return all_spaces

# ----------------------------
# Filters
# ----------------------------

def _icao_class_letter(raw) -> str:
    if raw is None:
        return ""
    try:
        s = str(raw).strip().upper()
    except Exception:
        return ""
    if s in {"A","B","C","D","E","F","G"}:
        return s
    try:
        n = int(s)
        mapping = {0:"A",1:"B",2:"C",3:"D",4:"E",5:"F",6:"G"}
        return mapping.get(n, "")
    except Exception:
        return ""

def is_controlled(a: Airspace) -> bool:
    # Controlled = A–D (E ist unterdrückt)
    cls = _icao_class_letter(a.icao_class)
    if cls in {"A","B","C","D"}:
        return True
    name = (a.name or "").upper()
    atype = str(a.as_type).upper()
    controlled_tokens = ["CTR", "CTA", "TMA", "CONTROL ZONE", "CONTROL AREA", "KONTROLLZONE", "KONTROLL"]
    if any(tok in name for tok in controlled_tokens):
        if cls == "E":
            return False
        return True
    if atype in {"CTR", "CTA", "TMA"}:
        if cls == "E":
            return False
        return True
    return False

def is_restricted(a: Airspace) -> bool:
    name = (a.name or "").upper()
    return name.startswith("ED-R") or ("RESTRICTED" in name) or name.startswith("R-")

def is_tmz(a: Airspace) -> bool:
    name = (a.name or "").upper()
    return ("TMZ" in name) or (str(a.as_type).upper()=="TMZ")

def is_rmz(a: Airspace) -> bool:
    name = (a.name or "").upper()
    return ("RMZ" in name) or ("RADIO MANDATORY ZONE" in name) or (str(a.as_type).upper()=="RMZ")

def class_label(a: Airspace) -> str:
    cls = _icao_class_letter(a.icao_class)
    if cls:
        return cls
    name = (a.name or "").upper()
    if "TMZ" in name or str(a.as_type).upper() == "TMZ":
        return "TMZ"
    if "RMZ" in name or "RADIO MANDATORY ZONE" in name or str(a.as_type).upper() == "RMZ":
        return "RMZ"
    if name.startswith("ED-R") or "RESTRICTED" in name or name.startswith("R-"):
        return "R"
    return str(a.as_type).upper()

# ----------------------------
# Geometry helpers
# ----------------------------

def interpolate_time(t0, t1, frac):
    delta=(t1-t0).total_seconds()
    return t0 + dt.timedelta(seconds=delta*max(0.0, min(1.0, frac)))

def segment_fraction_at_border(p0, p1, poly):
    line=LineString([p0,p1])
    inter=line.intersection(poly.boundary)
    if inter.is_empty: return None
    try:
        points=[inter] if inter.geom_type=="Point" else list(inter.geoms)
    except Exception:
        points=[]
    if not points: return None
    def dist(a,b):
        dx=a[0]-b[0]; dy=a[1]-b[1]
        return (dx*dx+dy*dy) ** 0.5
    best=min(points, key=lambda q: dist(p0,(q.x,q.y)))
    total=dist(p0,p1)
    if total==0: return 0.0
    return dist(p0,(best.x,best.y))/total

# ----------------------------
# Plotting
# ----------------------------

def _plot_event(airspace: Airspace, track, entry_time: dt.datetime, exit_time: dt.datetime,
                entry_ll: Tuple[float,float], exit_ll: Tuple[float,float], out_dir: Path,
                basemap: bool=True, basemap_provider: str="CartoDB.Positron"):
    import matplotlib.pyplot as plt

    # Segment ±5 Minuten
    t_start = entry_time - dt.timedelta(minutes=5)
    t_end   = exit_time + dt.timedelta(minutes=5)
    seg_points = [(lon,lat) for (t,lon,lat,alt) in track if t_start <= t <= t_end]

    g = airspace.geom

    if basemap:
        try:
            import contextily as ctx
            import pyproj

            def _resolve_provider(name):
                try:
                    p = ctx.providers
                    for part in str(name).split('.'):
                        p = getattr(p, part)
                    return p
                except Exception:
                    return ctx.providers.OSM

            to_3857 = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True).transform
            g3857 = shp_transform(to_3857, g)
            seg_points_3857 = [to_3857(x,y) for x,y in seg_points] if seg_points else None
            entry_m = to_3857(entry_ll[0], entry_ll[1])
            exit_m  = to_3857(exit_ll[0],  exit_ll[1])

            fig, ax = plt.subplots(figsize=(6,6))
            if g3857.geom_type == "Polygon":
                xs, ys = g3857.exterior.xy; ax.plot(xs, ys, label="Airspace boundary")
            elif g3857.geom_type == "MultiPolygon":
                for poly in g3857.geoms:
                    xs, ys = poly.exterior.xy; ax.plot(xs, ys, label="Airspace boundary")
            if seg_points_3857:
                xs, ys = zip(*seg_points_3857); ax.plot(xs, ys, label="Track")
            ax.plot(entry_m[0], entry_m[1], marker="o", label="Entry")
            ax.plot(exit_m[0],  exit_m[1],  marker="x", label="Exit")
            try:
                bx, by, tx, ty = g3857.bounds
                padx = (tx - bx) * 0.1; pady = (ty - by) * 0.1
                ax.set_xlim(bx - padx, tx + padx); ax.set_ylim(by - pady, ty + pady)
            except Exception:
                pass
            provider = _resolve_provider(basemap_provider)
            ctx.add_basemap(ax, source=provider, crs="EPSG:3857")
            ax.set_xlabel("Web Mercator X"); ax.set_ylabel("Web Mercator Y")
            cls = class_label(airspace)
            ax.set_title(f"{airspace.name} (Class {cls})\n{entry_time.strftime('%Y-%m-%d %H:%M:%S UTC')} → {exit_time.strftime('%Y-%m-%d %H:%M:%S UTC')}")
            ax.legend()
            safe_name = "".join(c if c.isalnum() or c in ("-","_") else "_" for c in airspace.name)[:60]
            fname = f"{entry_time.strftime('%Y%m%dT%H%M%S')}_{safe_name}.png"
            out_path = Path(out_dir) / fname
            fig.savefig(out_path, dpi=150); plt.close(fig)
            return out_path
        except Exception:
            pass

    fig, ax = plt.subplots(figsize=(6,6))
    if g.geom_type == "Polygon":
        xs, ys = g.exterior.xy; ax.plot(xs, ys, label="Airspace boundary")
    elif g.geom_type == "MultiPolygon":
        for poly in g.geoms:
            xs, ys = poly.exterior.xy; ax.plot(xs, ys, label="Airspace boundary")
    if seg_points:
        xs, ys = zip(*seg_points); ax.plot(xs, ys, label="Track")
    ax.plot(entry_ll[0], entry_ll[1], marker="o", label="Entry")
    ax.plot(exit_ll[0],  exit_ll[1],  marker="x", label="Exit")
    ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
    cls = class_label(airspace)
    ax.set_title(f"{airspace.name} (Class {cls})\n{entry_time.strftime('%Y-%m-%d %H:%M:%S UTC')} → {exit_time.strftime('%Y-%m-%d %H:%M:%S UTC')}")
    ax.legend()
    safe_name = "".join(c if c.isalnum() or c in ("-","_") else "_" for c in airspace.name)[:60]
    fname = f"{entry_time.strftime('%Y%m%dT%H%M%S')}_{safe_name}.png"
    out_path = Path(out_dir) / fname
    fig.savefig(out_path, dpi=150); plt.close(fig)
    return out_path

def _plot_event_detail(airspace: Airspace, track, entry_time: dt.datetime, exit_time: dt.datetime,
                       out_dir: Path, basemap: bool=True, basemap_provider: str="CartoDB.Positron"):
    """
    Detailplot der Verletzung:
    - zeigt den Flugweg VOR und NACH der Verletzung (blass, dünn)
    - legt den Fokus auf den Flugweg WÄHREND der Verletzung (dicker)
    - zeigt zusätzlich die Grenze des verletzten Luftraums
    - optionale OSM-Basemap wie in den Hauptplots
    """
    import matplotlib.pyplot as plt

    # Fenster um das Ereignis (±5 Minuten)
    t_pre = entry_time - dt.timedelta(minutes=5)
    t_post = exit_time + dt.timedelta(minutes=5)

    # Segmente aufteilen
    pre = [(lon, lat) for (t, lon, lat, alt) in track if t_pre <= t < entry_time]
    core = [(lon, lat) for (t, lon, lat, alt) in track if entry_time <= t <= exit_time]
    post = [(lon, lat) for (t, lon, lat, alt) in track if exit_time < t <= t_post]
    if len(core) < 2:
        return None

    if basemap:
        try:
            import contextily as ctx
            import pyproj

            def _resolve_provider(name):
                try:
                    p = ctx.providers
                    for part in str(name).split('.'):
                        p = getattr(p, part)
                    return p
                except Exception:
                    return ctx.providers.OSM

            to_3857 = pyproj.Transformer.from_crs("EPSG:4326","EPSG:3857",always_xy=True).transform
            pre3857  = [to_3857(x,y) for x,y in pre]  if pre  else []
            core3857 = [to_3857(x,y) for x,y in core] if core else []
            post3857 = [to_3857(x,y) for x,y in post] if post else []

            # Transform airspace geometry to 3857 and draw its boundary (thin)
            g = airspace.geom
            g3857 = shp_transform(to_3857, g)

            fig, ax = plt.subplots(figsize=(6,6))
            # Airspace boundary overlay (thin, alpha)
            try:
                if g3857.geom_type == "Polygon":
                    xs, ys = g3857.exterior.xy; ax.plot(xs, ys, linewidth=0.8, alpha=0.7, label="Luftraumgrenze")
                elif g3857.geom_type == "MultiPolygon":
                    for poly in g3857.geoms:
                        xs, ys = poly.exterior.xy; ax.plot(xs, ys, linewidth=0.8, alpha=0.7, label="Luftraumgrenze")
            except Exception:
                pass

            if len(pre3857) > 1:
                xs, ys = zip(*pre3857); ax.plot(xs, ys, linewidth=1.0, alpha=0.5, linestyle="--", label="Vor der Verletzung")
            if len(post3857) > 1:
                xs, ys = zip(*post3857); ax.plot(xs, ys, linewidth=1.0, alpha=0.5, linestyle="--", label="Nach der Verletzung")

            xs, ys = zip(*core3857); ax.plot(xs, ys, linewidth=2.5, label="Verletzungssegment")

            allx = [x for pts in (pre3857, core3857, post3857) for x,_ in pts]
            ally = [y for pts in (pre3857, core3857, post3857) for _,y in pts]
            bx, tx = min(allx), max(allx); by, ty = min(ally), max(ally)
            padx = max((tx - bx) * 0.15, 1000); pady = max((ty - by) * 0.15, 1000)
            ax.set_xlim(bx - padx, tx + padx); ax.set_ylim(by - pady, ty + pady)

            provider = _resolve_provider(basemap_provider)
            ctx.add_basemap(ax, source=provider, crs="EPSG:3857")
            ax.set_xlabel("Web Mercator X"); ax.set_ylabel("Web Mercator Y")
            ax.set_title(f"Detail: Luftraumverletzung\n{entry_time.strftime('%Y-%m-%d %H:%M:%S UTC')} → {exit_time.strftime('%Y-%m-%d %H:%M:%S UTC')}")
            ax.legend()

            fname = f"{entry_time.strftime('%Y%m%dT%H%M%S')}_detail.png"
            out_path = Path(out_dir) / fname
            fig.savefig(out_path, dpi=150); plt.close(fig)
            return out_path
        except Exception:
            pass

    fig, ax = plt.subplots(figsize=(6,6))
    # Airspace boundary in lon/lat (thin)
    try:
        g = airspace.geom
        if g.geom_type == "Polygon":
            xs, ys = g.exterior.xy; ax.plot(xs, ys, linewidth=0.8, alpha=0.7, label="Luftraumgrenze")
        elif g.geom_type == "MultiPolygon":
            for poly in g.geoms:
                xs, ys = poly.exterior.xy; ax.plot(xs, ys, linewidth=0.8, alpha=0.7, label="Luftraumgrenze")
    except Exception:
        pass
    if len(pre) > 1:
        xs, ys = zip(*pre); ax.plot(xs, ys, linewidth=1.0, alpha=0.5, linestyle="--", label="Vor der Verletzung")
    if len(post) > 1:
        xs, ys = zip(*post); ax.plot(xs, ys, linewidth=1.0, alpha=0.5, linestyle="--", label="Nach der Verletzung")
    xs, ys = zip(*core); ax.plot(xs, ys, linewidth=2.5, label="Verletzungssegment")
    allx = [x for pts in (pre, core, post) for x,_ in pts]
    ally = [y for pts in (pre, core, post) for _,y in pts]
    bx, tx = min(allx), max(allx); by, ty = min(ally), max(ally)
    padx = max((tx - bx) * 0.15, 0.01); pady = max((ty - by) * 0.15, 0.01)
    ax.set_xlim(bx - padx, tx + padx); ax.set_ylim(by - pady, ty + pady)
    ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
    ax.set_title(f"Detail: Luftraumverletzung\n{entry_time.strftime('%Y-%m-%d %H:%M:%S UTC')} → {exit_time.strftime('%Y-%m-%d %H:%M:%S UTC')}")
    ax.legend()
    fname = f"{entry_time.strftime('%Y%m%dT%H%M%S')}_detail.png"
    out_path = Path(out_dir) / fname
    fig.savefig(out_path, dpi=150); plt.close(fig)
    return out_path

# ----------------------------
# Analyzer
# ----------------------------

def analyze_track_with_index(track, airspaces, plots=False, plots_dir=None, filter_mode="both",
                             basemap=True, basemap_provider="CartoDB.Positron"):
    results = []
    if len(track) < 2: return results
    geoms = [a.geom for a in airspaces]
    tree = STRtree(geoms)
    geom_to_idx = {g.wkb: i for i, g in enumerate(geoms)}

    def pass_filter(a: Airspace) -> bool:
        if filter_mode == "none": return True
        if filter_mode == "controlled": return is_controlled(a)
        if filter_mode == "restricted": return is_restricted(a)
        if filter_mode == "tmzrmz": return is_tmz(a) or is_rmz(a)
        return is_controlled(a) or is_restricted(a) or is_tmz(a) or is_rmz(a)

    if plots:
        try:
            import matplotlib.pyplot as plt  # noqa: F401
        except Exception:
            print("Warnung: matplotlib nicht verfügbar – Plots werden übersprungen.", file=sys.stderr)
            plots = False
        if basemap:
            try:
                import contextily as ctx  # noqa: F401
                import pyproj  # noqa: F401
            except Exception:
                print("Hinweis: contextily/pyproj nicht verfügbar – Basemap wird ignoriert.", file=sys.stderr)
                basemap = False
        if plots and plots_dir:
            Path(plots_dir).mkdir(parents=True, exist_ok=True)

    active: Dict[int, Dict[str, Any]] = {}

    for i in range(1, len(track)):
        t0, lon0, lat0, alt0 = track[i-1]
        t1, lon1, lat1, alt1 = track[i]

        pt1 = Point(lon1, lat1)
        candidates_geoms = tree.query(pt1)

        candidate_idxs = set()
        for g in candidates_geoms:
            if hasattr(g, 'wkb'):
                idx = geom_to_idx.get(g.wkb)
                if idx is None:
                    try:
                        idx = geoms.index(g)
                    except ValueError:
                        idx = None
                if idx is not None:
                    candidate_idxs.add(idx)
            else:
                try:
                    idx_int = int(g)
                    if 0 <= idx_int < len(geoms):
                        candidate_idxs.add(idx_int)
                except Exception:
                    pass

        inside_now = set()
        for idx in candidate_idxs:
            a = airspaces[idx]
            if not pass_filter(a): continue
            if a.is_inside(lon1, lat1, alt1): inside_now.add(idx)

        # Entries
        for idx in inside_now - set(active.keys()):
            a = airspaces[idx]
            frac = segment_fraction_at_border((lon0,lat0),(lon1,lat1), a.geom)
            if frac is None:
                et = t1; el = (lon1, lat1)
            else:
                et = interpolate_time(t0, t1, frac)
                el = (lon0 + (lon1 - lon0) * frac, lat0 + (lat1 - lat0) * frac)
            active[idx] = {"entry": et, "entry_ll": el, "max_alt": (alt1 if alt1 is not None else float("-inf"))}

        # Update max altitude
        for idx in inside_now:
            if idx in active and alt1 is not None:
                prev = active[idx].get("max_alt", float("-inf"))
                if alt1 > prev:
                    active[idx]["max_alt"] = alt1

        # Exits
        to_close = []
        for idx in list(active.keys()):
            a = airspaces[idx]
            if not a.is_inside(lon1, lat1, alt1):
                frac = segment_fraction_at_border((lon0,lat0),(lon1,lat1), a.geom)
                if frac is None:
                    xt = t1; xl = (lon1, lat1)
                else:
                    xt = interpolate_time(t0, t1, frac)
                    xl = (lon0 + (lon1 - lon0) * frac, lat0 + (lat1 - lat0) * frac)
                entry = active[idx]["entry"]
                duration = (xt - entry).total_seconds()
                max_alt_val = active[idx].get("max_alt")
                res = {
                    "airspace_name": a.name,
                    "airspace_type": a.as_type,
                    "class": class_label(a),
                    "icao_class_raw": a.icao_class,
                    "entry": entry,
                    "exit": xt,
                    "duration_s": max(0.0, duration),
                    "entry_lon": active[idx]["entry_ll"][0],
                    "entry_lat": active[idx]["entry_ll"][1],
                    "exit_lon": xl[0],
                    "exit_lat": xl[1],
                    "lower_m": a.lower_m,
                    "upper_m": a.upper_m,
                    "max_alt_m": (max_alt_val if max_alt_val is not None and max_alt_val != float("-inf") else None)
                }
                results.append(res)
                if plots and plots_dir:
                    try:
                        _plot_event(airspace=a, track=track, entry_time=entry, exit_time=xt,
                                    entry_ll=(res["entry_lon"], res["entry_lat"]),
                                    exit_ll=(res["exit_lon"], res["exit_lat"]),
                                    out_dir=plots_dir, basemap=basemap, basemap_provider=basemap_provider)
                        _plot_event_detail(airspace=a, track=track, entry_time=entry, exit_time=xt,
                                           out_dir=plots_dir, basemap=basemap, basemap_provider=basemap_provider)
                    except Exception as e:
                        print(f"Plot-Fehler für {a.name}: {e}", file=sys.stderr)
                to_close.append(idx)
        for idx in to_close:
            active.pop(idx, None)

    # Close any still-active incursions
    if active:
        t_last, lon_last, lat_last, _ = track[-1]
        for idx, info in active.items():
            a = airspaces[idx]
            res = {
                "airspace_name": a.name,
                "airspace_type": a.as_type,
                "class": class_label(a),
                "icao_class_raw": a.icao_class,
                "entry": info["entry"],
                "exit": t_last,
                "duration_s": max(0.0, (t_last - info["entry"]).total_seconds()),
                "entry_lon": info["entry_ll"][0],
                "entry_lat": info["entry_ll"][1],
                "exit_lon": lon_last,
                "exit_lat": lat_last,
                "lower_m": a.lower_m,
                "upper_m": a.upper_m,
                "max_alt_m": (info.get("max_alt") if info.get("max_alt", float("-inf")) != float("-inf") else None)
            }
            results.append(res)
            if plots and plots_dir:
                try:
                    _plot_event(airspace=a, track=track, entry_time=info["entry"], exit_time=t_last,
                                entry_ll=(res["entry_lon"], res["entry_lat"]),
                                exit_ll=(res["exit_lon"], res["exit_lat"]),
                                out_dir=plots_dir, basemap=basemap, basemap_provider=basemap_provider)
                    _plot_event_detail(airspace=a, track=track, entry_time=info["entry"], exit_time=t_last,
                                       out_dir=plots_dir, basemap=basemap, basemap_provider=basemap_provider)
                except Exception as e:
                    print(f"Plot-Fehler für {a.name}: {e}", file=sys.stderr)

    return results

# ----------------------------
# Main
# ----------------------------

def main():
    # Defaults
    default_geojson_dir = Path("ressources")
    default_geojson = default_geojson_dir / "de_asp.geojson"
    default_igc_dir = Path("input")

    ap = argparse.ArgumentParser(description="Prüft IGC-Dateien auf Luftraum-Einflüge; lädt alle GeoJSON aus ./ressources; CSVs nach ./output; OSM-Basemap; Detailplots.")
    ap.add_argument("geojson", nargs="?", type=str, default=str(default_geojson),
                    help=f"Pfad zur GeoJSON-Datei ODER Ordner mit *.geojson. Default: {default_geojson}")
    ap.add_argument("--dir", type=str, default=str(default_igc_dir),
                    help=f"Verzeichnis mit *.igc (Default: {default_igc_dir}).")
    ap.add_argument("--no-alt", action="store_true", help="Höhenlimits ignorieren (nur 2D).")
    ap.add_argument("--filter", type=str, choices=["both","controlled","restricted","tmzrmz","none"], default="both",
                    help="Filter: kontrollierte Lufträume (A–D; E unterdrückt), Sperrgebiete (ED-R/RESTRICTED), TMZ/RMZ, beides oder kein Filter.")
    ap.add_argument("--plots", action="store_true", help="Für jedes Ereignis ein PNG erzeugen.")
    ap.add_argument("--plots-dir", type=str, default="event_plots", help="Zielordner für PNGs (Default: event_plots).")
    ap.add_argument("--no-basemap", dest="basemap", action="store_false", help="Basemap deaktivieren.")
    ap.add_argument("--basemap-provider", type=str, default="CartoDB.Positron", help="Basemap-Provider (z. B. OSM, Stamen.Terrain, CartoDB.Positron).")
    ap.add_argument("--geoid-correct", dest="geoid_correct", action="store_true", default=True,
                    help="Höhen pauschal um 50 m auf MSL korrigieren (Default: an).")
    ap.add_argument("--no-geoid-correct", dest="geoid_correct", action="store_false",
                    help="50-m-Korrektur deaktivieren.")
    args = ap.parse_args()

    # GeoJSONs bestimmen
    geojson_path = Path(args.geojson).resolve()
    geojson_files = []
    if geojson_path.is_dir():
        geojson_files = sorted(geojson_path.glob("*.geojson"))
    else:
        # Wenn Standardpfad genutzt, lade alle *.geojson aus ./ressources
        if geojson_path == default_geojson.resolve() and default_geojson_dir.exists():
            geojson_files = sorted(default_geojson_dir.glob("*.geojson"))
        elif geojson_path.exists() and geojson_path.is_file():
            geojson_files = [geojson_path]

    if not geojson_files:
        print(
            f"GeoJSON nicht gefunden: {geojson_path}\n"
            f"Hinweis: Lege Dateien unter '{default_geojson_dir.resolve()}' ab (z. B. *.geojson) oder gib einen Pfad/Ordner an.",
            file=sys.stderr
        )
        sys.exit(1)

    airspaces = load_airspaces_many(geojson_files)
    print(f"  → GeoJSON-Dateien geladen: {len(geojson_files)} | Lufträume insgesamt: {len(airspaces)}")

    if args.no_alt:
        for a in airspaces:
            a.lower_m = float("-inf"); a.upper_m = float("inf")

    igc_dir = Path(args.dir).resolve()
    candidates = list(igc_dir.glob('*.[iI][gG][cC]'))
    seen = {}
    for f in candidates:
        key = str(f.resolve()).lower()
        seen[key] = f
    igc_files = sorted(seen.values())
    if not igc_files:
        print(f"Keine IGC-Dateien in {igc_dir} gefunden.\n"
              f"Hinweis: Lege IGC-Dateien in '{default_igc_dir}' ab oder nutze --dir.",)
        sys.exit(0)

    print(f"Gefundene Lufträume: {len(airspaces)} | IGC-Dateien: {len(igc_files)} | Filter: {args.filter} | Plots: {args.plots}\n")

    plots_dir = Path(args.plots_dir).resolve() if args.plots else None
    if plots_dir:
        plots_dir.mkdir(parents=True, exist_ok=True)

    # Output dir for CSVs
    output_dir = Path("output").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    for igc in igc_files:
        print(f"=== {igc.name} ===")
        track, meta = read_igc_track_with_meta(igc)
        if args.geoid_correct:
            track = _apply_fixed_offset(track, 50.0)
        if len(track) < 2:
            print("  (Track zu kurz oder keine B-Records)"); continue

        results = analyze_track_with_index(track, airspaces, plots=args.plots, plots_dir=plots_dir, filter_mode=args.filter,
                                           basemap=args.basemap, basemap_provider=args.basemap_provider)

        if results:
            out_csv = output_dir / f"{igc.stem}.airspace_events.csv"
            with out_csv.open("w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow(["airspace_name","airspace_type","class","icao_class_raw","entry_utc","exit_utc","duration_s",
                                 "entry_lon","entry_lat","exit_lon","exit_lat","lower_m","upper_m","max_alt_m"])
                for r in results:
                    writer.writerow([
                        r["airspace_name"], r["airspace_type"], r["class"], r["icao_class_raw"],
                        r["entry"].strftime("%Y-%m-%d %H:%M:%S UTC"),
                        r["exit"].strftime("%Y-%m-%d %H:%M:%S UTC"),
                        int(round(r["duration_s"])),
                        f"{r['entry_lon']:.6f}", f"{r['entry_lat']:.6f}", f"{r['exit_lon']:.6f}", f"{r['exit_lat']:.6f}",
                        f"{r['lower_m'] if math.isfinite(r['lower_m']) else ''}",
                        f"{r['upper_m'] if math.isfinite(r['upper_m']) else ''}",
                        f"{r.get('max_alt_m','')}"
                    ])
            print(f"  → {len(results)} Ereignis(se). CSV: {out_csv}")
        else:
            print("  → Keine Einflüge gefunden. Keine CSV erzeugt.")
        if plots_dir:
            print(f"  → PNGs (falls vorhanden): {plots_dir}")

if __name__ == "__main__":
    main()
