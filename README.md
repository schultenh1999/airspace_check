# Features

- Datenquellen (Standardpfade)
    - GeoJSON: alle *.geojson in ./ressources/
    - IGC: alle *.IGC in ./input/
    - CSV-Ausgabe: ./output/
    - Ereignis-Plots (PNG): ./event_plots/ (wenn --plots aktiv)
- Filterlogik
    - Controlled: Klassen A–D (E wird unterdrückt)
    - Restricted: Namen wie ED-R…, RESTRICTED, R-…
    - TMZ und RMZ: Namens-/Typ-Erkennung
    - Kombinierte Modi (siehe --filter)
- Höhenlogik
    - Luftraum-Höhen (Lower/Upper) aus OpenAIP (GND/SFC, FL, ft, m)
    - Einfache Höhenkorrektur: fester Offset 50 m (GPS h → MSL H = h − 50 m); deaktivierbar
- Erkennung
    - Schnelle Kandidatensuche via Shapely STRtree
    - Eintritts-/Austrittszeitpunkt linear entlang Segment interpoliert
    - Pro Verletzung wird maximal erreichte Höhe im Event erfasst
- Ausgaben
    - CSV pro IGC (nur wenn Events gefunden):
      - Felder: airspace_name, airspace_type, class, icao_class_raw, entry_utc, exit_utc, duration_s, entry_lon, entry_lat, exit_lon, exit_lat, lower_m, upper_m, max_alt_m
    - PNG-Plots je Event:
      - Hauptplot: Airspace-Grenze + Track-Ausschnitt + Entry/Exit + Basemap
      - Detailplot: Flugweg vor/während/nach der Verletzung (Fokus auf Kernsegment), + Airspace-Grenze, Basemap
        → Dateinamen z. B. 20250430T101530_EDDL_CTA.png und 20250430T101530_detail.png
- Datumserkennung
    - Robuste Erkennung aus HFDTE…-Headern (Formate wie HFDTE300423, HFDTEDATE:30/04/23, HFDTE 30-04-2023, …)

# Installation

```
pip install shapely matplotlib contextily pyproj geojson fiona rasterio
```

# Nutzung

```
python3 airspace_check.py
```
- Lädt alle GeoJSONs aus ./ressources/
- Nimmt IGCs aus ./input/
- Schreibt CSVs nach ./output/
- Plots sind aus (aktivierbar mit --plots)

# Mögliche Argumente

```
python airspace_check.py [geojson | geojson_dir] [--dir IGC_DIR]
                         [--no-alt]
                         [--filter {both,controlled,restricted,tmzrmz,none}]
                         [--plots] [--plots-dir DIR]
                         [--no-basemap] [--basemap-provider NAME]
                         [--geoid-correct | --no-geoid-correct]
```
- geojson (Positionsargument, optional)
    - Datei ODER Ordner.
    - Wird nichts angegeben, nutzt das Skript ./ressources/de_asp.geojson als Marker und lädt in dem Fall alle *.geojson in ./ressources/.
    - Gibst du einen Ordner an, werden alle *.geojson darin geladen.
    - Gibst du eine einzelne Datei an, wird nur diese genutzt.
- --dir IGC_DIR (Default: ./input)
    - Ordner mit den IGC-Dateien. IGCs werden case-insensitive gefunden und dedupliziert.
- --no-alt
    - Höhenlimits ignorieren (nur 2D-Prüfung). Setzt lower_m = -∞, upper_m = +∞.
- --filter both|controlled|restricted|tmzrmz|none (Default: both)
    - controlled: nur A–D (E ist bewusst unterdrückt)
    - restricted: nur ED-R/RESTRICTED/R-…
    - tmzrmz: nur TMZ/RMZ
    - both: controlled ∪ restricted ∪ TMZ ∪ RMZ
    - none: keine Filterung (alle Airspaces)
- --plots
    - Erstellt PNG-Plots pro Event (Hauptplot und Detailplot).
    - Zielordner via --plots-dir.
- --plots-dir DIR (Default: event_plots)
    - Ablageort für PNGs; wird bei Bedarf angelegt.
- --no-basemap
    - OSM-Basemap ausschalten (falls contextily/pyproj fehlen, wird automatisch ohne Basemap geplottet).
- --basemap-provider NAME (Default: CartoDB.Positron)
    - Provider-Bezeichner für contextily, z. B.:
    - OSM, CartoDB.Positron, Stamen.Terrain, Stamen.TonerLite …
    - (Dot-Pfad ist möglich, analog ctx.providers.X.Y.)
- --geoid-correct / --no-geoid-correct (Default: --geoid-correct)
    - Schaltet die feste 50-m-Höhenkorrektur ein/aus (H = h − 50 m).
    - Das ist eine vereinfachte Annäherung für GPS-→MSL.


