"""Site list shared by the TROPOMI and ERA5 download scripts.

Coordinates are approximate city/industrial-area centres (lat, lon).
Add the background control site here AFTER choosing it from a land-use map
(before looking at any NO2 result), for example:
    "control": (lat, lon),
"""
import math

SITES = {
    "algiers": (36.75, 3.06),
    "oran": (35.70, -0.63),
    "arzew": (35.85, -0.31),
    "constantine": (36.37, 6.61),
    "setif": (36.19, 5.41),
    "hassi_messaoud": (31.68, 6.07),
    # Background control sites (chosen from land cover, population and night lights before any NO2 data
    # were opened; amended selection, outputs\\control_sites_selected_amended.csv). They are NOT in
    # DEFAULT_SITES: select them with --sites control_tell,control_highland.
    "control_tell": (35.0, 5.25),
    "control_highland": (34.5, -1.0),
}

# Sites that passed (or are kept as exploratory) in the screening step.
DEFAULT_SITES = ["algiers", "oran", "arzew", "constantine", "setif", "hassi_messaoud"]

# Half-size of the box cut around each site, expressed as a radius in degrees of
# latitude (1.5 deg is about 167 km). Longitude is widened by 1/cos(lat) so the
# box is roughly square on the ground.
HALF_BOX_DEG = 1.5


def site_box(name, half=HALF_BOX_DEG):
    """Return (lon_min, lat_min, lon_max, lat_max) for a site."""
    lat, lon = SITES[name]
    half_lon = half / math.cos(math.radians(lat))
    return (lon - half_lon, lat - half, lon + half_lon, lat + half)


def parse_site_list(text):
    if not text or text.lower() == "default":
        return list(DEFAULT_SITES)
    names = [t.strip().lower() for t in text.split(",") if t.strip()]
    bad = [n for n in names if n not in SITES]
    if bad:
        raise SystemExit("Unknown site(s): %s. Known: %s" % (", ".join(bad), ", ".join(SITES)))
    return names
