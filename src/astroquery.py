from astroquery.jplhorizons import Horizons
from astropy.time import Time
from typing import Tuple, Dict, List
import numpy as np
from src.utils import load_parameters

PARAMS = load_parameters()


def fetch_planet(body_id: int, time: str) -> Tuple[np.ndarray, np.ndarray]:
    obj = Horizons(
        id=body_id, location="@sun", epochs=Time(time).jd, id_type="id"
    ).vectors()
    p = np.array([float(obj[i]) for i in ["x", "y", "z"]])
    v = np.array([float(obj[i]) for i in ["vx", "vy", "vz"]])
    return p, v


def fetch_vectors(names: List[str], time: str) -> Tuple[Dict[str, np.ndarray]]:
    positions = dict()
    velocities = dict()
    for name in names:
        body_id = PARAMS["id"][name]
        p, v = fetch_planet(body_id=body_id, time=time)
        positions[name] = p
        velocities[name] = v

    return positions, velocities
