import numpy as np
from astropy.time import Time
from sunpy.net import hek, attrs
from sunpy.time import parse_time
from sunpy.coordinates.sun import carrington_rotation_number
from datetime import datetime, date
import requests
import urllib.request
import shutil
from pathlib import Path


def load():
    BASE_DIR = Path(__file__).resolve().parent.parent
    print(BASE_DIR / 'some_string' + 'string')
load()