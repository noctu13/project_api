import requests
import shutil
from datetime import date, datetime, timedelta, timezone

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr

import astropy.units as u
import astropy.constants as const
from astropy.time import Time
from astropy.io import fits
from astropy.visualization import HistEqStretch, ImageNormalize
from astropy.coordinates import SkyCoord

from sunpy.net import hek, attrs
from sunpy.coordinates.sun import (
    carrington_rotation_number, carrington_rotation_time)
from sunpy.map import Map
from sunpy.coordinates.frames import HeliographicStonyhurst as HGS

from pfsspy import tracing, utils

from django.http import HttpResponse
from django.conf import settings

from .common import *
from .models import Event, Polyline, Point


def load_SPOCA_CH_from_HEK():
    g_CH_load_status = False

    def ast_utc(obj):
        return obj.datetime.astimezone(timezone.utc)

    try:
        last_hole = Event.objects.filter(
            type=short_type_dict['CH']).latest('start_time')
        load_start_time = last_hole.start_time.replace(
            hour=0, minute=0, second=0)
    except Event.DoesNotExist:
        load_start_time = Time('2023-01-01T00:00:00', scale='utc', format='isot')
    load_end_time = Time(datetime.now())
    hek_client = hek.HEKClient()
    responses = hek_client.search(
        attrs.Time(load_start_time, load_end_time),
        attrs.hek.CH,
        attrs.hek.FRM.Name == 'SPoCA',
    )
    responses.sort(['event_starttime'])
    for ch in responses:
        event_start_time = ast_utc(ch['event_starttime'])
        event_end_time = ast_utc(ch['event_endtime'])
        event, created = Event.objects.get_or_create(
            type=short_type_dict['CH'],
            spec_id=ch['frm_specificid'],
        )
        if not event.start_time:
            event.start_time = event_start_time
        event.end_time = event_end_time
        event.save()
        polyline = Polyline.objects.create(
            event=event,
            start_time=event_start_time,
            end_time=event_end_time,
        )
        polygon_string = ch['hgs_boundcc'][9:-2] # strip POLYGON(())
        polygon_list = [item.split(' ') for item in polygon_string.split(',')]
        for point in polygon_list:
            lon, lat = float(point[0]), float(point[1])
            Point.objects.create(
                phi = lon,
                theta = lat,
                polyline=polyline,
            )
    g_CH_load_status = True
    return HttpResponse(g_CH_load_status)


def load_STOP_maps():
    g_PML_load_status = False

    def fits_url(cr):
        url = 'http://158.250.29.123:8000/web/Stop/Synoptic%20maps/'
        return f'{url}stop_{cr}.fits'
    
    def fits_exists(cr):
        return requests.head(fits_url(cr)).status_code == 200

    try:
        last_pml = Event.objects.filter(
            type=short_type_dict['PML']).latest('start_time')
        start_cr = int(carrington_rotation_number(last_pml.start_time)) + 1
    except Event.DoesNotExist:
        start_cr = int(carrington_rotation_number(date(2023,1,1)))
    end_cr = int(carrington_rotation_number(date.today()))
    while not fits_exists(end_cr): end_cr -= 1
    for cr_ind in range(start_cr, end_cr + 1):
        path = settings.BASE_DIR / f'maps/synoptic/photospheric/stop/{cr_ind}.fits' # add carrington folder?
        if not path.exists() and fits_exists(cr_ind):
            response = requests.get(fits_url(cr_ind), stream=True)
            with open(path, 'wb') as out_file:
                response.raw.decode_content = True
                shutil.copyfileobj(response.raw, out_file)
            plot_STOP_PML(path)
    g_PML_load_status = True
    return HttpResponse(g_PML_load_status)

def load_daily_STOP_maps():

    def date_range(start_date: date, end_date: date):
        days = int((end_date - start_date).days)
        for n in range(days):
            yield start_date + timedelta(n)

    def fits_url(fits_date):
        url = 'http://158.250.29.123:8000/web/Stop/Synoptic%20maps%20daily/'
        return f'{url}stop_blr0_{fits_date:%y%m%d}.fits'
    
    def fits_exists(fits_date):
        return requests.head(fits_url(fits_date)).status_code == 200

    start_date = date(2024, 1, 1)
    end_date = date.today()
    for fits_date in date_range(start_date, end_date):
        path = settings.BASE_DIR / f'maps/synoptic/daily/photospheric/stop/{fits_date:%y%m%d}.fits'
        if not path.exists() and fits_exists(fits_date):
            response = requests.get(fits_url(fits_date), stream=True)
            with open(path, 'wb') as out_file:
                response.raw.decode_content = True
                shutil.copyfileobj(response.raw, out_file)
            plot_STOP_PML(path, fits_date=fits_date)

def plot_STOP_PML(path, fits_date=None):

    def polarity_convector(string):
        return None if string == '0' else int(string) > 0
    
    with fits.open(path) as hdul:
        data = hdul[0].data
        header = hdul[0].header
    data = np.flip(data, 0)
    data_ratio = data.shape[0]/data.shape[1]
    header['CUNIT1'] = 'deg'
    header['CUNIT2'] = 'deg'
    header['CDELT1'] = 1 / 2
    header['CDELT2'] = 1 / 2
    header['CTYPE1'] = 'CRLN-CAR'
    header['CTYPE2'] = 'CRLT-CAR'
    header['CRVAL1'] = 180
    if fits_date: # daily event end time (next event start) problem!
        start_time = datetime(fits_date, tzinfo=timezone.utc)
        end_time = start_time + timedelta(days=1)
    else:
        cr_ind = header['CARROT']
        start_time = carrington_rotation_time(cr_ind).to_datetime(timezone=timezone.utc)
        end_time = carrington_rotation_time(cr_ind + 1).to_datetime(timezone=timezone.utc)
    stop_map = Map(data, header)

    norm = ImageNormalize(stretch=HistEqStretch(data))
    fig = plt.figure()
    ax = fig.add_subplot(projection=stop_map)
    stop_map.plot(cmap='bwr', norm=norm, axes=ax)
    norm = clr.SymLogNorm(linthresh=1, vmin=stop_map.min(), vmax=stop_map.max())
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='bwr'),
        ax=ax, fraction=0.047*data_ratio)
    title = f'Photospheric magnetogram, '
    path = settings.BASE_DIR / 'media/synoptic/'
    if fits_date:
        title += f'date {fits_date:%Y-%m-%d}' # read from FITS obs_time?
        path = path / 'daily'
        ph_name = f'PH_{fits_date:%y%m%d}.png'
    else:
        title += f'CR {cr_ind}'
        ph_name = f'PH_{cr_ind}.png'
    ph_path = path / 'photospheric/stop/'
    plt.title(title)
    plt.savefig(ph_path / ph_name, bbox_inches='tight')

    nrho, rss, r, divisor = 35, 2.5, 2.5 * const.R_sun, 16
    stop_map = utils.car_to_cea(stop_map)
    pfss_in = utils.pfsspy.Input(stop_map, nrho, rss)
    pfss_out = utils.pfsspy.pfss(pfss_in)

    ss_br = pfss_out.source_surface_br
    fig = plt.figure()
    ax = plt.subplot(projection=ss_br)
    ss_br.plot(cmap='bwr')
    for item in pfss_out.source_surface_pils:
        ax.plot_coord(item, 'k')
    plt.colorbar(fraction=0.047*data_ratio)
    title = f'Source surface magnetogram, '
    if fits_date:
        ss_name = f'SS_{fits_date:%y%m%d}.png'
    else:
        title += f'CR {cr_ind}'
        ss_name = f'SS_{cr_ind}.png'
    ax.set_title(title)
    ss_path = path / 'source_surface/stop/'
    plt.savefig(ss_path / ss_name, bbox_inches='tight')

    tracer = tracing.FortranTracer()
    lat = np.linspace(-np.pi / 2, np.pi / 2, divisor, endpoint=False)
    lon = np.linspace(0, 2 * np.pi, 2 * divisor, endpoint=False)
    lat, lon = np.meshgrid(lat, lon, indexing='ij')
    lat, lon = lat.ravel() * u.rad, lon.ravel() * u.rad
    seeds = SkyCoord(lon, lat, r, frame=pfss_out.coordinate_frame)
    field_lines = tracer.trace(seeds, pfss_out)
    event = Event.objects.create(
        type=short_type_dict['dPML' if fits_date else 'PML'],
        start_time = start_time,
        end_time = end_time,
    )
    for field_line in field_lines.open_field_lines:
        coords = field_line.coords
        coords.representation_type = 'spherical'
        coords.transform_to(HGS)
        if len(coords) > 1:
            polyline = Polyline.objects.create(
                event=event,
                start_time=start_time,
                end_time=end_time,
                polarity=polarity_convector(field_line.polarity),
            )
            for item in coords:
                lat = float(item.lat.degree / u.deg)
                lon = float(item.lon.degree / u.deg)
                Point.objects.create( # PML order by id
                    phi=lon,
                    theta=lat,
                    r=item.radius / const.R_sun,
                    polyline=polyline,
                )