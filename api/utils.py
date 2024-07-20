import requests
import shutil
from datetime import date, datetime, timedelta, timezone

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as cls

from skimage import measure
from skimage.morphology import binary_closing, disk
from sklearn.cluster import DBSCAN

from reproject import reproject_adaptive, reproject_exact, reproject_interp

import astropy.units as u
import astropy.constants as const
from astropy.time import Time
from astropy.io import fits
from astropy.visualization import HistEqStretch, ImageNormalize
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

from sunpy.net import hek, attrs
from sunpy.coordinates.sun import (
    carrington_rotation_number, carrington_rotation_time)
from sunpy.map import Map
from sunpy.coordinates.frames import HeliographicStonyhurst as HGS

from pfsspy import tracing, utils

from django.http import HttpResponse
from django.conf import settings

from .common import *
from .models import (HEKCoronalHole, HEKCoronalHoleContourPoint,
    CoronalHole, CoronalHoleContour, CoronalHoleContourPoint, CoronalHolePoint, 
    MagneticLineSet, MagneticLine, MagneticLinePoint)


def load_SPOCA_CH_from_HEK():
    g_CH_load_status = False

    def ast_utc(obj):
        return obj.datetime.astimezone(timezone.utc)

    try:
        last_hole = HEKCoronalHole.objects.latest('start_time')
        load_start_time = last_hole.start_time + timedelta(seconds=1)
    except HEKCoronalHole.DoesNotExist:
        load_start_time = Time('2024-01-01T00:00:00', scale='utc', format='isot')
    load_end_time = Time(datetime.now())
    hek_client = hek.HEKClient()
    responses = hek_client.search(
        attrs.Time(load_start_time, load_end_time),
        attrs.hek.CH,
        attrs.hek.FRM.Name == 'SPoCA',
    )
    responses.sort(['event_starttime'])
    for ch in responses:
        event = HEKCoronalHole.objects.create(
            spec_id=ch['frm_specificid'],
            sol=ch['solar_object_locator'],
            start_time=ast_utc(ch['event_starttime']),
            end_time=ast_utc(ch['event_endtime']),
        )
        polygon_string = ch['hgc_boundcc'][9:-2]
        polygon_list = [item.split(' ') for item in polygon_string.split(',')]
        for point in polygon_list:
            lon, lat = float(point[0]), float(point[1])
            HEKCoronalHoleContourPoint.objects.create(
                lon=lon,
                lat=lat,
                ch=event,
            )
    g_CH_load_status = True
    return HttpResponse(g_CH_load_status)


def load_STOP_maps():
    g_SML_load_status = False

    def fits_url(cr):
        url = 'http://158.250.29.123:8000/web/Stop/Synoptic%20maps/'
        return f'{url}stop_{cr}.fits'
    
    def fits_exists(cr):
        return requests.head(fits_url(cr)).status_code == 200

    try:
        last_pml = MagneticLineSet.objects.filter(
            type='SML').latest('start_time')
        start_cr = int(carrington_rotation_number(last_pml.start_time)) + 1
    except MagneticLineSet.DoesNotExist:
        start_cr = int(carrington_rotation_number(date(2024,7,1))) - 1 # test date
    end_cr = int(carrington_rotation_number(date.today()))
    while not fits_exists(end_cr): end_cr -= 1
    for cr_ind in range(start_cr, end_cr + 1):
        path = settings.BASE_DIR / f'maps/synoptic/photospheric/stop/{cr_ind}.fits' # add carrington folder?
        if not path.exists() and fits_exists(cr_ind):
            response = requests.get(fits_url(cr_ind), stream=True)
            with open(path, 'wb') as out_file:
                response.raw.decode_content = True
                shutil.copyfileobj(response.raw, out_file)
            plot_STOP_ML(path)
    g_SML_load_status = True
    return HttpResponse(g_SML_load_status)

def load_daily_STOP_maps():
    g_load_daily_STOP_status = False
    print(datetime.now())

    def date_range(start_date: date, end_date: date):
        days = int((end_date - start_date).days)
        for n in range(days):
            yield start_date + timedelta(n)

    def fits_url(fits_date):
        url = 'http://158.250.29.123:8000/web/Stop/Synoptic%20maps%20daily/'
        return f'{url}stop_blr0_{fits_date:%y%m%d}.fits'
    
    def fits_exists(fits_date):
        return requests.head(fits_url(fits_date)).status_code == 200

    try:
        last_sch = CoronalHole.objects.filter(
            type='SCH').latest('start_time')
        start_date = last_sch.start_time + timedelta(days=1)
    except CoronalHole.DoesNotExist:
        start_date = date(2024, 7, 19) # test date
    end_date = date.today()
    for fits_date in date_range(start_date, end_date):
        ph_path = f'maps/synoptic/daily/photospheric/stop/{fits_date:%y%m%d}.fits'
        fits_fname = settings.BASE_DIR / ph_path
        if not fits_fname.exists() and fits_exists(fits_date):
            response = requests.get(fits_url(fits_date), stream=True)
            with open(fits_fname, 'wb') as out_file:
                response.raw.decode_content = True
                shutil.copyfileobj(response.raw, out_file)
        fits_time = datetime(fits_date.year, fits_date.month, fits_date.day,
            tzinfo=timezone.utc)
        sch = CoronalHole.objects.filter(
            type='SCH', start_time=fits_time)
        if not sch:
            print(f'processing {fits_date:%y%m%d}.fits')
            plot_STOP_ML(fits_fname, fits_date=fits_date)
    
    g_load_daily_STOP_status = True
    print(datetime.now())
    return HttpResponse(g_load_daily_STOP_status)

def plot_STOP_ML(fits_fname, fits_date=None):

    def cea_to_car(m, method='interp'):

        def is_cea_map(m, error=False):
            return utils._check_projection(m, 'CEA', error=error)
        
        methods = {'adaptive': reproject_adaptive,
                'interp': reproject_interp,
                'exact': reproject_exact}
        if method not in methods:
            raise ValueError(f'method must be one of {methods.keys()} '
                            f'(got {method})')
        reproject = methods[method]

        from pfsspy.utils import is_full_sun_synoptic_map
        is_full_sun_synoptic_map(m, error=True)
        is_cea_map(m, error=True)

        header_out = m.wcs.to_header()
        header_out['CTYPE1'] = header_out['CTYPE1'][:5] + 'CAR'
        header_out['CTYPE2'] = header_out['CTYPE2'][:5] + 'CAR'
        header_out['CDELT2'] = 1 / 2
        #header_out['CRVAL1'] = xxx
        wcs_out = WCS(header_out, fix=False)

        data_out = reproject(m, wcs_out, shape_out=m.data.shape,
                            return_footprint=False)

        return Map(data_out, header_out)

    def source_surface_pils(ss_br):
        contours = measure.find_contours(ss_br.data, 0)
        contours = [ss_br.wcs.pixel_to_world(c[:, 1], c[:, 0]) for c in contours]
        return contours

    def polarity_convector(string):
        return None if string == '0' else int(string) > 0
    
    def get_seeds(r_coef, divisor, frame):
        r = r_coef * const.R_sun
        lon_1d = np.linspace(0, 2 * np.pi, 2 * divisor)
        lat_1d = np.linspace(-np.pi / 2, np.pi / 2, divisor)
        lon, lat = np.meshgrid(lon_1d, lat_1d, indexing='ij')
        lon, lat = (lon * u.rad).ravel(), (lat * u.rad).ravel()
        return SkyCoord(lon, lat, r, frame=frame)
    
    data, header = fits.getdata(fits_fname, header=True)
    data = np.flip(data, 0)
    data_ratio = data.shape[0]/data.shape[1]
    header['CUNIT1'] = 'deg'
    header['CUNIT2'] = 'deg'
    header['CDELT1'] = 1 / 2
    header['CDELT2'] = 1 / 2
    header['CTYPE1'] = 'CRLN-CAR'
    header['CTYPE2'] = 'CRLT-CAR'
    header['CRVAL1'] = 180
    if fits_date:
        start_time = datetime(fits_date.year, fits_date.month, fits_date.day,
            tzinfo=timezone.utc)
        end_time = start_time + timedelta(days=1)
    else:
        cr_ind = header['CARROT']
        start_time = carrington_rotation_time(
            cr_ind).to_datetime(timezone=timezone.utc)
        end_time = carrington_rotation_time(
            cr_ind + 1).to_datetime(timezone=timezone.utc)
    stop_map = Map(data, header)

    norm = ImageNormalize(stretch=HistEqStretch(stop_map.data))
    fig = plt.figure()
    ax = fig.add_subplot(projection=stop_map)
    stop_map.plot(cmap='bwr', norm=norm, axes=ax)
    norm = cls.SymLogNorm(linthresh=1, 
        vmin=stop_map.min(), vmax=stop_map.max())
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='bwr'),
        ax=ax, fraction=0.047*data_ratio)
    title_prefix = 'Photospheric magnetogram, '
    path = settings.BASE_DIR / 'media/synoptic/'
    if fits_date:
        title = f'date {fits_date:%Y-%m-%d}' # read from FITS obs_time?
        path = path / 'daily'
        ph_name = f'PH_{fits_date:%y%m%d}.png'
    else:
        title = f'CR {cr_ind}'
        ph_name = f'PH_{cr_ind}.png'
    ph_path = path / 'photospheric/stop/'
    plt.title(f'{title_prefix}{title}')
    plt.savefig(ph_path / ph_name, bbox_inches='tight')

    nrho, rss = 35, 2.5
    stop_map = utils.car_to_cea(stop_map)
    pfss_in = utils.pfsspy.Input(stop_map, nrho, rss)
    pfss_out = utils.pfsspy.pfss(pfss_in)
    ss_br = pfss_out.source_surface_br
    ss_br = cea_to_car(ss_br)

    fig = plt.figure()
    ax = plt.subplot(projection=ss_br)
    for item in source_surface_pils(ss_br):
        ax.plot_coord(item, 'k')
    norm = cls.CenteredNorm()
    ss_br.plot(cmap='bwr', norm=norm)
    plt.colorbar(fraction=0.047*data_ratio)

    title_prefix = 'Source surface magnetogram, '
    if fits_date:
        ss_name = f'SS_{fits_date:%y%m%d}.png'
    else:
        ss_name = f'SS_{cr_ind}.png'
    ss_path = path / 'source_surface/stop/'
    plt.title(f'{title_prefix}{title}')
    plt.savefig(ss_path / ss_name, bbox_inches='tight')

    tracer = tracing.FortranTracer(max_steps=3000)

    ss_seeds = get_seeds(2.5, 16, pfss_out.coordinate_frame)
    ss_field_lines = tracer.trace(ss_seeds, pfss_out)
    lineset = MagneticLineSet.objects.create(
        type='dSML' if fits_date else 'SML',
        start_time=start_time,
    )  
    for field_line in ss_field_lines.open_field_lines:
        coords = field_line.coords
        if len(coords) > 1:
            line = MagneticLine.objects.create(
                lineset=lineset,
                polarity=polarity_convector(field_line.polarity),
            )
            for item in coords:
                lat = float(item.lat.degree)
                lon = float(item.lon.degree)
                MagneticLinePoint.objects.create( # PML order by id
                    lon=lon,
                    lat=lat,
                    r=item.radius / const.R_sun,
                    line=line,
                )

    height, width = stop_map.data.shape    
    ph_seeds = get_seeds(1, height, pfss_out.coordinate_frame)
    ph_field_lines = tracer.trace(ph_seeds, pfss_out)

    pol = np.zeros((height, width)) # pol[lat][lon] = polarity
    ratio = height / 180
    for field_line in ph_field_lines.open_field_lines:
        ph_footpoint = field_line.solar_footpoint
        lon = round(ph_footpoint.lon.degree * ratio) % width
        lat = round((ph_footpoint.lat.degree + 90) * ratio) % height
        pol[lat][lon] = field_line.polarity # not unipolar for stop_map.data!
    
    Z = np.array([(lat, lon, pol[lat][lon]) for lat in range(height) for lon in range(width)])
    mask = (Z[:, 2] == 1) | (Z[:, 2] == -1)
    X = Z[mask][:, :2]
    db = DBSCAN(eps=10, min_samples=5).fit(X) #! eps=5, min_samp=10
    labels = db.labels_

    unique_labels = set(labels)
    core_samples_mask = np.zeros_like(labels, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True

    selem = disk(1)
    one_px_area = 1.4743437*10**14 # m^2 sun_surface/full_solid_angle
    for k in unique_labels:
        class_member_mask = labels == k
        xy = X[class_member_mask & core_samples_mask]
        if k != -1 and len(xy):
            ch = CoronalHole.objects.create(
                start_time=start_time,
                type='SCH',
            )

            cluster_image = np.zeros((height, width))
            mag_sum, max_flux = 0, 0
            for lat, lon in xy:
                lat, lon = int(lat), int(lon)
                Br = stop_map.data[lat][lon]
                cluster_image[lat, lon] = 1
                chp = CoronalHolePoint.objects.create(
                    lon=lon,
                    lat=lat,
                    Br=Br,
                    ch=ch,
                )
                mag_sum += Br
                max_flux = max(max_flux, abs(Br)) # no sign!
            cluster_image = binary_closing(cluster_image, selem)

            center = np.array(measure.centroid(cluster_image)) / 2
            ch.location = f'{center[0]:>6.2f} {center[1] - 90:>7.2f}'
            ch.sol = f'SOL{fits_date:%Y-%m-%dT%H:%M}L{center[0]:.2f}C{center[1]:.2f}'
            ch.area = len(xy) * one_px_area * 10**-12
            ch.mag_flux = mag_sum * one_px_area * 10**4 # cm^2 correction
            ch.avg_flux = mag_sum / len(xy) # cm^2 correction
            ch.max_flux = max_flux
            ch.save()

            contours = measure.find_contours(cluster_image, 0.5)
            chc = CoronalHoleContour.objects.create(ch=ch)
            for contour in contours:
                for point in contour:
                    lat, lon = point
                    chcp = CoronalHoleContourPoint.objects.create(
                        lon=lon,
                        lat=lat,
                        contour=chc,
                    )
            # border point processing!

# def load_GONG_maps !
# def plot_GONG_ML !