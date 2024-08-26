import random
import shutil
import requests
from datetime import date, datetime, time, timedelta, timezone

import numpy as np
from scipy.spatial.distance import pdist, squareform

import matplotlib as mpl
from matplotlib.figure import Figure
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
from sunpy.coordinates.frames import HeliographicCarrington
from sunpy.map import Map

from pfsspy import tracing, utils
from pfsspy.utils import is_full_sun_synoptic_map

from django.http import HttpResponse
from django.conf import settings

from .common import *
from .models import (CoronalHole, CoronalHolePoint,
    CoronalHoleContour, CoronalHoleContourPoint,
    MagneticLineSet, MagneticLine, MagneticLinePoint)


test_date = date(2024, 8, 1)
zero_time = time(0, 0, 0, tzinfo=timezone.utc)

def load_HEK_CH():

    def ast_utc(obj):
        return obj.datetime.astimezone(timezone.utc)

    try:
        last_hole = CoronalHole.objects.filter(
            s_type='HCH').latest('start_time')
        load_start_time = last_hole.start_time + timedelta(seconds=1)
    except CoronalHole.DoesNotExist:
        load_start_time = Time(datetime.combine(test_date, zero_time))
    load_end_time = Time(datetime.now())
    hek_client = hek.HEKClient()
    responses = hek_client.search(
        attrs.Time(load_start_time, load_end_time),
        attrs.hek.CH,
        attrs.hek.FRM.Name == 'SPoCA',
    )
    responses.sort(['event_starttime'])
    for event in responses:
        lon, lat = map(float, event['hgc_coord'][6:-1].split())
        ch = CoronalHole.objects.create(
            sol=event['solar_object_locator'],
            start_time=ast_utc(event['event_starttime']),
            end_time=ast_utc(event['event_endtime']),
            lon=lon, lat=lat, s_type='HCH')
        polygon_string = event['hgc_boundcc'][9:-2]
        polygon_list = [item.split(' ') for item in polygon_string.split(',')]
        contour = CoronalHoleContour.objects.create(ch=ch)
        chc_pts_batch = []
        for point in polygon_list:
            lon, lat = float(point[0]), float(point[1])
            chc_point = CoronalHoleContourPoint(
                lon=lon, lat=lat, contour=contour)
            chc_pts_batch.append(chc_point)
        CoronalHoleContourPoint.objects.bulk_create(chc_pts_batch)
    return HttpResponse(True)


def load_STOP():
    session = requests.Session()

    def fits_url(cr):
        url = 'http://158.250.29.123:8000/web/Stop/Synoptic%20maps/'
        return f'{url}stop_{cr}.fits'
    
    def fits_exists(cr):
        return session.head(fits_url(cr)).status_code == 200

    try:
        last_pml = MagneticLineSet.objects.filter(
            s_type='SML').latest('start_time')
        start_cr = int(carrington_rotation_number(last_pml.start_time)) + 1
    except MagneticLineSet.DoesNotExist:
        start_cr = int(carrington_rotation_number(test_date)) - 1
    end_cr = int(carrington_rotation_number(date.today()))
    ph_path = settings.BASE_DIR / f'maps/synoptic/photospheric/stop'
    ph_path.mkdir(parents=True, exist_ok=True)
    while not fits_exists(end_cr): end_cr -= 1
    for cr_ind in range(start_cr, end_cr + 1):
        fits_fname = ph_path / f'{cr_ind}.fits'
        if not fits_fname.exists():
            if not fits_exists(cr_ind): continue
            response = session.get(fits_url(cr_ind), stream=True)
            with open(fits_fname, 'wb') as out_file:
                response.raw.decode_content = True
                shutil.copyfileobj(response.raw, out_file)
        full_plot(fits_fname, 'stop', carrot=cr_ind)
    return HttpResponse(True)

def load_STOP_daily():
    session = requests.Session()

    def date_range(start_date: date, end_date: date):
        days = int((end_date - start_date).days)
        for n in range(days + 1):
            yield start_date + timedelta(n)

    def fits_url(fits_date):
        url = 'http://158.250.29.123:8000/web/Stop/Synoptic%20maps%20daily/'
        return f'{url}stop_blr0_{fits_date:%y%m%d}.fits'
    
    def fits_exists(fits_date):
        return session.head(fits_url(fits_date)).status_code == 200

    try:
        last_sch = CoronalHole.objects.filter(
            s_type='SCH').latest('start_time')
        start_date = last_sch.start_time.date() + timedelta(days=1)
    except CoronalHole.DoesNotExist:
        start_date = test_date
    end_date = date.today()
    ph_path = settings.BASE_DIR / 'maps/synoptic/daily/photospheric/stop'
    ph_path.mkdir(parents=True, exist_ok=True)
    for fits_date in date_range(start_date, end_date):
        fits_fname = ph_path / f'{fits_date:%y%m%d}.fits'
        if not fits_fname.exists():
            if not fits_exists(fits_date): continue
            response = session.get(fits_url(fits_date), stream=True)
            with open(fits_fname, 'wb') as out_file:
                response.raw.decode_content = True
                shutil.copyfileobj(response.raw, out_file)
        full_plot(fits_fname, 'stop', fits_date, True)
    return HttpResponse(True)

def full_plot(fits_fname, m_type, fits_date=None, plot_CH=False, carrot=None):
    uid = random.randint(100000, 999999) # Time execution info
    exec_time = datetime.now()
    print(f'plot ML {uid:10d} started at ', exec_time)
        
    def fits2map(fits_fname): # m_type, fits_date
        data, header = fits.getdata(fits_fname, header=True)
        if m_type == 'stop':
            data = np.flip(data, 0)
            header['CUNIT1'] = 'deg'
            header['CUNIT2'] = 'deg'
            header['CDELT1'] = 1 / 2
            header['CDELT2'] = 1 / 2
            header['CTYPE1'] = 'CRLN-CAR'
            header['CTYPE2'] = 'CRLT-CAR'
            header['CRVAL1'] = 180
            if fits_date:
                header['DATE'] = f'{fits_date:%Y-%m-%d}'
        return Map(data, header)
    
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
        is_full_sun_synoptic_map(m, error=True)
        is_cea_map(m, error=True)
        header_out = m.wcs.to_header()
        header_out['CTYPE1'] = header_out['CTYPE1'][:5] + 'CAR'
        header_out['CTYPE2'] = header_out['CTYPE2'][:5] + 'CAR'
        header_out['CDELT2'] = 1 / 2
        wcs_out = WCS(header_out, fix=False)
        data_out = reproject(m, wcs_out, shape_out=m.data.shape,
            return_footprint=False)
        return Map(data_out, header_out)

    def source_surface_pils(ss_map):
        contours = measure.find_contours(ss_map.data, 0)
        contours = [ss_map.wcs.pixel_to_world(c[:, 1], c[:, 0]) for c in contours]
        return contours

    def polarity_convector(value):
        return None if value == 0 else int(value) > 0
    
    def get_seeds(r_coef, divisor, frame):
        r = r_coef * const.R_sun
        lon_1d = np.linspace(0, 2 * np.pi, 2 * divisor)
        lat_1d = np.linspace(-np.pi / 2, np.pi / 2, divisor)
        lon, lat = np.meshgrid(lon_1d, lat_1d, indexing='ij')
        lon, lat = (lon * u.rad).ravel(), (lat * u.rad).ravel()
        return SkyCoord(lon, lat, r, frame=frame)
    
    def spherical_metric(point1, point2):
        lam = np.array([point1[0], point2[0]])
        phi = np.array([point1[1], point2[1]])
        lam *= np.pi / 180 # radian convertion
        phi *= np.pi / 180
        d_phi = phi[1] - phi[0]
        d_lam = lam[1] - lam[0]
        hav_theta = .5 * (1 - np.cos(d_phi) + np.cos(phi[0]) * np.cos(phi[1]) * (1 - np.cos(d_lam)))
        return 2 * np.arcsin(hav_theta**.5) * 180 / np.pi
    
    def median(cluster):
        v = pdist(cluster, metric=spherical_metric)
        m = squareform(v)
        min_dist = np.inf
        for i in range(len(cluster)):
            dist_sum = sum(m[i])
            if dist_sum < min_dist:
                min_dist, ind = dist_sum, i
        return cluster[ind]
    
    def prob_reduce(data, limit=np.inf):
        size = len(data)
        if size > limit:
            reduction_factor = size / limit
            probability = 1 / reduction_factor
            mask = np.random.rand(size) < probability
            data = data[mask] # not uniform
        return data
    
    def expand_data(data, ratio, unpack=False):
        if data.ndim > 1:
            exp_data = np.array(data)
        else:
            exp_data = np.array([data])
            unpack = True
        exp_data[:, 0] += 180 # longitude
        exp_data[:, 1] += 90 # latitude
        exp_data *= ratio
        if unpack:
            exp_data = exp_data[0]
        return exp_data.astype(int)
    
    def narrow_data(data, ratio):
        data /= ratio
        data[:, 1] -= 180 # longitude
        data[:, 0] -= 90 # latitude

    ph_map = fits2map(fits_fname)
    height, width = ph_map.data.shape
    data_ratio = height / width
    path = settings.BASE_DIR / 'media/synoptic/'

    # plot photospheric figure
    fig = Figure()
    ax = fig.add_subplot(projection=ph_map)
    norm = ImageNormalize(stretch=HistEqStretch(ph_map.data))
    ph_map.plot(axes=ax, cmap='bwr', norm=norm)
    norm = cls.SymLogNorm(vmin=ph_map.min(), vmax=ph_map.max(),
        linthresh=1)
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap='bwr'),
        ax=ax, fraction=0.047*data_ratio)
    cbar.ax.set_ylabel(r'$B_{r}$, G', rotation=-90)
    title = 'Photospheric magnetogram, '
    if fits_date:
        title += f'date {fits_date:%Y-%m-%d}'
        ph_name = f'PH_{fits_date:%y%m%d}.png'
        path = path / 'daily'
    if carrot:
        title += f'CR {carrot}'
        ph_name = f'PH_{carrot}.png'
    ax.set_title(title)
    ax.set_xlabel('Carrington longitude')
    ax.set_ylabel('Carrington latitude')
    ph_path = path / f'photospheric/{m_type}/'
    ph_path.mkdir(parents=True, exist_ok=True)
    fig.savefig(ph_path / ph_name, bbox_inches='tight')
    print(f'PH fig {uid:10d} saved in ', datetime.now() - exec_time)
    exec_time = datetime.now()

    nrho, rss = 35, 2.5
    cea_ph_map = utils.car_to_cea(ph_map)
    pfss_in = utils.pfsspy.Input(cea_ph_map, nrho, rss)
    pfss_out = utils.pfsspy.pfss(pfss_in)
    cea_ss_map = pfss_out.source_surface_br
    ss_map = cea_to_car(cea_ss_map)
    print(f'PFSS {uid:10d} calculated in ', datetime.now() - exec_time)
    exec_time = datetime.now()
    
    # plot source surface figure
    fig = Figure()
    ax = fig.add_subplot(projection=ss_map)
    for item in source_surface_pils(ss_map):
        ax.plot_coord(item, 'k')
    norm = cls.CenteredNorm()
    ss_map.plot(axes=ax, cmap='bwr', norm=norm)
    cbar = fig.colorbar(mpl.cm.ScalarMappable(cmap='bwr', norm=norm),
        ax=ax, fraction=0.047*data_ratio)
    cbar.ax.set_ylabel(r'$B_{r}$, G', rotation=-90)
    title = 'Source surface magnetogram, '
    if fits_date:
        title += f'date {fits_date:%Y-%m-%d}'
        ss_name = f'SS_{fits_date:%y%m%d}.png'
    if carrot:
        title += f'CR {carrot}'
        ss_name = f'SS_{carrot}.png'
    ax.set_title(title)
    ax.set_xlabel('Carrington longitude')
    ax.set_ylabel('Carrington latitude')
    ss_path = path / f'source_surface/{m_type}/'
    ss_path.mkdir(parents=True, exist_ok=True)
    fig.savefig(ss_path / ss_name, bbox_inches='tight')
    print(f'SS fig {uid:10d} saved in ', datetime.now() - exec_time)
    exec_time = datetime.now()

    # plot solar wind figure
    W_s = 6.25 * (ss_map.data / ph_map.data) ** 2
    sw_data = 393.2 + 192.9 * W_s + 3.94 * abs(ss_map.data) - 0.019 * abs(ph_map.data)
    sw_map = Map(sw_data, ph_map.meta)

    fig = Figure()
    ax = fig.add_subplot(projection=sw_map)
    norm = ImageNormalize(stretch=HistEqStretch(sw_map.data))
    sw_map.plot(axes=ax, cmap='RdBu_r', norm=norm)
    norm = cls.Normalize(vmin=250, vmax=750)
    cbar = fig.colorbar(mpl.cm.ScalarMappable(cmap='RdBu_r', norm=norm),
        ax=ax, fraction=0.047*data_ratio)
    cbar.ax.set_ylabel(r'$V_{r}$, Km/s', rotation=-90, labelpad=15)
    title = 'Source surface solar wind, '
    if fits_date:
        title += f'date {fits_date:%Y-%m-%d}'
        sw_name = f'SW_{fits_date:%y%m%d}.png'
    if carrot:
        title += f'CR {carrot}'
        sw_name = f'SW_{carrot}.png'
    ax.set_title(title)
    ax.set_xlabel('Carrington longitude')
    ax.set_ylabel('Carrington latitude')
    sw_path = path / f'solar_wind/{m_type}/'
    sw_path.mkdir(parents=True, exist_ok=True)
    fig.savefig(sw_path / sw_name, bbox_inches='tight')
    print(f'SW fig {uid:10d} saved in ', datetime.now() - exec_time)
    exec_time = datetime.now()

    tracer = tracing.FortranTracer(max_steps=3000)
    ss_seeds = get_seeds(2.5, 16, pfss_out.coordinate_frame)
    ss_field_lines = tracer.trace(ss_seeds, pfss_out)
    print(f'SS seed ML {uid:10d} traced in ', datetime.now() - exec_time)
    exec_time = datetime.now()

    if fits_date:
        start_time = datetime.combine(fits_date, zero_time) # datetime.fromisoformat(ph_map.meta['DATE'])
        end_time = start_time + timedelta(days=1)
        obstime = start_time + timedelta(hours=8) # no obstime!
    if carrot:
        start_time = carrington_rotation_time(
            carrot).to_datetime(timezone=timezone.utc)
        end_time = carrington_rotation_time(
            carrot + 1).to_datetime(timezone=timezone.utc)
        obstime = end_time
    if m_type == 'stop':
        ML_type = 'dSML' if fits_date else 'SML'
    elif m_type == 'gong':
        ML_type = 'dGML' if fits_date else 'GML'
    lineset = MagneticLineSet.objects.create(
        s_type=ML_type, start_time=start_time)
    
    lines_batch, points_batch = [], []
    for field_line in ss_field_lines.open_field_lines:
        coords = field_line.coords
        if len(coords) > 1:
            car_coords = coords.transform_to(HeliographicCarrington(obstime=obstime))
            line = MagneticLine(lineset=lineset,
                polarity=polarity_convector(field_line.polarity))
            lines_batch.append(line)
            for item in car_coords:
                lon = item.lon.wrap_at('180d').degree
                lat = item.lat.degree
                point = MagneticLinePoint(
                    lon=lon,
                    lat=lat,
                    r=item.radius / const.R_sun,
                    line=line,
                )
                points_batch.append(point)
        
    MagneticLine.objects.bulk_create(lines_batch)
    MagneticLinePoint.objects.bulk_create(points_batch)
    lines_batch, points_batch = [], []
    print(f'SS seed ML {uid:10d} saved in ', datetime.now() - exec_time)
    exec_time = datetime.now()
    
    if plot_CH:
        ph_seeds = get_seeds(1, height, pfss_out.coordinate_frame)
        ph_field_lines = tracer.trace(ph_seeds, pfss_out)
        print(f'PH seed ML {uid:10d} traced in ', datetime.now() - exec_time)
        exec_time = datetime.now()

        open_field_coords = []
        ratio = height / 180
        for field_line in ph_field_lines.open_field_lines:
            coords = field_line.solar_footpoint
            car_coords = coords.transform_to(HeliographicCarrington(obstime=obstime))
            lon = car_coords.lon.wrap_at('180d').degree
            lat = car_coords.lat.degree
            open_field_coords.append((lon, lat))
        open_field_coords = np.array(open_field_coords)
        print(f'polarity matix {uid:10d} calculated in ', datetime.now() - exec_time)
        exec_time = datetime.now()
        
        db = DBSCAN(eps=5, min_samples=10,
            metric=spherical_metric).fit(open_field_coords)
        labels = db.labels_
        print(f'DBSCAN {uid:10d} calculated in ', datetime.now() - exec_time)
        exec_time = datetime.now()

        unique_labels = set(labels)
        core_samples_mask = np.zeros_like(labels, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True

        fig = Figure(figsize=(24, 12), dpi=600)
        ax = fig.add_subplot()
        ax.set_facecolor('#808080')

        one_px_area = 1.4743437*10**14 # m^2 sun_surface/full_solid_angle
        for k in unique_labels:
            class_member_mask = labels == k
            cluster = open_field_coords[class_member_mask & core_samples_mask]
            size = len(cluster)
            if k >= 0 and size:
                ch = CoronalHole.objects.create(
                    start_time=start_time, end_time=end_time, s_type='SCH')
                cluster_image = np.zeros((height, width), dtype=int)
                expanded_cluster = expand_data(cluster, ratio)
                cluster_image[expanded_cluster[:,1], expanded_cluster[:,0]] = 1
                Br = ph_map.data[expanded_cluster[:,1], expanded_cluster[:,0]]
                mag_sum = sum(Br)
                first = expanded_cluster[0]
                col = 'r' if ss_map.data[first[1], first[0]] > 0 else 'b'
                ax.scatter(cluster[:, 0], cluster[:, 1], s=1, color=col)
                ch_pts_batch = []
                batch_limit = 20000
                for ind in range(size):
                    lon, lat = cluster[ind]
                    point = CoronalHolePoint(
                        lon=lon, lat=lat, Br=Br[ind], ch=ch)
                    ch_pts_batch.append(point)
                    if len(ch_pts_batch) >= batch_limit:
                        CoronalHolePoint.objects.bulk_create(ch_pts_batch)
                        ch_pts_batch = []
                if ch_pts_batch:
                    CoronalHolePoint.objects.bulk_create(ch_pts_batch)
                    ch_pts_batch = []
                print(f'CH point {uid:10d} for {k}-label saved in ', datetime.now() - exec_time)
                exec_time = datetime.now()
                
                reduced_cluster = prob_reduce(cluster, limit=1024)
                center = median(reduced_cluster)
                ch.lon, ch.lat = center
                sol_lon, sol_lat = expand_data(center, ratio)
                ch.sol = f'SOL{obstime:%Y-%m-%dT%H:%M}L{sol_lon:.2f}C{sol_lat:.2f}'
                ch.area = size * one_px_area * 10**-12 # Mm^2
                ch.mag_flux = mag_sum * one_px_area * 10**4 # Mx
                ch.avg_flux = mag_sum / size # G
                ch.max_flux = max(Br, key=abs) # G
                ch.save()

                cluster_image = binary_closing(cluster_image, disk(3))
                contours = measure.find_contours(cluster_image, 0.5)
                for contour in contours:
                    chc = CoronalHoleContour.objects.create(ch=ch)
                    contour = prob_reduce(contour, limit=1024)
                    narrow_data(contour, ratio)
                    chc_pts_batch = []
                    for lat, lon in contour:
                        point = CoronalHoleContourPoint(
                            lon=lon, lat=lat, contour=chc)
                        chc_pts_batch.append(point)
                    CoronalHoleContourPoint.objects.bulk_create(chc_pts_batch)
                    ax.plot(contour[:, 1], contour[:, 0], linewidth=1, color='k')
                chc_pts_batch = []
                print(f'CH contours {uid:10d} for {k}-label saved in ', datetime.now() - exec_time)
                exec_time = datetime.now()
        
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
        ax.set_xticks(np.arange(-180, 181, 10))
        ax.set_yticks(np.arange(-90, 91, 10))
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.grid(True, color='k', linewidth=0.5)
        ch_path = path / f'coronal_hole/{m_type}/'
        ch_path.mkdir(parents=True, exist_ok=True)
        ch_name = f'CH_{fits_date:%y%m%d}.png'
        fig.savefig(ch_path / ch_name, bbox_inches='tight', pad_inches=0)

    print(f'CR{carrot} date:{fits_date} {uid:10d} finished at ', datetime.now())
