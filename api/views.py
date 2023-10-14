import pytz
import shutil
import requests
import urllib.request
from pathlib import Path
from datetime import date, datetime, timedelta, timezone

from rest_framework import generics
from django.conf import settings
from django.http import HttpResponse
from django.db.models import Prefetch, Q

import numpy as np
import astropy.units as u
import astropy.constants as const
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.time import Time
from sunpy.map import Map
from sunpy.net import hek, attrs
from sunpy.time import parse_time
from sunpy.coordinates.sun import (
    carrington_rotation_number, carrington_rotation_time)
from sunpy.coordinates.frames import HeliographicStonyhurst as HGS
from pfsspy import coords, tracing, utils

from .common import *
from .models import Event, Polyline, Point
from .serializers import EventSerializer, PolylineSerializer, PointSerializer
from .permissions import IsAdminOrReadOnly


class EventsView(generics.ListAPIView):
    serializer_class = EventSerializer
    
    def get_queryset(self):
        short_type = short_type_dict[self.kwargs['short_type']]
        queryset = Event.objects.filter(type=short_type)
        year = self.kwargs.get('year')
        month = self.kwargs.get('month')
        day = self.kwargs.get('day')
        event_dict = {}
        for key in short_type_dict.keys():
             event_dict[key] = short_type == short_type_dict[key]
        if year and month and day:
            query_day = datetime(year, month, day, 0, 0, 0, 
                tzinfo=timezone.utc)
            query1 = Q(start_time__lt=query_day + timedelta(days=1))
            query2 = Q(end_time__gt=query_day)
            into_query_day = Polyline.objects.filter(
                query1 & query2).order_by('start_time')
            if event_dict['PML']:
                noon = query_day + timedelta(hours=12)
                query1 = Q(start_time__lte=noon)
                query2 = Q(end_time__gt=noon)
            queryset = queryset.filter(query1 & query2).order_by(
                'start_time').prefetch_related(
                Prefetch(
                    'polyline',
                    queryset=into_query_day,
                    to_attr='into_query_day'
                )
            )
        return queryset # много элементов > 150k - ошибка

def ast_utc(obj):
    return obj.datetime.astimezone(timezone.utc)

#2 требует повышенные привелегии/ нет контроля состояния
def load_HEK_CH():
    g_CH_load_status = False
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


def load_STOP_PFSS_lines():
    g_PML_load_status = False
    nrho, rss, r, divisor = 35, 2.5, 2.5 * const.R_sun, 16

    def fits_url(cr):
        url = 'http://158.250.29.123:8000/web/Stop/Synoptic%20maps/'
        return f'{url}stop_{cr}.fits'

    def cr_exists(cr):
        return requests.head(fits_url(cr)).status_code == 200
    
    def polarity_convector(string):
        return None if string == '0' else int(string) > 0
    
    try:
        last_pml = Event.objects.filter(
            type=short_type_dict['PML']).latest('start_time')
        start_cr = int(carrington_rotation_number(last_pml.start_time)) + 1
    except Event.DoesNotExist:
        start_cr = int(carrington_rotation_number(date(2023,1,1)))
    end_cr = int(carrington_rotation_number(date.today()))
    while not cr_exists(end_cr): end_cr -= 1
    for cr_ind in range(start_cr, end_cr + 1):
        path = settings.BASE_DIR / f'media/stop_{cr_ind}.fits'
        if not path.exists():
            with urllib.request.urlopen(fits_url(cr_ind)) as response, open(
                path, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)
        with fits.open(path) as hdul:
            data = hdul[0].data
            data = np.flip(data, 0)
            header = hdul[0].header
        header['CUNIT1'] = 'deg'
        header['CUNIT2'] = 'deg'
        header['CDELT2'] = 1 / np.pi
        header['CTYPE1'] = 'CRLN-CEA'
        header['CTYPE2'] = 'CRLT-CEA'
        header['CRVAL1'] = 180
        car_data = carrington_rotation_time(cr_ind).to_datetime(timezone=timezone.utc)
        next_car_data = carrington_rotation_time(cr_ind + 1).to_datetime(timezone=timezone.utc)
        stop_map = Map(data, header)
        stop_map = stop_map.resample([360, 180] * u.pix)
        pfss_in = utils.pfsspy.Input(stop_map, nrho, rss)
        pfss_out = utils.pfsspy.pfss(pfss_in)
        tracer = tracing.FortranTracer()
        lat = np.linspace(-np.pi / 2, np.pi / 2, divisor, endpoint=False)
        lon = np.linspace(0, 2 * np.pi, divisor, endpoint=False)
        lat, lon = np.meshgrid(lat, lon, indexing='ij')
        lat, lon = lat.ravel() * u.rad, lon.ravel() * u.rad
        seeds = SkyCoord(lon, lat, r, frame=pfss_out.coordinate_frame)
        field_lines = tracer.trace(seeds, pfss_out)
        event = Event.objects.create(
            type=short_type_dict['PML'],
            start_time = car_data,
            end_time = next_car_data,
        )   
        for field_line in field_lines:
            coords = field_line.coords
            coords.representation_type = 'spherical'
            coords.transform_to(HGS)
            if len(coords) > 1:
                polyline = Polyline.objects.create(
                    event=event,
                    start_time=car_data,
                    end_time=next_car_data,
                    polarity=polarity_convector(field_line.polarity),
                )
                for item in coords:
                    lat = float(str(item.lat * u.deg)[:-5])
                    lon = float(str(item.lon * u.deg)[:-5])
                    Point.objects.create( # PML order by id
                        phi=lon,
                        theta=lat,
                        r=item.radius / const.R_sun,
                        polyline=polyline,
                    )
    g_PML_load_status = True
    return HttpResponse(g_PML_load_status)
    
#2.1 создание и обработка токенов
