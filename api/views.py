from dataclasses import field
from logging import Filter
import pytz
import shutil
import requests
import urllib.request
from pathlib import Path
from datetime import date, datetime, timedelta, timezone

from rest_framework import generics
from rest_framework.filters import SearchFilter, OrderingFilter
from django_filters import rest_framework as filters
from django.conf import settings
from django.http import HttpResponse
from django.db.models import Prefetch

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
    queryset = Event.objects.all()
    serializer_class = EventSerializer
    filter_backends = [filters.DjangoFilterBackend, SearchFilter, OrderingFilter]
    filterset_fields = ['type', 'start_time', 'end_time']
    search_fields = ['spec_id']
    ordering_fields = ['start_time', 'end_time']
    ordering = ['start_time']
    

#2 ������� ���������� ����������/ ��� �������� ���������
def load_HEK_CH():
    g_CH_load_status = False

    def ast_utc(obj):
        return obj.datetime.astimezone(timezone.utc)
    
    load_start_time = Time('2023-01-01T00:00:00', scale='utc', format='isot')
    load_end_time = Time(datetime.now())
    hek_client = hek.HEKClient()
    responses = hek_client.search(
        attrs.Time(load_start_time, load_end_time), 
        attrs.hek.CH, 
        attrs.hek.FRM.Name == 'SPoCA',
    )
    responses = sorted(responses, key=lambda x: x['event_starttime'])
    for index in range(len(responses)):
        ch = responses[index] # JSON dict
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
            Point.objects.create(
                phi=float(point[0]),
                theta=float(point[1]), 
                polyline=polyline,
            )
    g_CH_load_status = True
    return HttpResponse(g_CH_load_status)


def load_STOP_PFSS_lines():
    g_PML_load_status = False

    def fits_url(cr):
        url = 'http://158.250.29.123:8000/web/Stop/Synoptic%20maps/'
        return f'{url}stop_{cr}.fits'
    start_cr = int(carrington_rotation_number(date(2023,1,1)))
    end_cr = int(carrington_rotation_number(date.today()))
    end_cr_exists = requests.head(fits_url(end_cr)).status_code == 200
    end_cr += 1 if end_cr_exists else 0
    nrho, rss, r, divisor = 35, 2.5, 2.5 * const.R_sun, 16
    for cr_ind in range(start_cr, end_cr):
        path = settings.BASE_DIR / f'Media/stop_{cr_ind}.fits'
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
        car_data = carrington_rotation_time(cr_ind).to_datetime(timezone=pytz.timezone(settings.TIME_ZONE))
        next_car_data = carrington_rotation_time(cr_ind + 1).to_datetime(timezone=pytz.timezone(settings.TIME_ZONE))
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
        def polarity_convector(string):
            return None if string == '0' else int(string) > 0
            
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
                        phi=lat,
                        theta=lon,
                        r=item.radius / const.R_sun,
                        polyline=polyline,
                    )
    g_PML_load_status = True
    return HttpResponse(g_PML_load_status)
    
#2.1 �������� � ��������� �������
