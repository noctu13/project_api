import numpy as np
from astropy.time import Time
from sunpy.net import hek, attrs
from sunpy.time import parse_time
from datetime import datetime

from .models import Event, Polyline, Point
from .common import short_type_dict


def load():
    load_start_time = Time('2023-01-01T00:00:00', scale='utc', format='isot')
    load_end_time = Time('2023-01-01T00:00:00', scale='utc', format='isot')
    #Time(datetime.now()) # timedelta UTC+3
    hek_client = hek.HEKClient()
    responses = hek_client.search(
        attrs.Time(load_start_time, load_end_time), 
        attrs.hek.CH, 
        attrs.hek.FRM.Name == 'SPoCA',
    )
    responses = sorted(responses, key=lambda x: x['event_starttime'])
    for index in range(len(responses)):
        ch = responses[index] # JSON dict
        event_start_time = ch['event_starttime']
        event_end_time = ch['event_endtime']
        event, created = Event.objects.get_or_create(
            type=short_type_dict['CH'],
            spec_id=ch['frm_specificid'],
        )
        if not event.start_time:
            event.start_time = event_start_time
        event.end_time = event_end_time
        event.save()
        polyline = Polyline.objects.create(
            event=event.id,
            start_time=event_start_time,
            end_time=event_end_time,
        )
        polygon_string = ch["hgs_boundcc"][9:-2] # strip POLYGON(())
        polygon_list = polygon_string.split(',') # point strings
        polygon_list = [item.split(" ") for item in polygon_list] # split items to pairs
        point_list = []
        for point in polygon_list:
            point_list.append(
                Point.objects.create( # PML order by id
                    theta=float(point[0]), 
                    phi=float(point[1]),
                    polyline=polyline,
                )
            )
load()