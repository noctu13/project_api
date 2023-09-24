from datetime import date, datetime
from rest_framework import generics
from django.http import HttpResponse

import numpy as np
from astropy.time import Time, TimeDelta
from sunpy.net import hek, attrs
from sunpy.time import parse_time


from .models import Event, Polyline, Point
from .serializers import EventSerializer, PolylineSerializer, PointSerializer
from .permissions import IsAdminOrReadOnly
from .common import *


class EventsView(generics.ListAPIView):
    serializer_class = EventSerializer
    
    def get_queryset(self):
        short_type = short_type_dict[self.kwargs['short_type']]
        year = self.kwargs.get('year')
        month = self.kwargs.get('month')
        day = self.kwargs.get('day')
        if year and month and day:
            query_day = date(year, month, day)
            return Event.objects.filter(
                type=short_type, 
                start_time__lte=query_day,
                end_time__gt=query_day,
            )
        return Event.objects.filter(type=short_type)

#2 get ������ �� ���������� ������ � HEK (���������� ����������)
def load_HEK_CH(request):
    result = None
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
        event_start_time = ch['event_starttime']
        event_end_time = ch['event_endtime']
        event, created = Event.objects.get_or_create(
            type=short_type_dict['CH'],
            spec_id=ch['frm_specificid'],
        )
        if not event.start_time:
            event.start_time = event_start_time.datetime
        event.end_time = event_end_time.datetime
        event.save()
        polyline = Polyline.objects.create(
            event=event,
            start_time=event_start_time.datetime,
            end_time=event_end_time.datetime,
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
    result = True
    return HttpResponse(result)

#2.1 �������� � ��������� �������
