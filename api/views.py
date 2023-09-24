from rest_framework import generics
from django.shortcuts import render
from datetime import date, timedelta

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

#2 get запрос на обновление данных с HEK (повышенные привелегии)
#2.1 создание и обработка токенов
