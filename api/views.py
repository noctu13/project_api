from datetime import datetime, timedelta, timezone

from rest_framework import generics
from django.db.models import Q

from .common import *
from .models import HEKCoronalHole, CoronalHole, MagneticLineSet
from .serializers import (
    HEKCoronalHoleSerializer, CoronalHoleSerializer, MagneticLineSetSerializer)


class EventsView(generics.ListAPIView):

    def get_queryset(self):
        short_type = self.kwargs.get('short_type')
        year = self.kwargs.get('year')
        month = self.kwargs.get('month')
        day = self.kwargs.get('day')
        if short_type == 'HCH':
            events_query = HEKCoronalHole.objects.all()
        elif short_type in ch_dict:
            events_query = CoronalHole.objects.filter(type=short_type)
        elif short_type in ml_dict:
            events_query = MagneticLineSet.objects.filter(type=short_type)
        if year and month and day:
            query_day = datetime(year, month, day, 0, 0, 0,
                tzinfo=timezone.utc)
            query1 = Q(start_time__gt=query_day)
            query2 = Q(start_time__lte=query_day + timedelta(days=1))
            queryset = events_query.filter(query1 & query2).order_by(
                'start_time')
        return queryset
    
    def get_serializer_class(self):
        short_type = self.kwargs.get('short_type')
        if short_type == 'HCH':
            serializer = HEKCoronalHoleSerializer
        elif short_type in ch_dict:
            serializer = CoronalHoleSerializer
        elif short_type in ml_dict:
            serializer = MagneticLineSetSerializer
        return serializer


# написать эндпоинт - фильтрованная выдача media(ph/ss/sw date/cr stop/gong)?