from datetime import datetime, timedelta, timezone

from rest_framework import generics
from django.db.models import Prefetch, Q

from .common import *
from .models import Event, Polyline, Point
from .serializers import EventSerializer


class EventsView(generics.ListAPIView):
    serializer_class = EventSerializer

    def get_queryset(self):
        short_type = short_type_dict[self.kwargs['short_type']]
        events_query = Event.objects.filter(type=short_type)
        year = self.kwargs.get('year')
        month = self.kwargs.get('month')
        day = self.kwargs.get('day')
        if year and month and day:
            query_day = datetime(year, month, day, 0, 0, 0,
                tzinfo=timezone.utc)
            query1 = Q(start_time__lte=query_day + timedelta(days=1))
            query2 = Q(end_time__gt=query_day)
            into_query_day = Polyline.objects.filter(
                query1 & query2).order_by('start_time')
            queryset = events_query.filter(query1 & query2).order_by(
                'start_time').prefetch_related(
                Prefetch(
                    'polyline',
                    queryset=into_query_day,
                    to_attr='into_query_day'
                )
            )
            if not queryset and short_type == short_type_dict['PML']:
                queryset = events_query.latest('start_time')
                queryset.into_query_day = queryset.polyline
                queryset = [queryset]
        return queryset # много элементов > 150k - ошибка

# написать эндпоинт - фильтрованная выдача media(ph/ss/sw date/cr stop/gong)?