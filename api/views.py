from datetime import date, datetime, timedelta

from rest_framework import generics
from django.db.models import Q

from .common import *
from .utils import zero_time
from .models import CoronalHole, MagneticLineSet
from .serializers import (
    CoronalHoleSerializer, MagneticLineSetSerializer)


class EventsView(generics.ListAPIView):

    def get_queryset(self):
        short_type = self.kwargs.get('short_type')
        year = self.kwargs.get('year')
        month = self.kwargs.get('month')
        day = self.kwargs.get('day')
        if short_type in ch_dict:
            events_query = CoronalHole.objects.filter(s_type=short_type)
        elif short_type in ml_dict:
            events_query = MagneticLineSet.objects.filter(s_type=short_type)
        if year and month and day:
            query_day = datetime.combine(date(year, month, day), zero_time)
            queryset = None
            while not queryset:
                query1 = Q(start_time__gte=query_day)
                query2 = Q(start_time__lt=query_day + timedelta(days=1))
                queryset = events_query.filter(query1 & query2).order_by('start_time')
                query_day -= timedelta(days=1)
        return queryset
    
    def get_serializer_class(self):
        short_type = self.kwargs.get('short_type')
        if short_type in ch_dict:
            serializer = CoronalHoleSerializer
        elif short_type in ml_dict:
            serializer = MagneticLineSetSerializer
        return serializer
    
    def get_serializer_context(self):
        context = super().get_serializer_context()
        param = self.request.query_params.get('param')
        if param:
            context.update({param: True})
        return context


# написать эндпоинт - фильтрованная выдача media(ph/ss/sw date/cr stop/gong)?
