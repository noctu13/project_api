from rest_framework import generics
from django.shortcuts import render
from datetime import date, timedelta

from .models import Event, Polyline, Point
from .serializers import EventSerializer, PolylineSerializer, PointSerializer
from .permissions import IsAdminOrReadOnly


class EventsView(generics.ListAPIView):
    serializer_class = EventSerializer
    short_type_dict = {
        'CH': 'Coronal Hole',
        'PML': 'PFSS Magnetic Line',
    }
    
    def get_queryset(self):
        short_type = self.short_type_dict[self.kwargs['short_type']]
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
#2.1 �������� � ��������� �������
