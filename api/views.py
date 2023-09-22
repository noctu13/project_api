from rest_framework import generics
from django.shortcuts import render

from .models import Event, Polyline, Point
from .serializers import EventSerializer, PolylineSerializer, PointSerializer
from .permissions import IsAdminOrReadOnly


class CoronalHolesView(generics.ListAPIView):
    queryset = Event.objects.filter(type='Coronal Hole')
    serializer_class = EventSerializer

class PFSSMagneticLinesView(generics.ListAPIView):
    queryset = Event.objects.filter(type='PFSS magnetic line')
    serializer_class = EventSerializer

#2 get запрос на обновление данных с HEK (повышенные привелегии)
#2.1 создание и обработка токенов
