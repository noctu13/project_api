from rest_framework import viewsets
from django.shortcuts import render
from .models import Event, Polyline, Point
from .serializers import EventSerializer, PolylineSerializer, PointSerializer

#1 Создать get запрос на полигоны с корональными дырами
class EventViewSet(viewsets.ModelViewSet):
    queryset = Event.objects.all()
    serializer_class = EventSerializer

class PointViewSet(viewsets.ModelViewSet):
    queryset = Point.objects.all()
    serializer_class = PointSerializer

class PolylineViewSet(viewsets.ModelViewSet):
    queryset = Polyline.objects.all()
    serializer_class = PolylineSerializer

#2 get запрос на обновление данных с HEK (повышенные привелегии)
#2.1 создание и обработка токенов
