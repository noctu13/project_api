from rest_framework import viewsets
from django.shortcuts import render

from .models import Event, Polyline, Point
from .serializers import EventSerializer, PolylineSerializer, PointSerializer
from .permissions import IsAdminOrReadOnly


#1 Создать get запрос на полигоны с корональными дырами
class EventViewSet(viewsets.ModelViewSet):
    queryset = Event.objects.all()
    serializer_class = EventSerializer
    permission_classes = (IsAdminOrReadOnly,)

#2 get запрос на обновление данных с HEK (повышенные привелегии)
#2.1 создание и обработка токенов
