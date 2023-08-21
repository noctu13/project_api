from rest_framework import viewsets
from django.shortcuts import render
from .models import Events
from .serializers import EventsSerializer

# Create your views here.
#1 Создать get запрос на полигоны с корональными дырами
class EventsViewSet(viewsets.ModelViewSet):
    queryset = Events.objects.all()
    serializer_class = EventsSerializer


#2 get запрос на обновление данных с HEK (повышенные привелегии)
#2.1 создание и обработка токенов
