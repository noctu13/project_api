from rest_framework import viewsets
from django.shortcuts import render
from .models import Event, Polyline, Point
from .serializers import EventSerializer, PolylineSerializer, PointSerializer

#1 ������� get ������ �� �������� � ������������ ������
class EventViewSet(viewsets.ModelViewSet):
    queryset = Event.objects.all()
    serializer_class = EventSerializer

class PointViewSet(viewsets.ModelViewSet):
    queryset = Point.objects.all()
    serializer_class = PointSerializer

class PolylineViewSet(viewsets.ModelViewSet):
    queryset = Polyline.objects.all()
    serializer_class = PolylineSerializer

#2 get ������ �� ���������� ������ � HEK (���������� ����������)
#2.1 �������� � ��������� �������
