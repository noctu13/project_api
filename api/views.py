from rest_framework import viewsets
from django.shortcuts import render

from .models import Event, Polyline, Point
from .serializers import EventSerializer, PolylineSerializer, PointSerializer
from .permissions import IsAdminOrReadOnly


#1 ������� get ������ �� �������� � ������������ ������
class EventViewSet(viewsets.ModelViewSet):
    queryset = Event.objects.all()
    serializer_class = EventSerializer
    permission_classes = (IsAdminOrReadOnly,)

#2 get ������ �� ���������� ������ � HEK (���������� ����������)
#2.1 �������� � ��������� �������
