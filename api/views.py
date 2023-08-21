from rest_framework import viewsets
from django.shortcuts import render
from .models import Events
from .serializers import EventsSerializer

# Create your views here.
#1 ������� get ������ �� �������� � ������������ ������
class EventsViewSet(viewsets.ModelViewSet):
    queryset = Events.objects.all()
    serializer_class = EventsSerializer


#2 get ������ �� ���������� ������ � HEK (���������� ����������)
#2.1 �������� � ��������� �������
