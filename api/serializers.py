from rest_framework import serializers
from .models import Event, Polyline, Point


class PointSerializer(serializers.ModelSerializer):
    
    class Meta:
        model = Point
        fields = ("theta", "phi", "r")
        

class PolylineSerializer(serializers.ModelSerializer):
    points = PointSerializer(
        many=True,
        read_only=True,
    )
    
    class Meta:
        model = Polyline
        fields = ("id", "start_time", "end_time", "points")


class EventSerializer(serializers.ModelSerializer):
    polyline = PolylineSerializer(
        many=True,
        read_only=True,
    )
    
    class Meta:
        model = Event
        fields = ("id", "type", "start_time", "end_time", "polyline")