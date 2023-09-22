from rest_framework import serializers
from .models import Event, Polyline, Point


class EventSerializer(serializers.ModelSerializer):
    polyline = serializers.PrimaryKeyRelatedField(
        queryset=Polyline.objects.all(),
        many=True,
    )

    class Meta:
        model = Event
        fields = ['type', 'start_time', 'end_time', 'polyline']


class PolylineSerializer(serializers.ModelSerializer):
    points = serializers.PrimaryKeyRelatedField(
        queryset=Point.objects.all(),
        many=True,
    )

    class Meta:
        model = Polyline
        fields = ['event', 'start_time', 'end_time', 'points']


class PointSerializer(serializers.ModelSerializer):

    class Meta:
        model = Point
        fields = '__all__'