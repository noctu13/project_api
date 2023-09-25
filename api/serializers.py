from rest_framework import ISO_8601
from rest_framework import serializers
from rest_framework.settings import api_settings

from .models import Event, Polyline, Point
from .common import *


class PointSerializer(serializers.ModelSerializer):
    
    class Meta:
        model = Point
        exclude = ('id', 'polyline', )


class PolylineSerializer(serializers.ModelSerializer):
    start_time = serializers.DateTimeField()
    end_time = serializers.DateTimeField()
    
    class Meta:
        model = Polyline
        fields = ('id', 'start_time', 'end_time','polarity', 'points')
        
    def to_representation(self, instance):
        data = super().to_representation(instance)
        point_list = []
        for point in instance.points.all():
            point_list.append(str(point))
        data['points'] = point_list
        if instance.event.type == short_type_dict['CH']:
            data.pop('polarity')
        return data


class EventSerializer(serializers.ModelSerializer):
    start_time = serializers.DateTimeField()
    end_time = serializers.DateTimeField()
    polyline = PolylineSerializer(
        source='into_query_day',
        many=True,
        read_only=True,
    )
    
    class Meta:
        model = Event
        fields = ('id', 'type', 'start_time', 'end_time', 'polyline')