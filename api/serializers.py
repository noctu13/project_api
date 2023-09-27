from rest_framework import ISO_8601
from rest_framework import serializers
from rest_framework.settings import api_settings

from .models import Event, Polyline, Point
from .common import *


class PointSerializer(serializers.ModelSerializer):
    
    class Meta:
        model = Point
        exclude = ('id', 'polyline',)


class PolylineSerializer(serializers.ModelSerializer):
    start_time = serializers.DateTimeField()
    end_time = serializers.DateTimeField()
    
    class Meta:
        model = Polyline
        fields = ('id', 'start_time', 'end_time','polarity',)
        
    def to_representation(self, instance):
        data = super().to_representation(instance)
        event_dict = {}
        for key in short_type_dict.keys():
             event_dict[key] = instance.event.type == short_type_dict[key]
        point_list = []
        for point in instance.points.all():
            if event_dict['CH']:
                point_line = f'{point.phi:.1f} {point.theta:.1f}'
            elif event_dict['PML']:
                point_line = f'{point.phi:.3f}'\
                    f' {point.theta:.3f} {point.r:.3f}'
            else: str(point)
            point_list.append(point_line)
        data['points'] = point_list
        if event_dict['CH']:
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
        fields = ('id', 'type', 'start_time', 'end_time', 'polyline',)