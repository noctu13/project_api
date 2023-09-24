from rest_framework import ISO_8601
from rest_framework import serializers
from rest_framework.settings import api_settings

from .models import Event, Polyline, Point
from .common import *


class CustomDateTimeField(serializers.DateTimeField):

    def to_representation(self, value):
        if not value:
            return None
        output_format = getattr(self, 'format', api_settings.DATETIME_FORMAT)
        if output_format is None or isinstance(value, str):
            return value
        value = self.enforce_timezone(value)
        if output_format.lower() == ISO_8601:
            value = value.isoformat()
            return value
        return value.strftime(output_format)

class PointSerializer(serializers.ModelSerializer):
    
    class Meta:
        model = Point
        exclude = ('id', 'polyline', )
    
    def to_internal_value(self, data):
        lat = data.pop('lat', None)
        lon = data.pop('lon', None)
        if lat or lon:
            data['phi'] = lat
            data['theta'] = lon
        return super().to_internal_value(data)
    
    def to_representation(self, instance):
        data = super().to_representation(instance)
        event_type = instance.polyline.event.type
        if event_type == short_type_dict['CH']:
            data['lat'] = data.pop('phi')
            data['lon'] = data.pop('theta')
            data.pop('r')
        return data
        

class PolylineSerializer(serializers.ModelSerializer):
    start_time = CustomDateTimeField()
    end_time = CustomDateTimeField()
    points = PointSerializer(
        many=True,
        read_only=True,
    )
    
    class Meta:
        model = Polyline
        fields = ("id", "start_time", "end_time", "points")


class EventSerializer(serializers.ModelSerializer):
    start_time = CustomDateTimeField()
    end_time = CustomDateTimeField()
    polyline = PolylineSerializer(
        many=True,
        read_only=True,
    )
    
    class Meta:
        model = Event
        fields = ("id", "type", "start_time", "end_time", "polyline")