from rest_framework import ISO_8601
from rest_framework import serializers
from rest_framework.settings import api_settings


from .models import Event, Polyline, Point

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
        fields = ("theta", "phi", "r")
        

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