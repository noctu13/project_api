from rest_framework import serializers
from rest_framework.settings import api_settings

from .models import (CoronalHole, CoronalHolePoint, 
    CoronalHoleContour, CoronalHoleContourPoint,
    MagneticLineSet, MagneticLine, MagneticLinePoint)
from .common import *


class CoronalHoleContourPointSerializer(serializers.ModelSerializer):

    class Meta:
        model = CoronalHoleContourPoint
        fields = '__all__'

class CoronalHoleContourSerializer(serializers.ModelSerializer):
    points = CoronalHoleContourPointSerializer(many=True)

    class Meta:
        model = CoronalHoleContour
        fields = '__all__'
    
    def to_representation(self, instance):
        data = super().to_representation(instance)
        point_list = [str(point) for point in instance.points.all()]
        data['points'] = point_list
        del data['ch']
        return data

class CoronalHolePointSerializer(serializers.ModelSerializer):

    class Meta:
        model = CoronalHolePoint
        fields = '__all__'

class CoronalHoleSerializer(serializers.ModelSerializer):
    contour = CoronalHoleContourSerializer(many=True)
    points = CoronalHolePointSerializer(many=True)
    
    class Meta:
        model = CoronalHole
        fields = '__all__'
    
    def to_representation(self, instance):
        data = super().to_representation(instance)
        point_list = [str(point) for point in instance.points.all()]
        data['points'] = point_list
        data['type'] = ch_dict[instance.s_type]
        del data['s_type']
        return data

class MagneticLinePointSerializer(serializers.ModelSerializer):

    class Meta:
        model = MagneticLinePoint
        fields = '__all__'

class MagneticLineSerializer(serializers.ModelSerializer):
    points = MagneticLinePointSerializer(many=True)

    class Meta:
        model = MagneticLine
        fields = '__all__'
    
    def to_representation(self, instance):
        data = super().to_representation(instance)
        point_list = [str(point) for point in instance.points.all()]
        data['points'] = point_list
        del data['lineset']
        return data

class MagneticLineSetSerializer(serializers.ModelSerializer):
    lines = MagneticLineSerializer(many=True)

    class Meta:
        model = MagneticLineSet
        fields = '__all__'
    
    def to_representation(self, instance):
        data = super().to_representation(instance)
        data['type'] = ml_dict[instance.s_type]
        del data['s_type']
        return data