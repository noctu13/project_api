from rest_framework import serializers
from rest_framework.settings import api_settings

from .models import (HEKCoronalHole, HEKCoronalHoleContourPoint,
    CoronalHole, CoronalHoleContour, CoronalHoleContourPoint, CoronalHolePoint,
    MagneticLineSet, MagneticLine, MagneticLinePoint)
from .common import *


class HEKCoronalHoleContourPointSerializer(serializers.ModelSerializer):

    class Meta:
        model = HEKCoronalHoleContourPoint
        fields = '__all__'

class HEKCoronalHoleSerializer(serializers.ModelSerializer):
    contour = HEKCoronalHoleContourPointSerializer(many=True)
    
    class Meta:
        model = HEKCoronalHole
        fields = '__all__'
    
    def to_representation(self, instance):
        data = super().to_representation(instance)
        point_list = [str(point) for point in instance.contour.all()]
        data['contour'] = point_list
        return data

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
        return data

class CoronalHoleSerializer(serializers.ModelSerializer):
    contour = CoronalHoleContourSerializer(many=True)
    
    class Meta:
        model = CoronalHole
        fields = '__all__'

class CoronalHolePointSerializer(serializers.ModelSerializer):

    class Meta:
        model = CoronalHolePoint
        fields = '__all__'

class MagneticLineSetSerializer(serializers.ModelSerializer):

    class Meta:
        model = MagneticLineSet
        fields = '__all__'

class MagneticLineSerializer(serializers.ModelSerializer):

    class Meta:
        model = MagneticLine
        fields = '__all__'

class MagneticLinePointSerializer(serializers.ModelSerializer):

    class Meta:
        model = MagneticLinePoint
        fields = '__all__'