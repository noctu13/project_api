from rest_framework import serializers
from .models import Events

class EventsSerializer(serializers.Serializers):
    class Meta:
        model = Events
        fields = ['name', 'date', 'type', 'polynome']