from django.db import models


class Event(models.Model):
    type = models.CharField(max_length=20)
    start_time = models.DateTimeField()
    end_time = models.DateTimeField()


class Polyline(models.Model):
    event = models.ForeignKey(
        Event, 
        on_delete=models.CASCADE,
        related_name='polyline',
    )
    start_time = models.DateTimeField()
    end_time = models.DateTimeField()
    


class Point(models.Model):
    theta = models.FloatField()
    phi = models.FloatField()
    r = models.FloatField(default=1.0)
    polyline = models.ForeignKey(
        Polyline, 
        on_delete=models.CASCADE,
        related_name='points',
    )
