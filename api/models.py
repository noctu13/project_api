from django.db import models


class Event(models.Model):
    type = models.CharField(max_length=20)
    start_time = models.DateTimeField(null=True)
    end_time = models.DateTimeField(null=True)
    spec_id = models.CharField(max_length=36, null=True, unique=True)


class Polyline(models.Model):
    event = models.ForeignKey(
        Event, 
        on_delete=models.CASCADE,
        related_name='polyline',
    )
    start_time = models.DateTimeField()
    end_time = models.DateTimeField()
    polarity = models.BooleanField(null=True)


class Point(models.Model):
    phi = models.FloatField()
    theta = models.FloatField()
    r = models.FloatField(null=True)
    polyline = models.ForeignKey(
        Polyline, 
        on_delete=models.CASCADE,
        related_name='points',
    )
    
    def __str__(self):
        tail = f' {self.r}' if self.r else ''
        return f'{self.phi} {self.theta}{tail}'
