from django.db import models

from .common import *


class HEKCoronalHole(models.Model):
    spec_id = models.CharField(max_length=24)
    sol = models.CharField(max_length=34)
    start_time = models.DateTimeField()
    end_time = models.DateTimeField()


class HEKCoronalHoleContourPoint(models.Model):
    lon = models.FloatField()
    lat = models.FloatField()
    ch = models.ForeignKey(
        HEKCoronalHole, 
        on_delete=models.CASCADE,
        related_name='contour',
    )

    def __str__(self):
        return f'{self.lon} {self.lat}'


class CoronalHole(models.Model):
    start_time = models.DateTimeField()
    type = models.CharField(max_length=3, choices=ch_dict.items())
    location = models.CharField(max_length=14, null=True)
    sol = models.CharField(max_length=34, null=True)
    area = models.FloatField(null=True)
    mag_flux = models.FloatField(null=True)
    avg_flux = models.FloatField(null=True)
    max_flux = models.FloatField(null=True)


class CoronalHolePoint(models.Model):
    lon = models.FloatField()
    lat = models.FloatField()
    Br = models.FloatField()
    ch = models.ForeignKey(
        CoronalHole,
        on_delete=models.CASCADE,
        related_name='points',
    )

    def __str__(self):
        return f'{self.lon} {self.lat} {self.B}'


class CoronalHoleContour(models.Model):
    ch = models.ForeignKey(
        CoronalHole,
        on_delete=models.CASCADE,
        related_name='contour',
    )

class CoronalHoleContourPoint(models.Model):
    lon = models.FloatField()
    lat = models.FloatField()
    contour = models.ForeignKey(
        CoronalHoleContour,
        on_delete=models.CASCADE,
        related_name='points',
    )

    def __str__(self):
        return f'{self.lon} {self.lat}'


class MagneticLineSet(models.Model):
    start_time = models.DateTimeField()
    type = models.CharField(max_length=4, choices=ml_dict.items())


class MagneticLine(models.Model):
    polarity = models.BooleanField(null=True)
    lineset = models.ForeignKey(
        MagneticLineSet,
        on_delete=models.CASCADE,
        related_name='lines',
    )


class MagneticLinePoint(models.Model):
    lon = models.FloatField()
    lat = models.FloatField()
    r = models.FloatField()
    line = models.ForeignKey(
        MagneticLine,
        on_delete=models.CASCADE,
        related_name='points',
    )
    
    def __str__(self):
        return f'{self.lon} {self.lat} {self.r}'