from django.db import models

class Events(model.Model):
    name = models.CharField(max_length=50)
    date = models.DateTimeField()
    type = models.CharField(max_length=20)
    polynome = models.JSONField()