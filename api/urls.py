from django.urls import path, include
from rest_framework.routers import SimpleRouter

from .views import CoronalHolesView, PFSSMagneticLinesView


urlpatterns = [
    path('events/CH/', CoronalHolesView.as_view()),
    path('events/ML/', PFSSMagneticLinesView.as_view())
]