from django.urls import path, include
from rest_framework.routers import DefaultRouter

from api.models import Event, Polyline, Point

from .views import EventViewSet, PolylineViewSet, PointViewSet


router = DefaultRouter()
router.register('events', EventViewSet, basename='events')
router.register('polylines', PolylineViewSet, basename='polylines')
router.register('points', PointViewSet, basename='points')

urlpatterns = [
    path('', include(router.urls)),
]