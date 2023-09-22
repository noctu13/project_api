from django.urls import path, include
from rest_framework.routers import SimpleRouter

from .views import EventViewSet


router = SimpleRouter()
router.register('events', EventViewSet)

urlpatterns = [
    path('', include(router.urls)),
]