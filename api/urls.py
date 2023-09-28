from django.urls import path, include

from .views import EventsView

urlpatterns = [
    path('', EventsView.as_view()),
]