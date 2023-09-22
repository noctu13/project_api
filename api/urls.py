from django.urls import path, include

from .views import EventsView


urlpatterns = [
    path('events/<str:short_type>/<int:year>/<int:month>/<int:day>/', EventsView.as_view()),
    path('events/<str:short_type>/', EventsView.as_view()),
]