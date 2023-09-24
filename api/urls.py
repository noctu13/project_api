from django.urls import path, include

from .views import EventsView, load_HEK_CH


urlpatterns = [
    path('events/<str:short_type>/<int:year>/<int:month>/<int:day>/', EventsView.as_view()),
    path('events/<str:short_type>/', EventsView.as_view()),
    path('control/load', load_HEK_CH),
]