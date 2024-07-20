from django.urls import path

from .views import EventsView


urlpatterns = [
    path('events/<str:short_type>/<int:year>/<int:month>/<int:day>/', EventsView.as_view()),
    
]