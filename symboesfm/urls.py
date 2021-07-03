
from django.contrib import admin
from django.urls import path
from symboesfm import views as home_view
from django.conf import settings
from django.conf.urls.static import static

from integracion import views as integracion_views

urlpatterns = [
    path('admin/', admin.site.urls),
    path('creditos/', home_view.creditos, name = 'creditos'),
    path('', home_view.home, name = 'home'),
    path('construccion/', home_view.construccion, name = 'construccion'),

    #Integracion
    path('integracion/view/', integracion_views.view, name = 'view'),
    path('integracion/simple/', integracion_views.simple, name =  "simple"),
    path('integracion/doble/', integracion_views.doble, name =  "doble"),
    path('integracion/submit/', integracion_views.submit, name = 'submit'),
    path('integracion/extrapolacion/', integracion_views.extrapolacion, name = 'extrapolacion'),
    path('integracion/indefinida', integracion_views.indefinida, name = 'indefinida')

]+ static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
