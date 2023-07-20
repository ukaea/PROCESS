from django.conf.urls import patterns, include, url
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
from django.contrib import admin

urlpatterns = (
    patterns(
        "",
        url(r"^$", "ui.views.index"),
        url(r"^admin/", include(admin.site.urls)),
    )
    + staticfiles_urlpatterns()
)
