# Project SunGazers

The SunGazers project is a service for providing data related to solar activity. The focus is on space weather, and specifically, it is historical and forecast data on solar wind and coronal mass ejection data. 
The source data - photospheric magnetogram are taken from the Kislovodsk Mountain Astronomical Station
of the Pulkovo observatory [Kislovodsk MAS PO](http://solarstation.ru/) site.
The api of the SunGazers project is maintained in this repository.

At this stage, the database looks like
![database_schema](./readme_static/database_schema.png?raw=true)
and stores 2 types of events: coronal holes (CH) and coronal magnetic field lines (PML).
The api implements issuance of certain events on the date
    https://{sitename}:{port}/api/v1/events/${event}/${year}/${month}/${day}/
