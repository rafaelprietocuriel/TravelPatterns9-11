# RoutesSeptember11
Data on the routes followed by the 19 people on the September 11 attacks.

To map the spatial patterns of the 9/11 hijackers, we collected openly accessible data from the 9/11 Commission Report, the 9/11 Travel Staff Report, the U.S. Congress and a redacted document from the FBI (2007) that provides a detailed chronology of events for hijackers and associates. In each document, we searched for the name of the 19 hijackers and reported their location, geographic coordinates, and date of each travel since they left their country

When the exact location of an event was unknown, we used the center of gravity of the region or country.  When the exact date was not available, we assumed that it took place at an equal number of days from two known events. When an event took place “early” in a month according to the sources, we used the 5th calendar day, and when it took place “late” we used the 25th. In total, our data contains detailed information about 364 events and 83 unique locations.

Some of the locations were merged when they belonged to the same metropolitan area, resulting in 48 locations.

For flights with a layover, the person is assigned to the metropolitan area where they spent the night. In total we have 231 movements between metropolitan areas.

# Data structure
The data is structured in four tables:
- Attributes
Contains information about each of the 19 individuals, including name, role, nationality and others. The column "n" is the ID for other tables.

- Routes
Table sorted chronologically and by actor. Each row contains a location where person X is known to be at some date. We assume that a journey between any two consecutive rows was perfomed by that person. The column date corresponds to the first moment in which the person was known to be there. It is formated in standard consecutive dates, so the number 36679 corresponds to June 2, 2000. The column "POI" is used in the table Locations. The column "Actor" corresponds to the column "n" in the Attributes table.

- Locations
Table with the 48 metropolitan areas that were visited by hijackers. They are identified by the column "POI". Contain the coordinates and the number of times that it was the begining or ending of a journey.

- SocialcontactMatrix
Corresponds to a 19 X 19 table, where cell i,j = 1 if i and j are connected in the social contact network. The table was obtained in the UCINET software datasets. Columns and rows are labeled with the name that corresponds to the column "n" on the Attributes table.
