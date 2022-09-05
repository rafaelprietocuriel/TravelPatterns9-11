# Travel Patterns 9-11
Data on the routes followed by the 19 Al Qaeda operatives responsible for the September 11, 2001 attacks against the United States.


To map the spatial patterns of the 9/11 hijackers, we collected open data from the 9/11 Commission Report, the 9/11 Travel Staff Report, the U.S. Congress and a redacted document from the FBI (2007) that provides a detailed chronology of events for hijackers and associates. In each document, we searched for the name of the 19 hijackers and reported their location, geographic coordinates, and date of each travel since they left their country.

When the exact location of an event was unknown, we used the center of gravity of the region or country. When the exact date was unavailable, we assumed it took place at an equal number of days from two known events. When an event took place “early” in a month, according to the sources, we used the 5th calendar day, and when it took place “late” we used the 25th. The data contains detailed information about 364 events and 83 unique locations.

Some locations were merged when they belonged to the same metropolitan area, resulting in 48 locations.

For flights with a layover, the person is assigned to the metropolitan area where they spent the night. In total, we have 231 movements between metropolitan areas.

# Data structure
The data is structured in four tables:
- Attributes
 
Contain information about each of the 19 individuals, including name, role, nationality and others. The column "n" is the ID for other tables.


- Routes

Table sorted chronologically and by the actor. Each row contains a location where person X is known to be at some date. We assume that that person performed a journey between any two consecutive rows. The column date corresponds to the first moment in which the person was known to be there. It is formatted in consecutive standard dates, so the number 36679 corresponds to June 2, 2000. The column "POI" is used in the table Locations. The column "Actor" corresponds to the column "n" in the Attributes table.


- Locations

Table with the 48 metropolitan areas that were visited by hijackers. They are identified by the column "POI". Contain the coordinates and the number of times that it was the beginning or end of a journey.


- SocialcontactMatrix

Corresponds to a 19 X 19 table, where cell i,j = 1 if i and j are connected in the social contact network. The data was collected by Krebs, V.E. (2002). Mapping networks of terrorist cells. Connections 24.3: 43-52 and is openly available on the UCINET website (https://sites.google.com/site/ucinetsoftware/datasets/covert-networks/911-hijackers). Columns and rows are labelled with the name corresponding to the column "n" on the Attributes table. The data was symmetrized.


# How to cite
If you use this database, please cite it as below: 

Walther, O., Prieto-Curiel R., Padron, J., Scheuer, J. 2022. Mapping the Travel Geography of the 9/11 Terrorist Network. Available at SSRN: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4206577
