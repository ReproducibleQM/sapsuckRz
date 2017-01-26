//Project// Aphid flight patterns in the US midwest

//Data source// Aphids are tiny soft-bodied herbivores with very interesting life cycles- 
many species alternate between host species at different times of the year, and will also migrate between 
hosts of the same species due to a variety of environmental or densitiy dependant triggers. 
The Midwestern Aphid Suction trap network was established in 2004, originally as a surveilance tool for the recent 
invader and economic pest soybean aphid. Its 42 sites stretch over a 1000 km x 1000 KM area in the US midwest, from Manhattan, 
KS in the south to Crookston, MN in the north, and Brookings, SD in the west to Lake Erie in the east. 
The network provides weekly data on aphid captures May-September. Aphids are idenified by species and sex, and tallied for each location 
for each week. Approximately 120 species have been captured by this network, with ~25 species common enough to support a meaningful 
analysis of some kind. The survey was initiated by David Voegtlin of the Illinois Natural History Survey, and since his retirement, 
has been maintained by Doris Lagos, an (amazing) aphid systematist. Dr. Christie Bahlai is their database manager, and holds 
the database at the KBS data repository.

//Contacts// Christie Bahlai, Doris Lagos

//Dataset (full raw data)// https://dl.dropboxusercontent.com/u/98197254/Mastersuction_may29_2014.csv

//Variables// 
V1= Line ID.
Year= The year that the specimens were collected
State= The U.S. state that the specimens were collected. (IA=Iowa, IL= Illinois, IN=Indiana, KS= Kansas, KY= Kentucky, MI= Michigan, 
MN= Minnesota, MO= Missouri, SD= South Dakota, WI= Wisconsin)
Site= Specific location of turbine trap where specimen were collected. See Appendix A for exact Lat/Long coordinates.
Date= The day of the year (out of 365, 366 on leapyears) on which the speciam was collected. 
Sex= The sex of the collected speciman. (M=Male, F=Female)
Captures= The number of individuals caught at that site on that date.
Aphid.Species= The genus and species of the captured specimens. See Appendix B for list of pests of economic significance.

//Appendix A: Trap Locations// Data is tab seperated. 
SITEID  LAT LONG
ARVA	43.08	-81.20
BEAN & BEAT	43.30	-84.13
BROOKINGS	44.31	-96.78
BROWNSTOWN	38.95	-88.96
COLUMBIA	38.94	-92.32
CROOKSTON	47.81	-96.48
DEKALB	41.84	-88.85
DIXON SPRING 	37.43	-88.67
DPAC	40.25	-85.15
EAST LANSING	42.72	-84.46
ELORA	43.64	-80.41
FREEPORT	42.28	-89.70
HANCOCK	44.11	-89.54
HARROW	42.03	-82.90
HURON	43.32	-81.50
KELLOGG	42.32	-85.38
LAMBERTON	44.24	-95.32
LANCASTER	42.83	-90.79
LEXINGTON	38.09	-84.53
MANHATTAN	39.14	-96.64
MCNAY	40.97	-93.42
METAMORA	40.77	-89.34
MONMOUTH	40.94	-90.72
MONROE	41.92	-83.39
MORRISIL	41.37	-88.43
MORRISMN	45.59	-95.90
NASHUA	42.93	-92.57
NEPAC	41.10	-85.40
OCEANA	43.56	-86.38
PERRY	39.81	-90.82
PIONEER	44.76	-91.56
PIT	40.45	-86.93
PORTAGEVILLE	36.41	-89.70
PPAC	41.44	-86.93
PRINCETON	37.10	-87.86
RIDGETOWN	42.45	-81.88
ROSEMOUNT	44.72	-93.10
SEPAC	39.04	-85.53
SEYMOUR	44.52	-88.33
SUTHERLAND	42.92	-95.54
URBANA	40.10	-88.23
WALWORTH	42.53	-88.60

//Appendix B: Pest Species// Data is tab seperated. 
Acyrthosiphon pisum "pea aphid"
Aphis craccivora  "black legume aphid"
Aphis glycines  "soybean aphid"
Aphis gossypii  "cotton-melon aphid"
Aphis helianthi "sunflower or dogwood aphid"
Aphis nasturtii "buckthorn-potato aphid"
Aphis spiraecola  "spiraea aphid"
Lipaphis pseudobrassicae  "turnip aphid"
Macrosiphum euphorbiae  "potato aphid"
Myzus persicae  "peach potato aphid"
Rhopalosiphum insertum  "apple grass aphid"
Rhopalosiphum maidis  "corn leaf aphid"
Rhopalosiphum padi  "bird cherry-oat aphid"
Rhopalosiphum rufiabdominale  "rice root aphid"
Schizaphis graminum "greenbug"
Sitobion avenae "english grain aphid"
Therioaphis trifolii  "spotted alfalfa aphid"
