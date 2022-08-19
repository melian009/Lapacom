#Summary of DB_Phorcus :
Using the BD_LIMPETS_MAD_1996-2018.csv file as reference, the information in the databases (.xlxs) in the Lapacom/Data/ToCheck folder has been restructured.
Common database fixes:
- 	Individual_Number renumbering from 1 to n in all files
-	Format change for Lat and Long to DMS (degrees, minutes, and seconds) in all files. Format example: Lat 6° 57' 06.5"N and Long 25° 06' 20.3"Wº.
-	Write correction for a level of Protective_regime: from "full acess" to "full access". This correction was made in all the files, to avoid problems in subsequent graph and table representations.
-	Converting .xlsx files to .csv
The non-modified databases were not in the same format as the reference file. These files are not included in this folder.
Files in “Lapacom/Data/ToCheck/”:
Folder Canaries: Non-modified

BD_Phorcus_Azores_2018:

-	Individual_number: n = 70
-	Species: _Phorcus sauciatus_
-	Locality: Azores
-	Year: 2018
-	Month: NaN
-	Total_length_mm: numeric
-	Weigth_g: numeric
-	Sampling_site: 1 Site
-	Lat: DMS 
-	Long: DMS 
-	Protective_regime: Full access
-	Proximity_human_settlements: Far
 
BD_Phorcus_Canaries_2017-2018
 
-	Individual_number: n = 1374
-	Species: _Phorcus sauciatus_ and _Phorcus atratus_
-	Locality: Canary Island
-	Year: 2017 - 2018
-	Month: numeric
-	Total_length_mm: numeric 
-	Weigth_g: numeric.
-	Sampling_site: 5 sites
-	Lat: DMS 
-	Long: DMS
-	Protective_regime: Full access/MPS
-	Proximity_human_settlements: Far/Near

BD_Phorcus_CapeVerde_2017
 
-	Individual_number: n = 126
-	Species: _Phorcus mariae_
-	Locality: Cape Verde
-	Year: 2017 
-	Month: 1 month (numeric)
-	Data: dd/mm/aaaa
-	Total_length_mm: numeric
-	Weigth_g: numeric
-	Sampling_site: 1 site
-	Lat: DMS 
-	Long: DMS
-	Protective_regime: Full access
-	Proximity_human_settlements: Near

BD_Phorcus_Madeira_2016-2018

-	Individual_number: n = 7069
-	Species: _Phorcus sauciatus_
-	Locality: Madeira archipelago
-	Year: 2016 - 2018
-	Month: numeric
-	Data: dd/mm/aaaa
- 	Total_length_mm: numeric
-	Weigth_g: numeric
-	Sampling_site: 8 sites 
-	Lat: DMS 
-	Long: DMS
-	Protective_regime: Full access/MPA
-	Proximity_human_settlements: Control/Far/Near
-	Accessibility: North/South/NaN

In the original database, 274 rows had incorrect date information (N_o=7343). These files were dated for Month: 1 and Year: 1900. As they did not present "Date" values that would verify these values, the samples from these rows were eliminated from the final .csv file.

BD_Phorcus_mainlandPortugal_2017-2018

-	Individual_number: n = 7343
-	Species: _Phorcus sauciatus_
-	Locality: mainland Portugal
-	Year: 2017 - 2018
-	Month: numeric
-	Data: dd/mm/aaaa
-	Total_length_mm: numeric
-	Weigth_g: numeric
-	Water_salinity: numeric/NaN
-	Water_temperature: numeric/NaN
-	Sampling_site: 9 sites 
-	Lat: DMS 
-	Long: DMS
-	Protective_regime: Full access
-	Proximity_human_settlements: Far/Near

N_America_mollusks: Non-modified

RLS para CMelian_02102019: Non-modified

North America: executable file empty

Phorcus spp (sea snail): executable file empty
