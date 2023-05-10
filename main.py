#See PyCharm help at https://www.jetbrains.com/help/pycharm/
import numpy as np

from calculations import create_surge_series
from calculations import calculate_inflow_rain, calculate_storm_surge
from calculations import calculate_infiltration, calculation_drainage





#inflow_volume_rain
storm_duration          =                               #[hours]   Heb ik
rain_intensity          =                               #[mm/hour] Heb ik


#inflow_storm_surge
max_surge               =                               #[m]       Heb ik
tidal_amplitude         =                               #[m]       Heb ik
mean_sea_level          =                               #[m]       Heb ik
outside_waterlevel      = create_surge_series(max_surge, storm_duration, tidal_amplitude, timestep, mean_sea_level)   #[m]

#General
timestep                =                               #[hours]   Heb ik
time                    = storm_duration/timestep       #Easy calculation to get a calculation for each timestep
waterlevel_before       =                               #[m]


LineOfDefense           =                               #[-] List of all Lines of Defense (even though there aren't any defenses, the storm surge is calculated as a broad weir or open channel flow)
listbasin               =                               #[-] List of all basins.
basin_water_level       =[]                             #empty list of water levels, to be filled in from within the loop
retention               =[]                             #empty list of retention at each step.
Volume                  =[]
for i in range(len(time)):
        for basin in listbasin:
#Set variables
        drain_to_water_level            = Basin.DrainsTo[basin_water_level] ## IDK if this is correct, but idea is to get water level for basin where this basin drains to.
        surrounding_water_levels        = Basin.                            ## I want to call the water levels in neighboring basins in a similar fashion as above

#Different inflows: rain, storm surge as hydraulic boundaries, surface flow and drainage from other basins
        inflow_volume_rain = calculate_inflow_rain(Basin, rain_intensity, timestep)
        inflow_storm_surge = calculate_storm_surge(Basin, basin_water_level, LineOfDefense, outside_waterlevel, timestep)

        ####HOW TO CALCULATE THESE TWO?? I have the inflow as a result of calculation_drainage and calculate_interbasin.
        # Iets met i ofzo? inflow_drainage = outflow_drainage[Basin.DrainsTo][i]??
        inflow_neighbors   =
        inflow_drainage    =


#Different outflows: retention, drainage, infiltration
        outflow_retention  = min(Volume, Basin.RetentionCapacity - retention[i])                                #Heel omslachtig in code Erik, nog een keer checken, maar volgens mij klopt het nu
        [outflow_drainage, retention[i+1]]  = calculation_drainage(Basin, Volume, timestep, drain_to_water_level, basin_water_level, retention)
        outflow_infiltrate = calculate_infiltration(Basin, timestep, basin_water_level)
        outflow_neighbors  = calculate_interbasin(Basin, basin_water_level, surrounding_water_levels)

#Hydrological balance (van Berchum et al. 2020)
        Volume[i+1] = Volume[i] + (inflow_volume_rain + inflow_storm_surge + inflow_neighbors + inflow_drainage
                                - outflow_infiltrate - outflow_drainage - outflow_neighbors- outflow_retention)
        basin_water_level = Basin.VolumeToHeight(Volume)

