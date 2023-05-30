'''
The formulas used to calculate the volume for RUN1.ipynb are located here. Additionally, two formulas to support the calculation of storm surge are defined.



def calculate_inflow_rain(Basin, rain_intensity, timestep):

    def calculate_open_channel_flow(channel, water_level_outside, basin_water_level):
    def formula_broad_weir(land, water_level_outside, basin_water_level):
def calculate_storm_surge(Basin, basin_water_level, LineOfDefense, outside_waterlevel, timestep):

def calculation_drainage(Basin, Volume, timestep, drain_to_water_level_abs, basin_water_level_abs, retention):

def calculate_infiltration(Basin, timestep, abs_basin_wl):

def calculate_interbasin(Basin, basin_water_level, surrounding_water_levels):

'''
import math as math
import numpy as np
from scipy import interpolate

def calculate_inflow_rain(Basin, rain_intensity, timestep):
    """
    INFLOW
      Calculates the volume of rain falling on a basin during 1 timestep.

      :param Basin: Basin in which the rain is falling
      :param rain_intensity: amount of rain falling per hour [mm/hour]
      :param timestep: amount of seconds per timestep

      :return inflow_volume_rain:  volume of rain entering the basin in 1 timestep [m^3]
    """

    rain_volume = Basin.Area * (rain_intensity / 1000) * (timestep)
    return rain_volume    #[m^3/timestep
def calculate_open_channel_flow(channel, water_level_outside, basin_water_level):
    """
    This functions calculates discharge past a barrier in case it is connected through an open channel

    NOTES/QUESTIONS:
    - He originally used water_level_inside rather than basin_water_level? I think this should be the same right?
    - CHECK formula!! It looks like he used Chezy, and just

    :param channel: connecting body of water between outside water level and inside water level
    :param water_level_outside: absolute height of water level on the outer (sea) side of the channel [m+MSL}
    :param basin_water_level: absolute height of water level on the inner side of the channel [m+MSL]

    :connections  calculate_storm_surge
    :return q: discharge through the channel [m^3/m/s]

    """

    d_water = channel.Depth + (water_level_outside + basin_water_level) / 2         #[m]
    diff = water_level_outside - basin_water_level                                  #[m]
    q = 0
    if diff != 0 and d_water > 0:
        sign = diff / abs(diff)
        surf_c = channel.Width * d_water                                            #[m^2]
        radius = surf_c / (2 * d_water + channel.Width)  # de Vries                 #[m]
        factor = 0.5
        tmp = abs(diff) / (factor + 10 / (channel.Chezy ** 2) * channel.Length / radius) * 10 * surf_c ** 2
        q = sign * math.sqrt(tmp)
    return q #[m^3/m/s]
def formula_broad_weir(land, water_level_outside, basin_water_level):
    """
    Calculates discharge past a barrier (measure/ location) in case of a broad weir.

    NOTES/QUESTIONS:
    - He still mentioned channel here, rather than basin. I think the notes are just bad

    more information on : http://cirpwiki.info/wiki/Weirs
    :param land: the barrier that the water has to pass
    :param water_level_outside: absolute height of water level on the outer (sea) side of the channel [m+MSL}
    :param basin_water_level: absolute height of water level on the inner side of the channel [m+MSL]

    :connections calculate_storm_surge
    :return q: discharge past the barrier [m^3/m/s]
    """

    q = 0
    h_upper = max(water_level_outside, basin_water_level)
    h_lower = min(water_level_outside, basin_water_level)

    if land.Type == 'Water':
        height = 0
    else:
        height = land.Height

    if h_upper > height:
        C_w = 0.55                                                                      # HEC2010 predicts 0.46-0.55 for broad-crested weirs
        signum = 1
        if h_lower / h_upper <= 0.67:
            C_df = 1
        else:
            C_df = 1 - 27.8 * (h_lower / h_upper - 0.67) ** 3

        if basin_water_level > water_level_outside:
            signum = -1                                                                 # water flows outwards
        q = signum * C_df * C_w * math.sqrt(9.81) * (h_upper - height) ** (3 / 2)

    return q
def calculate_storm_surge(Basin, basin_water_level, LineOfDefense, outside_waterlevel, timestep):
    """
    This function simulates storm surge hitting a Line of Defense. For 1 timestep, it calculates the volume of water
    passing the Line of Defense. The volume is the total amount of volume flowing into 1 basin. On the Line of Defense,
    a flood defense can be placed. The volume is calculated both for the situation where it holds and fails.
    :param: Basin: The drainage basin that receives the volume of water
    :param: basin_water_level: Absolute water level in the basin at the start of the timestep [m+MSL]
    :param: location: Part of Line of Defense where the storm surge hits. Can be Land (e.g. Coast) or Water (e.g. River)
    :param: outside_waterlevel: Water level on the outer side of the Line of Defense at the time of timestep. [m+MSL]
    :param: timestep: length of a time step [s]

    Important variables:
    Q_hold: Discharge of water entering the basin per second in case a measure holds [m^3/s]

    Important functions:
    formula_broad_weir: The present land (without defense) is modeled as a broad weir.
                        This formula calculates the discharge in this case [m^3/m/s]
    get_overtopping: Could be included. Calculates overtopping, but not very relevant in this case
    calculate_open_channel_flow: Calculates the discharge through the channel.
                                This can be both negative (out of the basin) or positive (into the basin) [m^3/m/s]
    :return V_hold: Volume of water entering the basin during the timestep in case a measure holds [m^3]
    """
    # Q_hold = None
    Q_open = 0
    if LineOfDefense.Type == 'Land':
        if basin_water_level > LineOfDefense.Height:
            Q_open += formula_broad_weir(LineOfDefense, outside_waterlevel, basin_water_level) * Basin.Width  # [m3/s]

    elif LineOfDefense.Type == 'Water':
        Q_open += calculate_open_channel_flow(LineOfDefense, outside_waterlevel, basin_water_level)

    # V_hold = Q_hold * timestep
    V_open = Q_open * 3600 * timestep               #change from per second to per timestep
    return V_open
def calculation_drainage(Basin, Volume, timestep, drain_to_water_level_abs, basin_water_level_abs, retention):
    """
    This function calculates the volume that is drained from a drainage basin, as well as the remaining volume. The
    water can drain to another basin (flow_interbasin) or out to sea/sewerage (outflow_drainage). The water level in the
    drainage basin is compared with the water level of the receiving body (basin/sea). If the receiving water level is
    higher, no water flows. In other cases, a volume of water drains from the basin to the receiving basin. This volume
    is less when the difference in water level is small (<1 m).

    Uses formula 3,4 from World Bank paper

    NOTES/QUESTIONS:
    - Why does this include both basin_water_level, and volume if there is an easy way to convert this through Basin.VolumeToHeight?

    :param Basin: The current drainage basin
    :param drain_to_basin: the number of the drainage basin that receives the water. '0' represents sea/sewerage/none
    :param Volume: The total volume of water inside the basin at the start of the calculation [m^3]
    :param timestep: length of time of a model timestep [s]
    :param drain_to_water_level: Water level in the drainage basin that receives the drained volume of water [m+MSL]
    :param basin_water_level: Water level in the drainage basin that drains [m+MSL]
    :param retention: Total volume of retention inside the basin at start of calculation [m^3]

    :return outflow_drainage: volume of water that drains out of the system [m^3]
    :return remaining_volume: volume of water that remains in the drainage basin [m^3]
    :return flow_interbasin: volume of water that flow to another drainage basin [m^3]
    :return retention: New volume of retention after calculations. [m^3]
    """
    hdiff = (basin_water_level_abs) - drain_to_water_level_abs                   #can go negative
    drain_factor = max(min(hdiff, 1), 0)                                                #Capped between 0 and 1

    max_drainage_capacity = Basin.DrainageDischarge * 3600 * timestep * drain_factor #m^3/timestp
    outflow_drainage = min(Volume + retention, max_drainage_capacity)           #if drainage > volume in basin, only the volume + retention in the basin drain.
    remaining_volume = Volume - outflow_drainage

    if remaining_volume < 0:                                                    # negative remaining volume, so retention decreases
        retention += remaining_volume
        remaining_volume = 0

    return [outflow_drainage, retention]
def calculate_infiltration(Basin, timestep, abs_basin_wl):
    """
    OUTFLOW
    Calculates the infiltration of each basin for each timestep.
    If there is no water, no infiltration. If basin is full, maximum infiltration.

    :param Basin: Basin where infiltration is taking place
    :param timestep: amount of seconds per timestep [hours]
    :param basin_water_level: the water level in the basin at the start of the calculation [m]
    :param infiltration_rate: the amount of water that infiltrates in the basin[mm/hour]

    :return: outflow_volume_infiltration: volume of water that infiltrates basin in 1 time step [m^3]
    """
# This was used to calculate infiltration in FLORES.
# I think a better method could be possible though as this was for an urban landscape
    min_x = Basin.Contours[0].MinHeight
    max_x = Basin.Contours[-1].MinHeight
    x = min(Basin.Contours[-1].MinHeight, abs_basin_wl)

    approx_surface_area = ((x - min_x) / (max_x - min_x)) * Basin.Area
    outflow_volume_infiltration = Basin.InfiltrationRate / 1000 * approx_surface_area * (timestep)
    return outflow_volume_infiltration      #[m^3/timestep]
def calculate_interbasin(Basin, basin_water_level, surrounding_water_levels):
    """
    In original FLORES model this was an aspect of the Basin class. I changed it to be a formula.
    Based on the water levels of the current and surrounding basins, and the border heights around the basin
    :param Basin: Current basin from which the flow is calculated
    :param basin_water_level: The water level at the current timestep in the basin. [m+MSL]
    :param surrounding_water_levels: Water levels of neighboring basins. [m+MSL]

    :return: new_absolute_water_level: Water level in basin after calculation [m+MSL]
    :return: flow_to_other_basins: Volume of water going into other basins [m^3]
    """
    flow_to_other_basins = {}
    #Calculates surrounding water levels and border heights, and uses the maximum of these two to create a threshold
    surrounding_thresholds = {b: max(surrounding_water_levels[b], Basin.BorderHeights[b]) for b in
                              surrounding_water_levels.keys()}

    #Not sure what the water levels are sorted by. After they're sorted a number of lists is made in the same order      Background info: https://www.programiz.com/python-programming/methods/built-in/sorted
    sorted_thresholds = sorted(surrounding_thresholds.items(), key=lambda kv: kv[1])
    basin_numbers_sorted = []
    surrounding_basins_threshold = []
    surface_areas_sorted = []
    for i in range(len(sorted_thresholds)):
        basin_numbers_sorted.append(sorted_thresholds[i][0])
        surrounding_basins_threshold.append(sorted_thresholds[i][1])
        surface_areas_sorted.append(Basin.SurroundingAreas[sorted_thresholds[i][0]])

    # Using the absolute water level from the last time step and the sorted thresholds, for each consecutive threshold the
    # flow between those two basins is calculated and saved, and the new water level in the source basin is changed for the next threshold.
    # Based on equation 5 from worldbank document.
    new_absolute_water_level = basin_water_level
    for j in range(len(surrounding_basins_threshold)):
        tmp_absolute_water_level = new_absolute_water_level
        if tmp_absolute_water_level > surrounding_basins_threshold[j]:
            area_factor = Basin.Area / surface_areas_sorted[j]
            flow_to_other_basins[basin_numbers_sorted[j]] = Basin.HeightToVolume(tmp_absolute_water_level) - Basin.HeightToVolume(new_absolute_water_level)
            new_absolute_water_level = (tmp_absolute_water_level * area_factor + surrounding_basins_threshold[j]) / (1 + area_factor)

        else:
            flow_to_other_basins[basin_numbers_sorted[j]] = 0

             #This can be done without missing a threshold because they've already been sorted

    return flow_to_other_basins



