from secondaryCalc import formula_broad_weir, calculate_open_channel_flow
from scipy import interpolate

'''
Classes are not filled in yet, except where mentioned "done".

Everything that is called in the rest of the formulas is mentioned, however. 
'''
class Basin:
    def __init__(self, name, surfacearea=0, contours=[], width=0):
        self.Area               = surfacearea               # Area of the basin.
        self.Width              = width                     # used for storm surge calculation
        self.BorderHeights      ={}                         # the heights of the borders with surrounding basins.
        self.InfiltrationRate   = None                      #done the infiltration rate in the basin. Currently set to 2.5, but can be calculated with a different technique as well
        self.Contours           = contours                  # Don't know how Erik linked this to Contours
        self.DrainageCapacity   = None
        self.DrainsTo           = None
        self.VolumeToHeight     = None                      #done
        self.HeightToVolume     = None                      #done
        self.BorderHeights      = None
        self.SurroundingAreas   = None
        self.RetentionCapacity  = 0
    def get_infiltration_rate(self):

        self.InfiltrationRate = 2.5  # use until we have better method
        return

    def get_volume_inundation_curve(self):
        heights = [self.Contours[0].MinHeight]
        surface_areas = [self.Contours[0].Area]
        volumes = [0]
        for contour in range(1, len(self.Contours)):
            heights.append(self.Contours[contour].MinHeight)
            surface_areas.append(round(surface_areas[contour - 1] + self.Contours[contour].SurfaceArea, 1))
            volumes.append(volumes[contour - 1] + (heights[contour] - heights[contour - 1]) * surface_areas[contour])
        heights.append(30)
        surface_areas.append(surface_areas[-1])
        volumes.append(volumes[-1] + (heights[-1] - heights[-2]) * surface_areas[-1])
        curve_volume_to_height = interpolate.interp1d(volumes, heights, kind='linear')
        curve_height_to_volume = interpolate.interp1d(heights, volumes, kind='linear')
        self.VolumeToHeight = curve_volume_to_height
        self.HeightToVolume = curve_height_to_volume
        return



class Contours:
    def __init__(self): # nog beter invullen.
        self.Code                                           #Not called anywhere
        self.MinHeight
        self.SurfaceArea


class LineOfDefense:
    def __init__(self):
        self.Type               # used to mention whether land or water based
        self.Height
        #Only valid in cases where water based (schematized as a channel)
        self.Depth
        self.Width
        self.Chezy
        self.Length             #Length of ...


def create_surge_series(max_surge, storm_duration, tidal_amplitude, timestep, mean_sea_level):
    """
    This function makes a time series of the surge, consisting of tide and a surge. The surge is schematized as a
    half sine. the tide is schematized to have its peak halfway through the storm duration (when surge is max)

    :param max_surge: the maximum additional surge of the water level because of the storm [m]
    :param storm_duration: length of the storm [hours]
    :param tidal_amplitude: normal amplitude of the tide [m]
    :param timestep: length of one timestep [s]
    :param mean_sea_level: absolute height of the mean sea level [m]

    :return: time series of water levels [m+MSL]
    """

    TIDE_INTERVAL = 12.417  # 12 hours and 25 minutes
    time1 = np.linspace(0, storm_duration, int(storm_duration / timestep + 1), endpoint=True)
    tide = tidal_amplitude * np.cos(np.pi / (TIDE_INTERVAL / 2) * time1 - 0.5 * storm_duration) + mean_sea_level
    surge = max_surge * np.sin(np.pi / storm_duration * time1) ** 5
    total_storm_surge_series = tide + surge

    return list(total_storm_surge_series)

def calculate_inflow_rain(Basin, rain_intensity, timestep):
    """
    INFLOW
      Calculates the volume of rain falling on a basin during 1 timestep.

      :param Basin: Basin in which the rain is falling
      :param rain_intensity: amount of rain falling per hour [mm/hour]
      :param timestep: amount of seconds per timestep

      :return inflow_volume_rain:  volume of rain entering the basin in 1 timestep [m^3]
    """

    rain_volume = Basin.Area * (rain_intensity / 1000) * (timestep/3600)
    return rain_volume

def calculate_storm_surge(Basin, basin_water_level, LineOfDefense, outside_waterlevel, timestep):
    """
    This function simulates storm surge hitting a Line of Defense. For 1 timestep, it calculates the volume of water
    passing the Line of Defense. The volume is the total amount of volume flowing into 1 basin. On the Line of Defense,
    a flood defense can be placed. The volume is calculated both for the situation where it holds and fails.

    :param Basin: The drainage basin that receives the volume of water
    :param basin_water_level: Absolute water level in the basin at the start of the timestep [m+MSL]
    :param location: Part of Line of Defense where the storm surge hits. Can be Land (e.g. Coast) or Water (e.g. River)
    :param outside_waterlevel: Water level on the outer side of the Line of Defense at the time of timestep. [m+MSL]
    :param timestep: length of a time step [s]

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
    Q_open = None

    if LineOfDefense.Type == 'Land':
        if basin_water_level > LineOfDefense.Height + 1:
            Q_open = formula_broad_weir(LineOfDefense, outside_waterlevel, basin_water_level) * Basin.Width  # [m3/s]
        # else:
        #     Q_open = location.get_overtopping(outside_waterlevel, location.Height, waveheight) * Basin.Width

        # if measure:
        #     if waterlevel_before > measure.Height + 1:  # bay level is higher than barrier, so it acts as broad weir
        #         Q_hold = formula_broad_weir(measure, outside_waterlevel, waterlevel_before) * Basin.Width
        #     else:  # bay level is lower than the barrier, so it overtops/overflows
        #         Q_hold = measure.get_overtopping(outside_waterlevel, measure.Height, waveheight) * Basin.Width
        #
        #     if time < time_fail_land:
        #         Q_open = Q_hold
        #     else:
        #         Q_open = Q_open * measure.BreachFactor + Q_hold * (1 - measure.BreachFactor)
        # else:
        #     Q_hold = Q_open

    elif LineOfDefense.Type == 'Water':
        Q_open = calculate_open_channel_flow(LineOfDefense, outside_waterlevel, basin_water_level)

        # if measure and outside_waterlevel > h_close:
        #     if waterlevel_before > measure.Height + 1:  # Acts like a weir in this situation
        #         Q_hold = formula_broad_weir(location, outside_waterlevel,
        #                                     waterlevel_before) * Basin.Width
        #     else:  # bay level is lower than surge and lower than barrier, so it overtops/overflows
        #         Q_hold = measure.get_overtopping(outside_waterlevel, measure.Height,
        #                                          waveheight) * Basin.Width
        #     if time < time_fail_water:
        #         Q_open = Q_hold
        # else:
        #     Q_hold = Q_open

    # V_hold = Q_hold * timestep
    V_open = Q_open * timestep

    return [V_open]



def calculation_drainage(Basin, Volume, timestep, drain_to_water_level, basin_water_level, retention):
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
    hdiff = basin_water_level - drain_to_water_level
    drain_factor = min(hdiff, 1)

    # drain_factor = 1                # if drain_drop_off is more than 1, nothing happens in following statement and factor = 1
    # if drain_drop_off < 0:          # if water level downstream is higher than in basin, no outflow drainage/surface
    #     return [0, volume, 0, retention]
    # elif drain_drop_off < 1.0:      # if
    #     drain_factor = min(0.25 + drain_drop_off*3/4, 1) #not sure where this formula comes from

    max_drainage_capacity = Basin.DischargeCapacity * timestep * drain_factor
    outflow_drainage = min(Volume + retention, max_drainage_capacity)           #if drainage > volume in basin, only the volume + retention in the basin drain.
    remaining_volume = Volume - outflow_drainage

    if remaining_volume < 0:                                                    # negative remaining volume, so retention decreases
        retention += remaining_volume
        remaining_volume = 0

    inflow_drainage = 0
    if Basin.DrainsTo != 0:                                                     # if water flows to another basin ( 0 represents 'sea' or 'none')
        inflow_drainage = outflow_drainage
        outflow_drainage = 0

    return [outflow_drainage, retention]
def calculate_infiltration(Basin, timestep, basin_water_level):
    """
    OUTFLOW
    Calculates the infiltration of each basin for each timestep.
    If there is no water, no infiltration. If basin is full, maximum infiltration.

    :param Basin: Basin where infiltration is taking place
    :param timestep: amount of seconds per timestep [hours]
    :param basin_water_level: the water level in the basin at the start of the calculation [m]
    :param infiltration_rate: the amount of water that infiltrates in the basin[m/hour]

    :return: outflow_volume_infiltration: volume of water that infiltrates basin in 1 time step [m^3]
    """
# This was used to calculate infiltration in FLORES.
# I think a better method could be possible though as this was for an urban landscape
    min_x = Basin.Contours[0].MinHeight
    max_x = Basin.Contours[-1].MinHeight
    x = min(Basin.Contours[-1].MinHeight, basin_water_level + Basin.Contours[0].MinHeight)

    approx_surface_area = ((x - min_x) / (
                max_x - min_x)) * Basin.Area
    outflow_volume_infiltration = Basin.InfiltrationRate / 1000 * approx_surface_area * (timestep / 3600)
    return outflow_volume_infiltration


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
    surrounding_thresholds = {Basin: max(surrounding_water_levels[Basin], Basin.BorderHeights[Basin]) for Basin in
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
    surrounding_basins_threshold.append(100)

    # Using the absolute water level from the last time step and the sorted thresholds, for each consecutive threshold the
    # flow between those two basins is calculated and saved, and the new water level in the source basin is changed for the next threshold.
    # Based on equation 5 from worldbank document.
    new_absolute_water_level = basin_water_level
    for j in range(len(surrounding_basins_threshold)):
        tmp_absolute_water_level = new_absolute_water_level
        if tmp_absolute_water_level > surrounding_basins_threshold[j]:
            area_factor = Basin.Area / surface_areas_sorted[j]
            new_absolute_water_level = (tmp_absolute_water_level * area_factor + surrounding_basins_threshold[j]) / (1 + area_factor)
            flow_to_other_basins[basin_numbers_sorted[j]] = Basin.HeightToVolume(tmp_absolute_water_level) - Basin.HeightToVolume(new_absolute_water_level)
        else:
            break #This can be done without missing a threshold because they've already been sorted

    return [flow_to_other_basins]



