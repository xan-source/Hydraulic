"""
The formulas used to calculate the volume for RUN1.ipynb are located here. Additionally, two formulas to support the calculation of storm surge are defined.



def calculate_inflow_rain(Basin, rain_intensity, timestep):

    def calculate_open_channel_flow(channel, water_level_outside, basin_water_level):
    def formula_broad_weir(land, water_level_outside, basin_water_level):
def calculate_storm_surge(Basin, basin_water_level, LineOfDefense, outside_waterlevel, timestep):

def calculation_drainage(Basin, Volume, timestep, drain_to_water_level_abs, basin_water_level_abs, retention):

def calculate_infiltration(Basin, timestep, abs_basin_wl):

def calculate_interbasin(Basin, basin_water_level, surrounding_water_levels):

"""
import math as math
import numpy as np
from scipy import interpolate


def calculate_inflow_rain(Basin, rain_intensity, TIMESTEP):
    """Calculates the inflow from rain for one timestep

    Parameters
    ----------
    Basin : class
        the basin who's rainfall is being calculated
    rain_intensity : float
        the rainfall for the storm's duration [mm/hour]
    TIMESTEP : float
        the timestep the model takes [hours]

    Returns
    -------
    inflow_volume_rain: float
        volume of rain entering the basin in 1 timestep[m3/timestep]
    """

    rain_volume = Basin.area * (rain_intensity / 1000) * (TIMESTEP)
    return rain_volume


"""
    This functions calculates discharge past a barrier in case it is connected through an open channel

    NOTES/QUESTIONS:
    - He originally used water_level_inside rather than basin_water_level? I think this should be the same right?
    - CHECK formula!! It looks like he used chezy, and just

    :param LineOfDefense_water: water based line of defense object
    :param water_level_outside: absolute height of water level on the outer (sea) side of the channel [m+MSL}
    :param absolute_basin_water_level: absolute height of water level on the inner side of the channel [m+MSL]

    :connections  calculate_storm_surge
    :return q: discharge through the channel [m^3/m/s]

    """


def calculate_open_channel_flow(
    LineOfDefense_water, water_level_outside, absolute_basin_water_level
):
    """This functions calculates discharge past a barrier in case it is connected through an open channel.

    This formula is a secondary calculation to the calculate_storm_surge function.

    Parameters
    ----------
    LineOfDefense_water : class
        The line of defense. Care should be taken that this LOD should be water-based
    water_level_outside : float
        The water level outside of the basin, on the unprotected side of the LOD [m]
    absolute_basin_water_level : float
        Absolute water level inside the basin [m]

    Returns
    -------
    q : float
        flow coming in or going out of the basin. [m3/m/s]

    See Also
    --------
    formula_broad_weir, calculate_storm_surge

    Notes
    -----
    I think this is a modified manning equation, but it's not clear where the formula comes from. Might have to change it

    References
    ----------
    Rapid screening and evaluation of flood risk reduction strategies: Exploratory study on the use of the FLORES modelling approach for World Bank projects - Erik van Berchum
    """

    d_water = (
        LineOfDefense_water.depth
        + (water_level_outside + absolute_basin_water_level) / 2
    )  # [m]
    diff = water_level_outside - absolute_basin_water_level  # [m]
    q = 0
    if diff != 0 and d_water > 0:
        sign = diff / abs(diff)
        surf_c = LineOfDefense_water.width * d_water  # [m^2]
        radius = surf_c / (
            2 * d_water + LineOfDefense_water.width
        )  # de Vries                 [m]
        factor = 0.5
        tmp = (
            abs(diff)
            / (
                factor
                + 10
                / (LineOfDefense_water.chezy**2)
                * LineOfDefense_water.length
                / radius
            )
            * 10
            * surf_c**2
        )
        q = sign * math.sqrt(tmp)


    return q  # [m^3/m/s]


def formula_broad_weir(
    LineOfDefense_land, water_level_outside, absolute_basin_water_level
):
    """Calculates the flow over a barrier when it acts as a broad weir.

    This function is a secondary calculation to the calculate_storm_surge function.

    Parameters
    ----------
    LineOfDefense_land : class
        The line of defense. Care should be taken that this LOD should be water-based
    water_level_outside : float
        absolute water level on the other side of the weir [m]
    absolute_basin_water_level : float
        the absolute water level in the basin [m]

    Returns
    -------
    q : float
        water that either flows into or out of the basin at each timestep [m3/m/s]

    See Also
    --------
    calculate_open_channel_flow, calculate_storm_surge

    Notes
    -----
    Have to double check the formula at the reference site.

    References
    ----------
    more information on : http://cirpwiki.info/wiki/Weirs
    Rapid screening and evaluation of flood risk reduction strategies: Exploratory study on the use of the FLORES modelling approach for World Bank projects - Erik van Berchum

    """

    q = 0
    h_upper = max(water_level_outside, absolute_basin_water_level)
    h_lower = min(water_level_outside, absolute_basin_water_level)

    if LineOfDefense_land.type == "Water":
        height = 0
    else:
        height = LineOfDefense_land.height

    if h_upper > height + 1:
        C_w = 0.55  # HEC2010 predicts 0.46-0.55 for broad-crested weirs
        signum = 1
        if h_lower / h_upper <= 0.67:
            C_df = 1
        else:
            C_df = 1 - 27.8 * (h_lower / h_upper - 0.67) ** 3

        if absolute_basin_water_level > water_level_outside:
            signum = -1  # water flows outwards
        q = signum * C_df * C_w * math.sqrt(9.81) * (h_upper - height) ** (3 / 2)

    return q


def calculate_storm_surge(
    Basin, absolute_basin_water_level, LineOfDefense, outside_waterlevel, TIMESTEP
):
    """This function simulates storm surge hitting a Line of Defense.

    For 1 timestep, it calculates the volume of water passing the Line of Defense. The volume is the total amount of
    volume flowing into 1 basin. On the Line of Defense, a flood defense can be placed. The volume is calculated both
    for the situation where it holds and fails.

    Parameters
    ----------
    Basin : class
        The drainage basin the storm surge is calculated for
    absolute_basin_water_level : float
        absolute water level in the basin [m]
    LineOfDefense : class
        the line of defense that protects the basin from the outside water level
    outside_waterlevel : float
        absolute waterlevel outside the basin [m]
    TIMESTEP : float
        a constant. This represents the timestep taken in the model [hours]

    Returns
    -------
    V_open : float
        volume of water entering basin [m3]

    """
    # Q_hold = None
    Q_open = 0
    if LineOfDefense.type == "Land":
        Q_open += (
            formula_broad_weir(
                LineOfDefense, outside_waterlevel, absolute_basin_water_level
            )
            * LineOfDefense.width
        )  # [m3/s]

    if LineOfDefense.type == "Water":
        Q_open += (
            calculate_open_channel_flow(
                LineOfDefense, outside_waterlevel, absolute_basin_water_level
            )
        )

    # V_hold = Q_hold * timestep
    V_open = Q_open * 3600 * TIMESTEP
    return V_open


def calculation_drainage(
    Basin,
    volume,
    TIMESTEP,
    drain_to_water_level_abs,
    absolute_basin_water_level,
    retention,
):
    """This function calculates the volume that is drained from a drainage basin, as well as the remaining volume.

    If the receiving water level is higher, no water flows. In other cases, a volume of water drains from the
    basin to the receiving basin. This volume is less when the difference in water level is small (<1 m).

    Parameters
    ----------
    Basin : class
        current basin
    volume : float
        current volume in the basin [m3]
    TIMESTEP : float
        a constant. This represents the timestep taken in the model [hours]
    drain_to_water_level_abs : float
        absolute water level in the basin where the current basin drains to [m]
    absolute_basin_water_level : float
        absolute water level in the basin [m]
    retention : float
        retention in the basin before drainage [m3]

    Returns
    -------
    outflow drainage : float
        volume of water leaving the basin through drainage [m3]
    retention : float
        the new retention, after drainage [m3]
    """
    hdiff = (absolute_basin_water_level) - drain_to_water_level_abs  # can go negative
    drain_factor = max(min(hdiff, 1), 0)  # Capped between 0 and 1


    max_drainage_capacity = (
        Basin.drainage_discharge * 3600 * TIMESTEP * drain_factor
    )  # m^3/timestp
    outflow_drainage = min(
        volume + retention, max_drainage_capacity
    )  # if drainage > volume in basin, only the volume + retention in the basin drain.
    remaining_volume = volume - outflow_drainage

    if remaining_volume < 0:  # negative remaining volume, so retention decreases
        retention += remaining_volume
        remaining_volume = 0

    return [outflow_drainage, retention]


def calculate_infiltration(Basin, TIMESTEP, absolute_basin_water_level):
    """Calculates the infiltration of each basin for each timestep.
    If there is no water, no infiltration. If basin is full, maximum infiltration.

    Parameters
    ----------
    Basin : class
        the basin at hand
    TIMESTEP : float
        a constant. This represents the length of a timestep in the model [hours]
    absolute_basin_water_level : float
        absolute water level in the basin [m]

    Returns
    -------
    outflow_volume_infiltration : float
        the volume of water infiltrating at a timestep [m3]
    """
    min_x = Basin.contours[0].min_height
    max_x = Basin.contours[-1].min_height
    x = min(Basin.contours[-1].min_height, absolute_basin_water_level)

    approx_surface_area = ((x - min_x) / (max_x - min_x)) * Basin.area
    outflow_volume_infiltration = (
        Basin.infiltration_rate / 1000 * approx_surface_area * (TIMESTEP)
    )
    return outflow_volume_infiltration  # [m^3/timestep]


def calculate_interbasin(
    Basin, absolute_basin_water_level, absolute_surrounding_water_levels
):
    """Calculates flow to other basins based on their respective water level

    First the function calculates surrounding water levels and border heights, and uses the maximum of these two to
    create a threshold
    Then the water levels are sorted, after which a number of lists is made in the same order
    Using the absolute water level from the last time step and the sorted thresholds, for each consecutive threshold the
    flow between those two basins is calculated and saved, and the new water level in the source basin is changed for
    the next threshold.

    Parameters
    ----------
    Basin : class
        basin for which calculation is done.
    absolute_basin_water_level : float
        absolute water level in basin [m]
    absolute_surrounding_water_levels : dict
        dictionary with the absolute water level for the neighboring basins [m]

    Returns
    -------
    flow_to_other_basins : dict
        a dict with the neighboring basins and how much water flows to them for the current timestep [m3]

    References
    ----------
    Rapid screening and evaluation of flood risk reduction strategies: Exploratory study on the use of the FLORES
    modelling approach for World Bank projects - Erik van Berchum

    """
    flow_to_other_basins = {}
    surrounding_thresholds = {
        b: max(absolute_surrounding_water_levels[b], Basin.border_heights[b])
        for b in absolute_surrounding_water_levels.keys()
    }

    sorted_thresholds = sorted(surrounding_thresholds.items(), key=lambda kv: kv[1])
    basin_numbers_sorted = []
    surrounding_basins_threshold = []
    surface_areas_sorted = []
    for i in range(len(sorted_thresholds)):
        basin_numbers_sorted.append(sorted_thresholds[i][0])
        surrounding_basins_threshold.append(sorted_thresholds[i][1])
        surface_areas_sorted.append(Basin.surrounding_areas[sorted_thresholds[i][0]])

    # Based on equation 5 from worldbank document.
    new_absolute_water_level = absolute_basin_water_level
    for j in range(len(surrounding_basins_threshold)):
        tmp_absolute_water_level = new_absolute_water_level

        if tmp_absolute_water_level > surrounding_basins_threshold[j]:
            area_factor = Basin.area / surface_areas_sorted[j]

            #to avoid negative water levels, the water level can't go below lower contour
            new_absolute_water_level = max((
                tmp_absolute_water_level * area_factor + surrounding_basins_threshold[j]
            ) / (1 + area_factor), Basin.contours[0].min_height)

            flow_to_other_basins[basin_numbers_sorted[j]] = Basin.height_to_volume(
                tmp_absolute_water_level
            ) - Basin.height_to_volume(new_absolute_water_level)



        else:
            flow_to_other_basins[basin_numbers_sorted[j]] = 0

            # This can be done without missing a threshold because they've already been sorted

    return flow_to_other_basins
