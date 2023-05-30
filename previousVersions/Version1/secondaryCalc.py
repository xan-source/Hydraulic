def calculate_open_channel_flow(channel, water_level_outside, basin_water_level):
    """
    This functions calculates discharge past a barrier in case it is connected through an open channel

    NOTES/QUESTIONS:
    - He originally used water_level_inside rather than basin_water_level? I think this should be the same right?

    :param channel: connecting body of water between outside water level and inside water level
    :param water_level_outside: absolute height of water level on the outer (sea) side of the channel [m+MSL}
    :param basin_water_level: absolute height of water level on the inner side of the channel [m+MSL]

    :connections  calculate_storm_surge
    :return q: discharge through the channel [m^3/m/s]
    """

    d_water = channel.Depth + (water_level_outside + basin_water_level) / 2
    diff = water_level_outside - basin_water_level
    q = 0
    if diff != 0 and d_water > 0:
        sign = diff / abs(diff)
        surf_c = channel.Width * d_water
        radius = surf_c / (2 * d_water + channel.Width)  # de Vries
        factor = 0.5
        tmp = abs(diff) / (factor + 10 / (channel.Chezy ** 2) * channel.Length / radius) * 10 * surf_c ** 2
        q = sign * math.sqrt(tmp)
    return q


def formula_broad_weir(land, water_level_outside, basin_water_level):
    """
    Calculates discharge past a barrier (measure/ location) in case of a broad weir.

    NOTES/QUESTIONS:
    - He still mentioned channel here, rather than basin. I think the notes are just bad
    - Why does he


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

# def get_drain_drop_off(Basin, outside_water_level, basins_dict, scenario, j):
#     """
#
#     :param self:
#     :param outside_water_level:
#     :param basins_dict:
#     :param scenario:
#     :param j:
#     :return:
#     """
#     if Basin.DrainsToBasin == 'sea':
#         drain_drop_off = Basin.ScenarioWaterLevels[scenario][j] + Basin.OutletElevation - outside_water_level
#         drain_to_basin = 0
#     else:
#         try:
#             drain_to_basin = int(Basin.DrainsToBasin)
#             drain_drop_off = Basin.ScenarioWaterLevels[scenario][j] + max(Basin.ToBasinRelativeElevation, 0) - \
#                              basins_dict[drain_to_basin].ScenarioWaterLevels[scenario][j]
#
#         except ValueError:
#             print('wrong DrainToBasin input', Basin.DrainsToBasin)
#             return
#     return [drain_to_basin, drain_drop_off]
