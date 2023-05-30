'''
In this file the schematization files are read in, after which the classes are defined.
Additionally, general formulas to mold data and get certain values are defined and explained.

General formulas:
def tolist(input_string, splitby, export_type=int)
def pad_or_truncate(some_list, target_len):
def create_surge_series(max_surge, storm_duration, tidal_amplitude, timestep, mean_sea_level):
def get_absolute_surrounding_water_levels(Basin, i, abs_basin_wl):

Load in basins:
def load_basins(file):
def load_basin_borders(basin_master, file):
def load_drainage_capacities(basin_master, file, small_channel, mid_channel, large_channel, low_drain, mid_drain, high_drain):
class Contour
class Basin
    get_infiltration_rate(self):
    get_volume_inundation_curve(self):
    get_volume_inundation_curve(self):

Load in Coast/Line of defense:
def load_layers(file, basins_master):
class LineOfDefense(object):
class ProtectedArea(object):
    create_list_surface_area(self, dict_basin):
class UnprotectedArea(ProtectedArea):

Load in Hydraulic boundary conditions:
def load_hydraulic_conditions(file_surge, file_rain, timestep):
class HydraulicBoundaryConditions (object):




'''
import pandas as pd
import numpy as np
from scipy import interpolate
def tolist(input_string, splitby, export_type=int):
    """
    This function is used to divide list that has been imported as a string into an actual list. For example, if the
    import is '1;2;3;4', the output is [1,2,3,4]

    :param input_string: string value, consisting of divided numbers
    :param splitby: punctuation mark used to divide the numbers
    'param expert_type: type of the variables in the final list. standard = integer

    :return: list of values
    """

    if type(input_string) == str:
        return list(map(export_type, input_string.split(splitby)))
    else:
        str_variable = str(int(input_string))
        return list(map(export_type, str_variable.split(splitby)))
def pad_or_truncate(some_list, target_len):
    """
    This function shortens or extends a list. When it is extended, extra 0 values are added. This can be helpful to
    simulate what happens after a storm has passed.

    :param some_list: The original list of values
    :param target_len: the wanted length of the list

    :return: a list of values with length of target_len
    """
    return some_list[:target_len] + [0]*(target_len - len(some_list))
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
#Basin
# def load_basins(file):
#     """
#     This function loads in the different basin information.
#     First the CSV is loaded in, and an empty dataframe is built for the basins information.
#     Second the dataframe is filled with information from the file for each
#     #The load_basins function reads a CSV file containing information about basins and their contours. The function takes
#     one argument, file, which is the path to the CSV file. The function uses the pandas library to read the CSV file
#     and store the data in a pandas dataframe called basins_source.
#
#     The function then initializes an empty dictionary called basins_master, which will be used to store the Basin
#     objects. The function iterates through each row in the basins_source dataframe and creates a tmp_basin object
#     containing the information for the current row.
#
#     The function then checks if the Basin_ID of the tmp_basin object is equal to the length of the basins_master
#     dictionary minus one, and if the length of the basins_master dictionary is not zero. If this condition is true,
#     the function updates the existing Basin object in the basins_master dictionary with the new contour information.
#     If the condition is false, the function creates a new Basin object and adds it to the basins_master dictionary.
#
#     After all the Basin objects have been created, the function iterates through each Basin object in the
#     basins_master dictionary and calculates the infiltration rate and volume inundation curve for each Basin object.
#
#     Finally, the function returns the basins_master dictionary containing all the Basin objects.
#
#     The file needs the following information (case sensitive)
#     :param: 'Contour_ID'    ID of contour in the csv
#     :param: 'Height_min'    min height of said contour
#     :param: 'area (m2)'     The area of the contour at hand
#     :param: 'Basin_ID'      which basin the contour is located in
#
#     Dependencies:
#     Basin.Area
#     Basin.Contours
#     :func: get_infiltration_rate()
#     :func: get_volume_inundation_curve()
#     """
#     basins_source = pd.read_csv(file, sep=',', header=0)
#     basins_master = {}
#
#     for row in range(0, len(basins_source)):
#         tmp_basin = basins_source.loc[row]
#
#         if int(tmp_basin['Basin_ID']) == len(basins_master) - 1 and len(
#                 basins_master) != 0:  ##Probably not necessary to do the (basin-1)=basin ID thing
#
#             basins_master[int(tmp_basin['Basin_ID'])].Area += tmp_basin['area (m2)']
#             basins_master[int(tmp_basin['Basin_ID'])].Contours.append(
#                 Contour((tmp_basin['Contour_ID']), tmp_basin['Height_min'], tmp_basin['area (m2)']))
#
#         else:
#         basins_master[int(tmp_basin['Basin_ID'])] = Basin(str(int(tmp_basin['Basin_ID'])),
#                                                               surfacearea=tmp_basin['area (m2)'], contours=[
#                     Contour((tmp_basin['Contour_ID']), tmp_basin['Height_min'], tmp_basin['area (m2)'])])
#
#     for basin in basins_master:
#         basins_master[basin].get_infiltration_rate()
#         basins_master[basin].get_volume_inundation_curve()
#     return basins_master
def load_basins(file):
    '''
    MADE THIS MYSELF!
    First value of each basin is taken out of the basins source file, and used to fill an empty dictionary.
    Afterwards a file without these rows is used to append to the original dictionary with the contours and area.
    '''
    basins_source = pd.read_csv(file, sep=',', header=0)
    basins_master = {}

    k = basins_source
    p = k.drop_duplicates(subset=['Basin_ID'], keep='first').reindex()
    q = pd.concat([basins_source, p]).drop_duplicates(keep=False)

    for row in p.index:
        tmp_basin = p.loc[row]
        basins_master[int(tmp_basin['Basin_ID'])] = Basin(str(int(tmp_basin['Basin_ID'])),
                                                          surfacearea=tmp_basin['area (m2)'], contours=[
                Contour((tmp_basin['Contour_ID']), tmp_basin['Height_min'], tmp_basin['area (m2)'])])

    for row in q.index:
        tmp_basin = q.loc[row]
        basins_master[int(tmp_basin['Basin_ID'])].Area += tmp_basin['area (m2)']
        basins_master[int(tmp_basin['Basin_ID'])].Contours.append(
            Contour((tmp_basin['Contour_ID']), tmp_basin['Height_min'], tmp_basin['area (m2)']))

    for basin in basins_master:
        basins_master[basin].get_infiltration_rate()
        basins_master[basin].get_volume_inundation_curve()
    return basins_master
def load_basin_borders(basin_master, file):
    """
    Loads in the borders of the basin as well as the areas

    The code reads a CSV file which uses a semicolon separator with info about each border, the basins it separates, and their heights.
    For each border, the code updates the BorderHeights dictionary of each basin with the mean border height between the two basins.

    Next, the code iterates over each surrounding basin of the current basin and calls up the surrounding area of each
    basin. For each surrounding basin, the code updates the SurroundingAreas dictionary of the current basin with the area of said basin.

    Input file
    :param Basin 1: The first basin that the border is adjacent to
    :param Basin 2: The second basin that the border is adjacent to
    :param Border_height_mean_lowest25% (m): the mean of the lowest 25% of the border. This is later used to calculate surface flow between basins.

    CLASS: Basin
    :param: BorderHeights
    :param: Area
    :param: SurroundingAreas
    """
    borders_source = pd.read_csv(file, sep=',', header=0)
    for row in range(0, len(borders_source)):
        tmp_border = borders_source.loc[row]
        basin_master[tmp_border['Basin 1']].BorderHeights[tmp_border['Basin 2']] = tmp_border[
            'Border_height_mean_lowest25% (m)']
        basin_master[tmp_border['Basin 2']].BorderHeights[tmp_border['Basin 1']] = tmp_border[
            'Border_height_mean_lowest25% (m)']
    for basin in basin_master:
        for surrounding_basin in basin_master[basin].BorderHeights:
            basin_master[basin].SurroundingAreas[surrounding_basin] = basin_master[surrounding_basin].Area
    return
def load_drainage_capacities(basin_master, file, small_channel, mid_channel, large_channel, low_drain, mid_drain, high_drain):
    """
    Input file
    :param Basin_ID:            [-]     ID of the basin
    :param Retention:           [m2]            Amount of retention capacity in basin
    :param Drains to basin:     [-]     The basin the current basin drains to
    :param Drainage channel:    ['yes', 'no']    Mentions whether there is a drainage channel present
    :param channel size:        ['small', 'mid', 'large']. I based this on the stream order of the drainage point.
    :param other drainage:      ['low', 'mid', 'high']
    :param To_basin_diff:       [
    :param Outlet_elevation:    [

    CLASS Basin
    :ExitChannel
    :RetentionCapacity
    :DrainsToBasin
    :DrainageChannel
    :DrainageDischarge
    :ToBasinRelativeElevation
    :OutletElevation
    """
    drainage_source = pd.read_csv(file, sep=';', header=0)
    for row in range(0, len(drainage_source)):
        tmp_basin = drainage_source.loc[row]
        basin_master[tmp_basin['Basin_ID']].ExitChannel = False                               ### What does this mean??
        basin_master[tmp_basin['Basin_ID']].RetentionCapacity = float(tmp_basin['Retention'])
        basin_master[tmp_basin['Basin_ID']].DrainsToBasin = tmp_basin['Drains to basin']

        if tmp_basin['Drainage channel'] == 'yes':
            basin_master[tmp_basin['Basin_ID']].DrainageChannel = True
            channel_size = tmp_basin['channel size']
            if channel_size == 'small':
                basin_master[tmp_basin['Basin_ID']].DrainageDischarge = small_channel
            elif channel_size == 'medium':
                basin_master[tmp_basin['Basin_ID']].DrainageDischarge = mid_channel
            elif channel_size == 'big':
                basin_master[tmp_basin['Basin_ID']].DrainageDischarge = large_channel
                basin_master[tmp_basin['Basin_ID']].ExitChannel = True                        ### What does this mean??
            else:
                print('wrong channel size inserted', channel_size)

        elif tmp_basin['Drainage channel'] == 'no':
            basin_master[tmp_basin['Basin_ID']].DrainageChannel = False
            drainage = tmp_basin['other drainage']
            if drainage == 'low':
                basin_master[tmp_basin['Basin_ID']].DrainageDischarge = low_drain
            elif drainage == 'mid':
                basin_master[tmp_basin['Basin_ID']].DrainageDischarge = mid_drain
            elif drainage == 'high':
                basin_master[tmp_basin['Basin_ID']].DrainageDischarge = high_drain
            else:
                print('wrong drainage capacity chosen', drainage)
        else:
            print('wrong drainage channel input', tmp_basin['Drainage channel'])
        #PART BELOW IS USED TO CALCULATE DRAIN DROP OFF IN ORIGINAL FILE. BUT THIS IS NOT NECESSARY, SINCE THAT CALCULATION WAS DONE WEIRDLY
        # try:
        #     basin_master[tmp_basin['Basin_ID']].ToBasinRelativeElevation = float(tmp_basin['To_basin_diff']) ### What's the point of this?
        # except ValueError:
        #     if tmp_basin['To_basin_diff'] == 'sea':
        #         basin_master[tmp_basin['Basin_ID']].ToBasinRelativeElevation = float(tmp_basin['To_basin_diff'])
        #         basin_master[tmp_basin['Basin_ID']].OutletElevation = tmp_basin['Outlet_elevation']
        #     else:
        #         print('wrong basin_diff input', tmp_basin['To_basin_diff'])

    return
def get_absolute_surrounding_water_levels(Basin, i, abs_basin_wl):

    absolute_surrounding_water_levels = {}
    for neighbor in Basin.BorderHeights:
        absolute_surrounding_water_levels[neighbor] = float(abs_basin_wl[neighbor][i])
    return absolute_surrounding_water_levels
class Contour(object):
    def __init__(self, code, min_height, surface_area): # nog beter invullen.
        self.Code               = code                     #Not called anywhere
        self.MinHeight          = min_height
        self.Area               = surface_area

class Basin(object):
    def __init__(self, name, surfacearea=0, contours=[], width=0):
        self.Area = surfacearea  # Area of the basin.
        self.Width = width  # used for storm surge calculation
        self.BorderHeights = {}  # the heights of the borders with surrounding basins.
        self.InfiltrationRate = None  # done the infiltration rate in the basin. Currently set to 2.5, but can be calculated with a different technique as well
        self.Contours = contours  # Don't know how Erik linked this to Contours
        self.DrainageDischarge = None
        self.DrainsToBasin = None
        self.VolumeToHeight = None  # done
        self.HeightToVolume = None  # done
        self.SurroundingAreas = {}
        self.RetentionCapacity = 0
        self.ExitChannel = 0
        # Most likely not necessary
        self.ToBasinRelativeElevation = None
        self.OutletElevation = None
    def get_infiltration_rate(self):

        self.InfiltrationRate = 2.5  # use until we have better method
        return
    def get_volume_inundation_curve(self):
        """

        """
        heights = [self.Contours[0].MinHeight]
        surface_areas = [self.Contours[0].Area]
        volumes = [0]
        for contour in range(1, len(self.Contours)):
            heights.append(self.Contours[contour].MinHeight)
            surface_areas.append(round(surface_areas[contour - 1] + self.Contours[contour].Area, 1))
            volumes.append(volumes[contour - 1] + (heights[contour] - heights[contour - 1]) * surface_areas[contour])
        heights.append(self.Contours[-1].MinHeight+0.25)
        surface_areas.append(surface_areas[-1])
        volumes.append(volumes[-1] + (heights[-1] - heights[-2]) * surface_areas[-1])
        curve_volume_to_height = interpolate.interp1d(volumes, heights, kind='linear')
        curve_height_to_volume = interpolate.interp1d(heights, volumes, kind='linear')
        self.VolumeToHeight = curve_volume_to_height
        self.HeightToVolume = curve_height_to_volume
        return

    def add_retention_capacity(self, capacity):

        self.RetentionCapacity = capacity
        return

def load_layers(file, basins_master):
    """
    Input file
    :param Type: ['Unprotected area', 'Protected area', 'Line of Defense']
    :param Number of flood defenses:
    :param Basin codes:
    :param Basin widths:
    code (FD1, FD2 etc.)
    :param Type: ['Land', 'Water']
    :param Incoming basin:
    :Param Basin codes:

    CLASS

    """
    #READ IN EXCEL FILE AND START WITH EMPTY FILE
    layers_source = pd.read_csv(file, sep=';', header=0)
    layers_master = {}

    #GO THROUGH EXCEL FILE AND PER ROW CHECK WHETHER LINE OF DEFENSE, PROTECTED AREA OR UNPROTECTED AREA
    for row in range(0, len(layers_source)):
        tmp_layer = layers_source.loc[row]

    #FIRST CHECK WHETHER LINE OF DEFENSE
        if tmp_layer['Type'] == 'Line of Defense':
            fd_locations = {}
    #CHECK WHETHER LAND OR WATER BASED, AND MAKES DICT WITH THE DIFFERENT FLOOD DEFENSES PER SEQUENCE NUMBER
            for defense in range(1, int(tmp_layer['Number of flood defenses']) + 1):
                code = 'FD' + str(defense) + ' '
                if tmp_layer[code + 'Type'] == 'Land':
                    basin_codes = tolist(tmp_layer[code + 'Basin codes'], ';')
                    basin_widths = tolist(tmp_layer[code + 'Basin widths'], ';')
                    basin_codes_widths_coupled = dict(zip(basin_codes, basin_widths))
                    incomingbasin = int(tmp_layer[code + 'Incoming basin'])
                    fd_locations[tmp_layer[code + 'Name']] = LineOfDefense(tmp_layer[code + 'Name'],
                                                                      tmp_layer[code + 'Type'],
                                                                      tmp_layer['Sequence'],
                                                                      basin_codes,
                                                                      basin_codes_widths_coupled,
                                                                      tmp_layer[code + 'Width'],
                                                                      tmp_layer[code + 'Height'],
                                                                      tmp_layer[code + 'Depth'],
                                                                      tmp_layer[code + 'Length'],
                                                                      incoming_basin=incomingbasin)

                elif tmp_layer[code + 'Type'] == 'Water':
                    basin_codes = tolist(tmp_layer[code + 'Basin codes'], ';')
                    basin_widths = tolist(tmp_layer[code + 'Basin widths'], ';')
                    basin_codes_widths_coupled = dict(zip(basin_codes, basin_widths))
                    fd_locations[tmp_layer[code + 'Name']] = LineOfDefense(tmp_layer[code + 'Name'],
                                                                       tmp_layer[code + 'Type'],
                                                                       tmp_layer['Sequence'],
                                                                       basin_codes,
                                                                       basin_codes_widths_coupled,
                                                                       tmp_layer[code + 'Width'],
                                                                       tmp_layer[code + 'Height'],
                                                                       tmp_layer[code + 'Depth'],
                                                                       tmp_layer[code + 'Length'],
                                                                       incoming_basin=tmp_layer[code + 'Incoming basin'])
                else:
                    print('wrong flood defense type chosen:' + str(tmp_layer['Name']))

            layers_master[tmp_layer['Sequence']] = LineOfDefense(tmp_layer['Name'],
                                                                        tmp_layer['Type'],
                                                                        tmp_layer['Sequence'],
                                                                        basin_codes,
                                                                        basin_codes_widths_coupled,
                                                                        tmp_layer[code + 'Width'],
                                                                        tmp_layer[code + 'Height'],
                                                                        tmp_layer[code + 'Depth'],
                                                                        tmp_layer[code + 'Length'],
                                                                        fdlocations=fd_locations,
                                                                        incoming_basin=tmp_layer[code + 'Incoming basin']
                                                                        )

        elif tmp_layer['Type'] == 'Protected area':
            basin_codes = tolist(tmp_layer['Basin codes'], ';')
            layers_master[tmp_layer['Sequence']] = ProtectedArea(tmp_layer['Name'], tmp_layer['Type'],
                                                                 tmp_layer['Sequence'], basin_codes)
        elif tmp_layer['Type'] == 'Unprotected area':
            basin_codes = tolist(tmp_layer['Basin codes'], ';')
            layers_master[tmp_layer['Sequence']] = UnprotectedArea(tmp_layer['Name'], tmp_layer['Type'],
                                                                   tmp_layer['Sequence'], basin_codes)
        else:
            print("Wrong layer type chosen: " + str(tmp_layer['Name']))

    # load basin widths into basin objects
    for sequence in layers_master:
        if layers_master[sequence].Type == 'Line of Defense':
            for location in layers_master[sequence].FDLocations:
                for basin in layers_master[sequence].FDLocations[location].BasinCodesWidths:
                    basins_master[basin].Width = layers_master[sequence].FDLocations[location].BasinCodesWidths[basin]

    return layers_master
class LineOfDefense(object):

    """
    example of object class
    """
    Version = "Ver 0.4"

    def __init__(self, name, typ, sequence, basin_codes, codes_widths, width, height, depth, length, slope=0.01, irribaren=1.25,  chezy=37,  fdlocations=None, incoming_basin=None):
        """

        :param name:
        :param width_unit:
        """
        self.Name = name
        self.Type = typ
        self.Sequence = sequence
        self.FDLocations = fdlocations
        self.IncomingBasin = incoming_basin     #The basin from which the storm surge comes into said basin. Basin 0 is the sea
        self.BasinCodes = basin_codes
        self.BasinCodesWidths = codes_widths

        self.Width = width                      #The width of the flood defense/coast/river

        #Land based
        self.Height = height                    #The height of the flood defense/coast/river (0 if river)
        self.Slope = slope                      #Slope of the
        self.Irribaren = irribaren              #Irribaren number for the shore


        #Water based
        self.Length = length                    # characteristic length of the inlet from sea-bay or bay-city
        self.Chezy = chezy                      #Chezy number for
        self.Depth = depth                      #Depth of river or inlet
# class LineOfDefense(object):
#
#     """
#     example of object class
#     """
#     Version = "Ver 0.4"
#
#     def __init__(self, name, typ, sequence, width_unit="m", height_unit='m+MSL', depth_unit='m', fdlocations=None, incoming_basin=None):
#         """
#
#         :param name:
#         :param width_unit:
#         """
#         self.Name = name
#         self.Type = typ
#         self.Sequence = sequence
#         self.WidthUnit = width_unit
#         self.HeightUnit = height_unit
#         self.DepthUnit = depth_unit
#         self.FDLocations = fdlocations
#         self.ScenarioProbabilities = None
#         self.Scenarios = None
#         self.IncomingBasin = incoming_basin
#         self.Active = False
#
#
#
# class LOD_Land(LineOfDefense):
#
#     """
#     example of object class
#
# https://www.google.com/search?q=TypeError%3A+__init__()+takes+1+positional+argument+but+4+were+given&rlz=1C1ONGR_nlNL973NL973&oq=TypeError%3A+__init__()+takes+1+positional+argument+but+4+were+given&aqs=chrome.0.69i59j69i64j69i60j69i58.284j0j9&sourceid=chrome&ie=UTF-8
#     """
#
#
#     def __init__(self, name, typ, sequence, basin_codes_widths, width, height, incoming_basin=None, width_unit='m', height_unit="m", slope=0.25, irribaren=1.25, used_measure=None):
#         """
#         :param name:
#         :param width:
#         :param slope:
#         :param irribaren:
#         """
#
#         LineOfDefense.__init__(self, name, typ, sequence)
#         self.IncomingBasin = incoming_basin                 #The basin from which the storm surge comes into said basin. Basin 0 is the sea
#         self.WidthUnit = width_unit                         #The unit the widths are measured in
#         self.HeightUnit = height_unit                       #The unit the heights are measured in
#         self.BasinCodesWidths = basin_codes_widths
#         self.Width = width                                  #The width of
#         self.Height = height                                #The height of the LOD (also necessary if no measure present)
#         self.Slope = slope                                  #The slope at the LOD location
#         self.Irribaren = irribaren                          #The irribaren number at the LOD
#         self.UsedMeasure = used_measure
#
#     # def __init__(self, basin_codes_widths, width, height, slope=0.25, irribaren=1.25, incoming_basin=None):
#     #     """
#     #     :param name:
#     #     :param width:
#     #     :param slope:
#     #     :param irribaren:
#     #     """
#     #
#     #   #  LineOfDefense.__init__(self, name, typ, sequence, incoming_basin, width_unit, height_unit)
#     #     LineOfDefense.__init__(self, name, typ, sequence)
#     #     self.BasinCodesWidths = basin_codes_widths
#     #     self.Width = width
#     #     self.Height = height
#     #     self.Slope = slope
#     #     self.Irribaren = irribaren
#     #     self.IncomingBasin = incoming_basin
#     #
#     #     #super().__init__(self)
#
# class LOD_Water(LineOfDefense):
#     """
#        example of object class
#        """
#
#     def __init__(self, name, typ, sequence, basincodes,basin_codes_widths, width, depth, length, incoming_basin=None, width_unit='m', height_unit="m", depth_unit="m", length_unit='m', slope=0.25, irribaren=1.25, chezy=37, used_measure=None):
#         """
#
#         :param name:
#         :param width:
#         :param depth:
#         :param length:
#         :param depth_unit:
#         :param slope:
#         :param irribaren:
#         :param chezy:
#         """
#
#         LineOfDefense.__init__(self, name, typ, sequence)
#         self.BasinCodes = basincodes
#         self.BasinCodesWidths= basin_codes_widths
#         self.Width = width
#         self.Depth = depth
#         self.Length = length  # characteristic length of the inlet from sea-bay or bay-city
#         self.IncomingBasin = incoming_basin
#         self.WidthUnit = width_unit
#         self.HeightUnit = height_unit
#         self.DepthUnit = depth_unit
#         self.LengthUnit = length_unit
#         self.Slope = slope
#         self.Irribaren = irribaren
#         self.Chezy = chezy
#         self.UsedMeasure = used_measure
class ProtectedArea(object):
    """
    example of object class
    """
    Version = "0.0"

    def __init__(self, name, typ, sequence, basincodes, list_surface_area=[], surface_area_zero=0):
        """

        :param name:
        :param list_surface_area:
        :param surface_area_zero:
        """

        self.Name = name
        self.Type = typ
        self.Sequence = sequence
        self.BasinCodes = basincodes
        self.ListSurfaceArea = list_surface_area
        self.SurfaceAreaZero = surface_area_zero
        self.PDLocations = []

    def create_list_surface_area(self, dict_basin):
        """

        :param dict_basin:
        :return:
        """

        self.ListSurfaceArea = []
        for contour in range(1, 31):
            extra_surf = []
            for basin in dict_basin:
                try:
                    extra_surf.append(dict_basin[basin].Contours[contour].SurfaceArea)
                except ValueError:
                    continue

            self.ListSurfaceArea.append(self.SurfaceAreaZero + np.sum(extra_surf))

        for i in range(1, len(self.ListSurfaceArea)):
            if self.ListSurfaceArea[i] == self.SurfaceAreaZero:
                self.ListSurfaceArea[i] = self.ListSurfaceArea[i - 1]

        return self.ListSurfaceArea
class UnprotectedArea(ProtectedArea):

    def __init__(self, name, typ, sequence, basincodes):

        ProtectedArea.__init__(self, name, typ, sequence, basincodes, list_surface_area=[], surface_area_zero=0)

def load_hydraulic_conditions(file_surge, file_rain, timestep):

    hydraulic_conditions_master = {}
    surge_master = pd.read_csv(file_surge, sep=',', header=0)
    rain_master = pd.read_csv(file_rain, sep=',', header=0)
    surge_master_dropped = surge_master.dropna()
    rain_master_dropped = rain_master.dropna()
    timestep_model = timestep  # hours

    for i in range(0, len(surge_master_dropped)):
        for j in range(0, len(rain_master_dropped)):
            tmp_surge = surge_master_dropped.loc[i]
            tmp_rain = rain_master_dropped.loc[j]
            key = str(int(tmp_surge['return period storm'])) + ',' + str(int(tmp_rain['return period rain']))             # name of the storm/rain combination, e.g. '100;10'
            length_series = int(tmp_surge['storm duration']/timestep_model + 1)
            rain_series = [int(tmp_rain['Rain intensity_1'])] * length_series

            hydraulic_conditions_master[key] = HydraulicBoundaryConditions(
                tmp_surge['return period storm'], tmp_rain['return period rain'],
                rain_series, tmp_surge['maximum surge'], tmp_rain['Rain intensity_3'],
                tmp_surge['storm duration'], timestep_model,
                tmp_surge['normal tidal amplitude'], tmp_surge['MSL'])

    return hydraulic_conditions_master
class HydraulicBoundaryConditions (object):
    def __init__(self, storm, rain, rain_series, max_surge, max_rain, storm_duration, timestep_model, tide_amplitude, msl):
        self.ReturnPeriodStorm = storm
        self.ReturnPeriodRain = rain
        self.SurgeSeries = create_surge_series(max_surge, storm_duration, tide_amplitude, timestep_model, msl)
        self.RainSeries = rain_series
        self.MaxLevel = max_surge
        self.MaxRainIntensity = max_rain
        self.StormDuration = storm_duration
        self.SurgeLayerCount = 1
        self.InlandSurge = {}
        self.TideAmplitude = tide_amplitude
        self.MeanSeaLevel = msl

    def __str__(self):
        return ('    STORM SURGE\n'
                '   Return period storm:    {ReturnPeriodStorm}\n'
                '   Maximum surge level:    {MaxLevel}\n'
                '   Outside wave height:    {WaveHeight}\n'
                '   Storm duration:         {StormDuration}\n'
                '   Maximum wind velocity:  {WindVelocity}\n'
                '   \n'
                '     RAIN\n'
                '   Return period rain:     {ReturnPeriodRain}\n'
                '   Maximum rain intensity: {MaxRainIntensity}\n'
                '\n'
                ''
                ).format(**self.__dict__)

