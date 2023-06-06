"""
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

"""
import pandas as pd
import numpy as np
from scipy import interpolate


def tolist(input_string, splitby, export_type=int):
    """Convert imported string into a list.

    This function is used to divide list that has been imported as a string into an actual list. For example, if the
    import is '1;2;3;4', the output is [1,2,3,4]

    Parameters
    ----------
    input_string : str
        string value, consisting of divided numbers
    splitby : str
        punctuation mark used to divide the numbers
    export_type: data-type, default=int

    Returns
    -------
    list : list
        list of values
    """

    if type(input_string) == str:
        return list(map(export_type, input_string.split(splitby)))
    else:
        str_variable = str(int(input_string))
        return list(map(export_type, str_variable.split(splitby)))


def pad_or_truncate(some_list, target_len, values=0):
    """Shorten or extend a list based on target length.

    This function shortens or extends a list. When it is extended, extra values are added. This can be helpful to
    simulate what happens after a storm has passed.

    Parameters
    ----------
    some_list : list
        The original list of values
    target_len : int
        The wanted length of the returned list
    values : int, default=0
        Value to pad onto list

    Returns
    -------
    list
        a list of values with length of target_len
    """
    return some_list[:target_len] + [values] * (target_len - len(some_list))


def create_surge_series(
    max_surge, storm_duration, tidal_amplitude, TIMESTEP, mean_sea_level
):
    """Creates time series of surge

    This function makes a time series of the surge, consisting of tide and a surge. The surge is schematized as a
    half sine. the tide is schematized to have its peak halfway through the storm duration (when surge is max)

    Parameters
    ----------
    max_surge : int
        the maximum additional surge of the water level because of the storm [m]
    storm_duration : int
        duration of the storm [hours]
    tidal_amplitude : int
        normal amplitude of the tide [m]
    TIMESTEP: int
        length of one timestep [hours]
    mean_sea_level: int
        absolute height of the mean sea level [m]

    Returns
    -------
    list
        time series of water levels [m+MSL]
    """

    TIDE_INTERVAL = 12.417  # 12 hours and 25 minutes
    time1 = np.linspace(
        0, storm_duration, int(storm_duration / TIMESTEP + 1), endpoint=True
    )
    tide = (
        tidal_amplitude
        * np.cos(np.pi / (TIDE_INTERVAL / 2) * time1 - 0.5 * storm_duration)
        + mean_sea_level
    )
    surge = max_surge * np.sin(np.pi / storm_duration * time1) ** 5
    total_storm_surge_series = tide + surge

    return list(total_storm_surge_series)


def load_basins(file):
    """Load in the basins.

    This functions loads in the Basin class with data
    First value of each basin is taken out of the basins source file, and used to fill an empty dictionary.
    Afterwards a file without these rows is used to append to the original dictionary with the contours and contour_area.

    Parameters
    ----------
    file : str, path object or file-like object
        csv file containing the following data: 'Basin_ID', 'contour_area (m2)', 'Height_min'


    Returns
    -------
    basins_master : dict
        Dictionary with different Basin objects in it

    See Also
    --------
    load_basin_borders, load_drainage_capacities

    Notes
    -----
    MADE THIS MYSELF!
    """
    basins_source = pd.read_csv(file, sep=",", header=0)
    basins_master = {}

    k = basins_source
    p = k.drop_duplicates(subset=["Basin_ID"], keep="first").reindex()
    q = pd.concat([basins_source, p]).drop_duplicates(keep=False)

    for row in p.index:
        tmp_basin = p.loc[row]
        basins_master[int(tmp_basin["Basin_ID"])] = Basin(
            str(int(tmp_basin["Basin_ID"])),
            surfacearea=tmp_basin["area (m2)"],
            contours=[
                Contour(
                    (tmp_basin["Contour_ID"]),
                    tmp_basin["Height_min"],
                    tmp_basin["area (m2)"],
                )
            ],
        )

    for row in q.index:
        tmp_basin = q.loc[row]
        basins_master[int(tmp_basin["Basin_ID"])].area += tmp_basin["area (m2)"]
        basins_master[int(tmp_basin["Basin_ID"])].contours.append(
            Contour(
                (tmp_basin["Contour_ID"]),
                tmp_basin["Height_min"],
                tmp_basin["area (m2)"],
            )
        )

    for basin in basins_master:
        basins_master[basin].get_infiltration_rate()
        basins_master[basin].get_volume_inundation_curve()
    return basins_master


def load_basin_borders(basins_master, file):
    """Loads in the borders of the basin as well as the areas of neighboring basins.

    The code reads a CSV file which uses a semicolon separator with info about each border, the basins it separates, and their heights.
    For each border, the code updates the border_heights dictionary of each basin with the mean border height between the two basins.
    Next, the code iterates over each surrounding basin of the current basin and calls up the surrounding contour_area of each
    basin. For each surrounding basin, the code updates the surrounding_areas dictionary of the current basin with the contour_area of said basin.

    Parameters:
    ----------
    basins_master : dict
        The basins_masters has been instantiated with the load_basins function and is added onto with this one.
        The following attributes are added onto with this function: border_heights, contour_area, surrounding_areas
    file : str, path object or file-like object
        The csv file has the following 'Basin 1', 'Basin 2', 'Border_height_mean_lowest25% (m)'

    Returns
    -------


    See Also
    --------
    load_basins, load_drainage capacities

    """
    borders_source = pd.read_csv(file, sep=",", header=0)
    for row in range(0, len(borders_source)):
        tmp_border = borders_source.loc[row]
        basins_master[tmp_border["Basin 1"]].border_heights[
            tmp_border["Basin 2"]
        ] = tmp_border["Border_height_mean_lowest25% (m)"]
        if tmp_border["Basin 2"] != 0:
            basins_master[tmp_border["Basin 2"]].border_heights[
                tmp_border["Basin 1"]
            ] = tmp_border["Border_height_mean_lowest25% (m)"]
    for basin in basins_master:
        for surrounding_basin in basins_master[basin].border_heights:
            if surrounding_basin == 0:
                basins_master[basin].surrounding_areas[surrounding_basin] = 10**20
            else:
                basins_master[basin].surrounding_areas[surrounding_basin] = basins_master[
                    surrounding_basin
                ].area
    return


def load_drainage_capacities(
    basin_master,
    file,
    small_channel,
    mid_channel,
    large_channel,
    low_drain,
    mid_drain,
    high_drain,
):
    """Loads in the drainage capacities for the Basin object class.

    Adds onto basins_master
    Drainage for different sized channels and different sized other drainage options are specified, which is used to determine the drainage discharge.
    Additionally, the exit channel, retention capacity, basin it drains to as well as the size of the drainage channel are specified.

    Parameters
    ----------
    basin_master : dict
        Adds the following attributes: retention_capacity, drains_to_basin, drainage_channel, drainage_discharge
    file : str, path object or file-like object
        This file has the following: 'Basin_ID', 'Retention' [m2], 'Drains to basin', 'Drainage channel' ['yes', 'no'], 'channel size' ['small', 'mid', 'large'], 'other drainage' ['low', 'mid', 'high'], 'To_basin_diff', 'Outlet_elevation'.
    small_channel : int
        Drainage for a small channel as defined in 'channel size' [m3/s]
    mid_channel : int
        Drainage for a mid-sized channel as defined in 'channel size' [m3/s]
    large_channel : int
        Drainage for a large channel as defined in 'channel size' [m3/s]
    low_drain : int
        Drainage for low drainage as defined in 'other drainage' [m3/s]
    mid_drain : int
        Drainage for medium drainage as defined in 'other drainage' [m3/s]
    high_drain : int
        Drainage for high drainage as defined in 'other drainage' [m3/s]

    Returns
    -------

    See Also
    --------
    load_basins, load_basin_borders
    """
    drainage_source = pd.read_csv(file, sep=";", header=0)
    for row in range(0, len(drainage_source)):
        tmp_basin = drainage_source.loc[row]

        basin_master[tmp_basin["Basin_ID"]].retention_capacity = float(
            tmp_basin["Retention"]
        )
        basin_master[tmp_basin["Basin_ID"]].drains_to_basin = tmp_basin[
            "Drains to basin"
        ]

        if tmp_basin["Drainage channel"] == "yes":
            # basin_master[tmp_basin["Basin_ID"]].DrainageChannel = True
            channel_size = tmp_basin["channel size"]
            if channel_size == "small":
                basin_master[tmp_basin["Basin_ID"]].drainage_discharge = small_channel
            elif channel_size == "medium":
                basin_master[tmp_basin["Basin_ID"]].drainage_discharge = mid_channel
            elif channel_size == "big":
                basin_master[tmp_basin["Basin_ID"]].drainage_discharge = large_channel
            else:
                print("wrong channel size inserted", channel_size)

        elif tmp_basin["Drainage channel"] == "no":
            # basin_master[tmp_basin["Basin_ID"]].DrainageChannel = False
            drainage = tmp_basin["other drainage"]
            if drainage == "low":
                basin_master[tmp_basin["Basin_ID"]].drainage_discharge = low_drain
            elif drainage == "mid":
                basin_master[tmp_basin["Basin_ID"]].drainage_discharge = mid_drain
            elif drainage == "high":
                basin_master[tmp_basin["Basin_ID"]].drainage_discharge = high_drain
            else:
                print("wrong drainage capacity chosen", drainage)
        else:
            print("wrong drainage channel input", tmp_basin["Drainage channel"])

    return


def get_absolute_surrounding_water_levels(Basin, i, abs_basin_wl, outside_water_level):
    """Finds the absolute surrounding water levels for a specified basin at a certain timestep.

    The function finds the neighboring basins based on the border heights, and then finds their water level [m] for the specified timestep i.
    It puts these values into a dictionary and returns it.

    Parameters
    ----------
    Basin : object
        Basin to find the surrounding water levels for, from the basins_master file.
    i : int
        Timestep for which calculation is done.
    abs_basin_wl : float
        The absolute water level in the basin [m]
    outside_water_level : float
        storm surge forcing [m]

    Returns
    -------
    absolute_surrounding_water_levels : dict
        dictionary with for each neighboring basin the absolute water level {neighboring basin: wl [m], ...}

    """

    absolute_surrounding_water_levels = {}

    for neighbor in Basin.border_heights:
        if neighbor !=0:
            absolute_surrounding_water_levels[neighbor] = float(abs_basin_wl[neighbor][i])
        else:
            absolute_surrounding_water_levels[neighbor] = outside_water_level[i]
    return absolute_surrounding_water_levels


class Contour(object):
    """
    Object that holds a contour line for a basin.
    This information is loaded in the load_basins function
    ...
    Parameters
    ----------
    code : int
        ID of the Contour
    min_height : float
        the lower bound of the ID'd contour
    surface_area : float
        the area of the ID'd contour

    Attributes
    ----------
    code : int
        the ID of the contour
    min_height : float
        the lower bound of the contour. [m]
    contour_area : float
        the area of a contour [m2]

    See Also
    --------
    load_basin

    """

    def __init__(self, code, min_height, surface_area):
        self.code = code
        self.min_height = min_height
        self.contour_area = surface_area


class Basin(object):
    """class that holds the information for a basin object
    ...
    Parameters
    ----------
    surfacearea : float
        surface area of the Basin [m2]
    contours : iterable
        list with contours present within the basin
    width : float
        width of the basin

    Attributes
    ----------
    area : float
        the total area of the basin [m2]
    width : float
        the width of the flood defenses within the basin [m]
    border_heights : dict
        dictionary of neighboring basins and their border heights
    infiltration_rate : float
        infiltration rate of the basin [mm/hour]
    contours : list
        list of contours that are present within the basin
    drainage_discharge : float
        The maximum discharge possible for the basin. [m3/s]
    drains_to_basin : int
        the basin to which the current basin drains
    volume_to_height :
        method to convert volume of water in the basin to absolute water level [m3] to [m]
    height_to_volume :
        method to convert absolute water level in the basin to a volume of water [m] to [m3]
    surrounding_areas : dict
        dictionary of the neighboring basins and their surface areas [m2]
    retention_capacity : float
        retention capacity of the basin [m2]

    See Also
    --------
    Class Contour
    load_basin, load_drainage_capacities, load_basin_borders, load_layers
    """

    def __init__(self, name, surfacearea=0, contours=[], width=0):
        self.area = surfacearea
        self.width = width
        self.border_heights = {}
        self.infiltration_rate = None
        self.contours = contours
        self.drainage_discharge = None
        self.drains_to_basin = None
        self.volume_to_height = None
        self.height_to_volume = None
        self.surrounding_areas = {}
        self.retention_capacity = 0
        # self.exit_channel = 0
        # # Most likely not necessary
        # self.ToBasinRelativeElevation = None
        # self.OutletElevation = None

    def get_infiltration_rate(self):
        """sets the infiltration rate for the Basin

        The infiltration rate for sand is 20-30 [mm/hour]
        Can potentially also be done with a better method
        ...

        References
        ----------
        https://www.fao.org/3/S8684E/s8684e0a.htm

        """
        self.infiltration_rate = 20.5
        return

    def get_volume_inundation_curve(self):
        """Used to add attributes volume_to_height and height_to_volume

        The function uses the heights of the contours as well as their areas to construct a volume-depth curve.
        This curve is then used to calculate a Basin's height based on the water volume inside it and vice versa
        """
        heights = [self.contours[0].min_height]
        surface_areas = [self.contours[0].contour_area]
        volumes = [0]
        for contour in range(1, len(self.contours)):
            heights.append(self.contours[contour].min_height)
            surface_areas.append(
                round(
                    surface_areas[contour - 1] + self.contours[contour].contour_area, 1
                )
            )
            volumes.append(
                volumes[contour - 1]
                + (heights[contour] - heights[contour - 1]) * surface_areas[contour]
            )
        heights.append(self.contours[-1].min_height + 0.25)
        surface_areas.append(surface_areas[-1])
        volumes.append(volumes[-1] + (heights[-1] - heights[-2]) * surface_areas[-1])
        curve_volume_to_height = interpolate.interp1d(volumes, heights, kind="linear")
        curve_height_to_volume = interpolate.interp1d(heights, volumes, kind="linear")
        self.volume_to_height = curve_volume_to_height
        self.height_to_volume = curve_height_to_volume
        return

    def add_retention_capacity(self, capacity):
        """used to set the retention capacity for the basin.


        Parameters
        ----------
        capacity : float
            retention capacity [m2]

        Returns
        -------

        """
        self.retention_capacity = capacity
        return


def load_layers(file, basins_master):
    """loads in the layers of defense and protected/unprotected areas

    Parameters
    ----------
    file : str, path object or file-like object
        the csv file has the following information: 'Type', 'Name', 'Sequence', 'Number of flood defenses', and various data on width, height, length and depth.
    basins_master : dict
        the dictionary with the different basins inside it

    See Also
    --------
    Class LineOfDefense, ProtectedArea, UnprotectedArea

    Notes
    -----
    The csv file has to be set up in a specific manner. The general "Type" field can only have one of three options: ['Unprotected area', 'Protected area', 'Line of Defense'].
    The description of the variables for each flood defense should be preceded by its code, so for example 'FD1 Width', where the 1 is the order in the number of flood defenses
    The Code + Type field (FD1 Type, FD2 Type) should be one of two string values: ['Land', 'Water']

    """
    # READ IN EXCEL FILE AND START WITH EMPTY FILE
    layers_source = pd.read_csv(file, sep=",", header=0)
    layers_master = {}

    # GO THROUGH EXCEL FILE AND PER ROW CHECK WHETHER LINE OF DEFENSE, PROTECTED AREA OR UNPROTECTED AREA
    for row in range(0, len(layers_source)):
        tmp_layer = layers_source.loc[row]

        # FIRST CHECK WHETHER LINE OF DEFENSE
        if tmp_layer["Type"] == "Line of Defense":
            fd_locations = {}
            # CHECK WHETHER LAND OR WATER BASED, AND MAKES DICT WITH THE DIFFERENT FLOOD DEFENSES PER SEQUENCE NUMBER
            for defense in range(1, int(tmp_layer["Number of flood defenses"]) + 1):
                code = "FD" + str(defense) + " "
                if tmp_layer[code + "Type"] == "Land":
                    basin_codes = tolist(tmp_layer[code + "Basin codes"], ";")
                    basin_widths = tolist(tmp_layer[code + "Basin widths"], ";")
                    basin_codes_widths_coupled = dict(zip(basin_codes, basin_widths))
                    incomingbasin = int(tmp_layer[code + "Incoming basin"])
                    fd_locations[tmp_layer[code + "Name"]] = LineOfDefense(
                        tmp_layer[code + "Name"],
                        tmp_layer[code + "Type"],
                        tmp_layer["Sequence"],
                        basin_codes,
                        basin_codes_widths_coupled,
                        tmp_layer[code + "Width"],
                        tmp_layer[code + "Height"],
                        tmp_layer[code + "Depth"],
                        tmp_layer[code + "Length"],
                        incoming_basin=incomingbasin,
                    )

                elif tmp_layer[code + "Type"] == "Water":
                    basin_codes = tolist(tmp_layer[code + "Basin codes"], ";")
                    basin_widths = tolist(tmp_layer[code + "Basin widths"], ";")
                    basin_codes_widths_coupled = dict(zip(basin_codes, basin_widths))
                    fd_locations[tmp_layer[code + "Name"]] = LineOfDefense(
                        tmp_layer[code + "Name"],
                        tmp_layer[code + "Type"],
                        tmp_layer["Sequence"],
                        basin_codes,
                        basin_codes_widths_coupled,
                        tmp_layer[code + "Width"],
                        tmp_layer[code + "Height"],
                        tmp_layer[code + "Depth"],
                        tmp_layer[code + "Length"],
                        incoming_basin=tmp_layer[code + "Incoming basin"],
                    )
                else:
                    print("wrong flood defense type chosen:" + str(tmp_layer["Name"]))

            layers_master[tmp_layer["Sequence"]] = LineOfDefense(
                tmp_layer["Name"],
                tmp_layer["Type"],
                tmp_layer["Sequence"],
                basin_codes,
                basin_codes_widths_coupled,
                tmp_layer[code + "Width"],
                tmp_layer[code + "Height"],
                tmp_layer[code + "Depth"],
                tmp_layer[code + "Length"],
                fdlocations=fd_locations,
                incoming_basin=tmp_layer[code + "Incoming basin"],
            )

        elif tmp_layer["Type"] == "Protected area":
            basin_codes = tolist(tmp_layer["Basin codes"], ";")
            layers_master[tmp_layer["Sequence"]] = ProtectedArea(
                tmp_layer["Name"], tmp_layer["Type"], tmp_layer["Sequence"], basin_codes
            )
        elif tmp_layer["Type"] == "Unprotected area":
            basin_codes = tolist(tmp_layer["Basin codes"], ";")
            layers_master[tmp_layer["Sequence"]] = UnprotectedArea(
                tmp_layer["Name"], tmp_layer["Type"], tmp_layer["Sequence"], basin_codes
            )
        else:
            print("Wrong layer type chosen: " + str(tmp_layer["Name"]))

    # load basin widths into basin objects
    for sequence in layers_master:
        if layers_master[sequence].type == "Line of Defense":
            for location in layers_master[sequence].location_flood_defenses:
                for basin in (
                    layers_master[sequence]
                    .location_flood_defenses[location]
                    .basin_codes_widths
                ):
                    basins_master[basin].width = (
                        layers_master[sequence]
                        .location_flood_defenses[location]
                        .basin_codes_widths[basin]
                    )

    return layers_master


class LineOfDefense(object):
    """There are two types of lines of defense, water based and land based.

    Parameters
    ----------
    name : str
        name of the line of defense. This can be used as a description
    typ : str
        either ['Water'] or ['Land']
    sequence : int
        the place in the coastal defense sequence for the line of defense.
    basin_codes : list
        the basin codes where the line of defense is present.
    codes_widths : dict
        the coupled basin codes with their widths for the line of defense
    width : float
        The width of the flood defense/coast/river [m]
    height : float
        The height of the flood defense/coast/river [m]. height of a river is 0
    depth : float
        The depth of the water based flood defense [m]
    length : float
        The length of the water based flood defense [m]
    slope : float
        The bed slope of the water based flood defense
    iribarren : float
        The iribarren number for the coast
    chezy : float
        the Chezy number for the open channel [m0.5/s]
    fdlocations : dict
        a dictionary containing all flood defenses and their position in sequences and basins.
    incoming_basin : int
        The basin from which the storm surge comes into said basin. Basin 0 is the sea

    Attributes
    ----------
    name : str
        name of the line of defense. This can be used as a description
    type : str
        either ['Water'] or ['Land']
    sequence : int
        the place in the coastal defense sequence for the line of defense.
    location_flood_defenses : dict
        a dictionary containing all flood defenses and their position in sequences and basins.
    incoming_basin : int
        The basin from which the storm surge comes into said basin. Basin 0 is the sea
    basin_codes : list
        the basin codes where the line of defense is present.
    basin_codes_widths : dict
        the coupled basin codes with their widths for the line of defense
    width : float
        The width of the flood defense/coast/river [m]
    height : float
        The height of the flood defense/coast/river [m]. height of a river is 0
    iribarren : float
        The iribarren number for the coast
    length : float
        The length of the water based flood defense [m]
    slope : float
        The bed slope of the water based flood defense
    chezy : float
        the Chezy number for the open channel [m0.5/s]
    depth : float
        The depth of the water based flood defense [m]
    """

    def __init__(
        self,
        name,
        typ,
        sequence,
        basin_codes,
        codes_widths,
        width,
        height,
        depth,
        length,
        slope=0.01,
        iribarren=1.25,
        chezy=37,
        fdlocations=None,
        incoming_basin=None,
    ):
        self.name = name
        self.type = typ
        self.sequence = sequence
        self.location_flood_defenses = fdlocations
        self.incoming_basin = incoming_basin
        self.basin_codes = basin_codes
        self.basin_codes_widths = codes_widths

        self.width = width

        # Land based
        self.height = height
        self.irribaren = iribarren

        # Water based
        self.length = length
        self.slope = slope
        self.chezy = chezy
        self.depth = depth


class ProtectedArea(object):
    """
    A Protected Area is an area behind a line of defense.

    Parameters
    ----------
    name : str
        name of protected area
    typ : str
        type of area [Protected area, Unprotected area]
    sequence : int
        the place of the area in the coastal defense (starts at 0)
    basincodes : list
        list with the codes of the basins within the area

    Attributes
    ----------
    name : str
        name of protected area
    type : str
        type of area [Protected area, Unprotected area]
    sequence : int
        the place of the area in the coastal defense (starts at 0)
    basincodes : list
        list with the codes of the basins within the area
    protected_locations : list
        list of the protected locations         chopping block
    """

    def __init__(self, name, typ, sequence, basincodes):
        self.name = name
        self.type = typ
        self.sequence = sequence
        self.basin_codes = basincodes
        # self.surface_area_list = list_surface_area
        # self.surface_area_zero = surface_area_zero
        self.protected_locations = []

    # def create_list_surface_area(self, dict_basin):
    #     """Appends to the surface area list of the protected area.
    #
    #     The function calculates the areas
    #
    #     ...
    #     Parameters
    #     ----------
    #     dict_basin : dict
    #         dictionary with all the basins in it
    #
    #     Returns
    #     -------
    #     self.surface_area_list : list
    #         list with the surface areas for each basin within the area
    #     """
    #
    #     self.surface_area_list = []
    #     for contour in range(1, 31):
    #         extra_surf = []
    #         for basin in dict_basin:
    #             try:
    #                 extra_surf.append(dict_basin[basin].contours[contour].SurfaceArea)
    #             except ValueError:
    #                 continue
    #
    #         self.surface_area_list.append(self.surface_area_zero + np.sum(extra_surf))
    #
    #     for i in range(1, len(self.surface_area_list)):
    #         if self.surface_area_list[i] == self.surface_area_zero:
    #             self.surface_area_list[i] = self.surface_area_list[i - 1]
    #
    #     return self.surface_area_list


class UnprotectedArea(ProtectedArea):
    """An unprotected Area is a subclass of a protected area that is not behind a line of defense
    Parameters
    ----------
    name : str
        name of protected area
    typ : str
        type of area [Protected area, Unprotected area]
    sequence : int
        the place of the area in the coastal defense (starts at 0)
    basincodes : list
        list with the codes of the basins within the area

    Attributes
    ----------
    name : str
        name of protected area
    type : str
        type of area [Protected area, Unprotected area]
    sequence : int
        the place of the area in the coastal defense (starts at 0)
    basincodes : list
        list with the codes of the basins within the area
    """

    def __init__(self, name, typ, sequence, basincodes):
        ProtectedArea.__init__(
            self,
            name,
            typ,
            sequence,
            basincodes,
        )


def load_hydraulic_conditions(file_surge, file_rain, TIMESTEP):
    """Loads in the hydraulic boundary conditions from a surge file and a rain file.

    This function loads in two csv's: one with information about the tides and storm surge, and one with rain.
    These files are loaded into the HydraulicBoundaryConditions class, and are returned as a dict.

    Parameters
    ----------
    file_surge : str, path object or file-like object
        This file has the following: 'return period storm', 'storm duration', 'maximum surge', 'normal tidal amplitude', 'MSL' .
    file_rain : str, path object or file-like object
        This file has the following: 'return period rain', 'Rain_intensity_1', 'Rain intensity_3' .
    TIMESTEP : float
        a constant. This represents the timestep taken in the model [hours]

    Returns
    -------
    hydraulic_conditions_master : dict

    References
    ----------
    Rapid screening and evaluation of flood risk reduction strategies: Exploratory study on the use of the FLORES modelling approach for World Bank projects - Erik van Berchum


    """

    hydraulic_conditions_master = {}
    surge_master = pd.read_csv(file_surge, sep=",", header=0)
    rain_master = pd.read_csv(file_rain, sep=",", header=0)
    surge_master_dropped = surge_master.dropna()
    rain_master_dropped = rain_master.dropna()
    timestep_model = TIMESTEP

    for i in range(0, len(surge_master_dropped)):
        for j in range(0, len(rain_master_dropped)):
            tmp_surge = surge_master_dropped.loc[i]
            tmp_rain = rain_master_dropped.loc[j]
            key = (
                str(int(tmp_surge["return period storm"]))
                + ","
                + str(int(tmp_rain["return period rain"]))
            )  # name of the storm/rain combination, e.g. '100,10'
            length_series = int(tmp_surge["storm duration"] / timestep_model + 1)
            rain_series = [int(tmp_rain["Rain intensity_1"])] * length_series

            hydraulic_conditions_master[key] = HydraulicBoundaryConditions(
                tmp_surge["return period storm"],
                tmp_rain["return period rain"],
                rain_series,
                tmp_surge["maximum surge"],
                tmp_rain["Rain intensity_3"],
                tmp_surge["storm duration"],
                timestep_model,
                tmp_surge["normal tidal amplitude"],
                tmp_surge["MSL"],
            )

    return hydraulic_conditions_master


class HydraulicBoundaryConditions(object):
    """

    Parameters
    ----------
    storm : int
        return period for the storm [years]
    rain : int
        return period for the rainfall [years]
    rain_series : list
        list with rainfall for the storm duration [mm/timestep]
    max_surge : float
        the maximum skew surge for the return period [m]
    max_rain : float
        the maximal rainfall for the return period [mm/hour]
    storm_duration : int
        the duration of the storm. [hours]
    timestep_model : float
        the timestep of the model
    tide_amplitude : float
        the amplitude of the tide. This is added on top of the maximum surge to get a surge series [m]
    msl : float
        mean sea level. [m]

    Attributes
    ----------
    return_period_storm : int
        return period for the storm [years]
    return_period_rain : int
        return period for the rainfall [years]
    surge_series : list
        list with the water level at each timestep of the storm's duration. This is based on the tide, max surge and MSL.
    rain_series : list
        list with rainfall for each timestep of the storm's duration [mm/timestep]
    max_surge : float
        the maximum skew surge for the return period [m]
    max_rain_intensity : float
        the maximal rainfall for the return period [mm/hour]
    storm_duration : int
        the duration of the storm. [hours]
    tidal_amplitude : float
        the amplitude of the tide. This is added on top of the maximum surge to get a surge series [m]
    mean_sea_level : float
        mean sea level. [m]

    """

    def __init__(
        self,
        storm,
        rain,
        rain_series,
        max_surge,
        max_rain,
        storm_duration,
        timestep_model,
        tide_amplitude,
        msl,
    ):
        self.return_period_storm = storm
        self.return_period_rain = rain
        self.surge_series = create_surge_series(
            max_surge, storm_duration, tide_amplitude, timestep_model, msl
        )
        self.rain_series = rain_series
        self.max_surge = max_surge
        self.max_rain_intensity = max_rain
        self.storm_duration = storm_duration
        # self.surge_layer_count = 1
        # self.inland_surge = {}
        self.tidal_amplitude = tide_amplitude
        self.mean_sea_level = msl

    def __str__(self):
        return (
            "    STORM SURGE\n"
            "   Return period storm:    {return_period_storm}\n"
            "   Maximum surge level:    {max_surge}\n"
            "   Outside wave height:    {WaveHeight}\n"
            "   Storm duration:         {storm_duration}\n"
            "   Maximum wind velocity:  {WindVelocity}\n"
            "   \n"
            "     RAIN\n"
            "   Return period rain:     {return_period_rain}\n"
            "   Maximum rain intensity: {max_rain_intensity}\n"
            "\n"
            ""
        ).format(**self.__dict__)
