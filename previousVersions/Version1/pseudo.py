for each timestep
    for each basin
        calculate inflow rain at time
        if basin is behind line of defense:
            calculate inflow storm surge
        calculate inflow from other basins

        calculate outflow into retention
        calculate outflow into drainage
        calculate outflow infiltration
        if water level > elevation drainage basin border
            calculate flow between basins

        calculate water level at t+1

class Basin
    area
    lowest point (lowest 25%)
    Lines of defense
    Retention area
    volumetric curve = contours
    drainage
    neighbors etc.
    infiltration rate?

class Line of Defense_water

