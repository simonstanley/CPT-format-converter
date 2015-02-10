"""
Convert NETCDF/grib/pp files into the format required for the IRI climate 
prediction tool (CPT). 
Important, this converter is currently for gridded data sets only. CPT also 
takes station data, index data and probability data, this does not convert 
these.

"""
import iris
import numpy
import datetime
import os

cpt_field_name_dict = {'precipitation_flux' : 'prcp'}

def check_filepath(filepath):
    """
    Check if file exists, if so, ask the user if they wish to overwrite it.
    
    Args:
    
    * filepath: string
        Full path to the file to check
    
    """
    if os.path.isfile(filepath):
        response = raw_input("This file already exists, would you like to "\
                             "overwrite?")
        if response.lower() in ['y', 'yes']:
            os.remove(filepath)
        elif response.lower() in ['n', 'no']:
            raise UserWarning('Program stopped.')
        else:
            print "Please respond yes or no"
            check_filepath(filepath)

def cube_time_converter(time, time_unit):
    """
    Convert time between datetime object and number value (standard format in
    an iris cube).

    Args:
    
    * time: datetime, or float, int
        Which ever type is given, it is converted to the other.
    
    * time_unit: iris.unit.Unit
        This describes the number so it can be converted. For example, the time
        unit may be 'hours since 1970'.

    Returns:
        datetime or float

    """
    assert type(time_unit) == iris.unit.Unit, 'time_unit must be iris.unit.Unit, not %s' % type(time_unit)
    if type(time) == datetime.datetime:
        converted_time = time_unit.date2num(time)
    else:
        converted_time = time_unit.num2date(time)
        # Convert to proper datetime so operations work.
        converted_time = datetime.datetime(converted_time.year, 
                                           converted_time.month, 
                                           converted_time.day,
                                           converted_time.hour,
                                           converted_time.minute,
                                           converted_time.second)    
    return converted_time

def iso_time_converter(time):
    """
    Convert time between datetime object and iso string format.

    Args:
    
    * time: datetime

    Returns:
        iso string

    """
    if type(time) == datetime.datetime:
        iso_time = time.isoformat()
    else:
        raise UserWarning('Given time must be datetime object.')
    return iso_time

def convert_time(time_coord):
    """
    Callback function to get iris cube time coordinate points, converted into 
    datetime format.
    
    Args:
    
    * time_coord: iris.coords
        Must be a time coordinate
    
    Returns:
        list of ISO time strings
    
    """
    time_unit = time_coord.units
    converted_time_points = []
    for time in time_coord.points:
        time = cube_time_converter(time, time_unit)
        converted_time_points.append(time)
    return converted_time_points

def sort_date(date, timestep):
    """
    Edit the date string according to the time step of the data.
    
    """
    if timestep in ['daily', 'weekly']:
        return iso_time_converter(date)[:10]
    elif timestep in ['monthly', 'seasonal', 'yearly']:
        return iso_time_converter(date)[:7]

def detect_timestep(cube):
    """
    Work out the time step of the data.
    
    """
    time_unit = cube.coord('time').units
    previous_date = None
    previous_diff = None
    timestep = None
    err_message = 'Inconsistent time steps in data, can not convert to valid '\
                  'CPT format.'
    for time_point in cube.coord('time').points:
        date = cube_time_converter(time_point, time_unit)
        if previous_date is not None:
            diff = (date - previous_date).days
            if diff == 1:
                if timestep == None:
                    timestep = 'daily'
                elif timestep != 'daily':
                    raise UserWarning(err_message)
                    
            elif diff == 7:
                if timestep == None:
                    timestep = 'weekly'
                elif timestep != 'weekly':
                    raise UserWarning(err_message)
                
            elif 28 <= diff <= 31:
                if timestep == None:
                    timestep = 'monthly'
                elif timestep != 'monthly':
                    raise UserWarning(err_message)
                
            elif 365 <= diff <= 366:
                if timestep == None:
                    timestep = 'yearly'
                elif timestep != 'yearly':
                    raise UserWarning(err_message)
                
            else:
                if timestep == None:
                    timestep = 'seasonal'
                else:
                    if (previous_diff - diff) > 4:
                        raise UserWarning(err_message)
                    elif timestep != 'seasonal':
                        raise UserWarning(err_message)
            previous_diff = diff
        previous_date = date
    return timestep
    
def transpose_data(xy_slice, x_coord, y_coord):
    """
    Transpose data so y is the first dimension and x is the second.
    
    Args:
    
    * xy_slice: iris.cube.Cube
        Must be a 2 dimensional iris cube.
    
    * x_coord: iris.coords.DimCoord
        Coordinate points in the x dimension.

    * y_coord: iris.coords.DimCoord
        Coordinate points in the y dimension.
    
    Returns
        iris.cube.Cube
    
    """
    assert xy_slice.ndim == 2, 'Can only transpose data with two dimensions.'
    
    x_index, = xy_slice.coord_dims(x_coord)
    y_index, = xy_slice.coord_dims(y_coord)
    xy_slice.transpose([y_index, x_index])
    return xy_slice

def look_for_coord(cube, coord_names, axis=None):
    """
    Look to find the given coordinate names in the cube. The first one found is
    returned so the coord_names list must be in order of preference. An 
    exception is raised if none are found.
    
    Args:
    
    * cube: iris.cube.Cube
    
    * coord_names: list
    
    * axis: string
    
    Returns:
        iris.coords instance
    
    """
    cube_coords = [coord.name() for coord in cube.coords()]
    found = False
    for coord_name in coord_names:
        if coord_name in cube_coords:
            coord = cube.coord(coord_name)
            found = True
            break
        
    if not found:
        if axis is None:
            axis = ''
        else:
            axis += ' '
        raise UserWarning('Could not find suitable %scoordinate.' % axis)
    else:
        return coord

def get_xy_coords(cube):
    """
    Get the coordinates describing the x and y dimensions of the data set.
    
    Args:
    
    * cube: iris.cube.Cube

    Returns:
        iris.coords instances
    
    """
    try:
        x_coord = cube.coord(axis='X')
    except iris.exceptions.CoordinateNotFoundError:
        # Try find an x coordinate manually.
        x_coord = look_for_coord(cube, ['longitude', 'grid_longitude'], 'x')
    
    try:
        y_coord = cube.coord(axis='Y')
    except iris.exceptions.CoordinateNotFoundError:
        # Try find the y coordinate manually.
        y_coord = look_for_coord(cube, ['latitude', 'grid_latitude'], 'y')
        
    return x_coord, y_coord

def make_coords_dimensions(cube, *coords):
    """
    Make sure all the given coordinates are dimension coordinates. I.e. if a 
    coordinate is a scalar (coordinate with length 1 but not a dimension) add a
    dimension to to the data.
    
    Args:
    
    * cube: iris.cube.Cube
    
    * coords: string or iris.coords instance
        Provide the coordinates to check as separate arguments.
    
    Returns:
        iris.cube.Cube
    
    """
    for coord in coords:
        if len(coord.points) == 1:
            cube = iris.util.new_axis(cube, scalar_coord=coord)
    return cube

def get_main_header_metadata(cubelist):
    """
    Get the main header metadata which must be printed at the top of every 
    file. 
    
    """
    header_dict = {}
    
    var_names = []
    for cube in cubelist:
        if cube.name() not in var_names:
            var_names.append(cube.name())
            
    header_dict['nfields'] = len(var_names)
    
    return header_dict

def get_slice_header_dict(xy_slice, x_coord, y_coord):
    """
    Get the specific header metadata which must be printed at the top of each 
    data chunk (xy slice of data).
    
    """
    header_dict = {}
    for coord in xy_slice.coords():
        # This metadata comes from coordinates which are not x or y, these are
        # coordinates with length one.
        if len(coord.points) == 1 and coord not in [x_coord, y_coord]:
            # Check it is a CPT recognised coordinate.
            coord_dict = cpt_convert_dict.get(coord.name())
            if coord_dict:
                if coord_dict.get('callback'):
                    coord_point, = coord_dict['callback'](coord)
                else:
                    coord_point, = coord.points
                header_dict[coord_dict['cpt_name']] = coord_point
    
    # Add the compulsory CPT tags.
    header_dict['nrow']  = len(xy_slice.coord(y_coord).points)
    header_dict['ncol']  = len(xy_slice.coord(x_coord).points)
    header_dict['row']   = cpt_convert_dict[y_coord.name()]['cpt_name']
    header_dict['col']   = cpt_convert_dict[x_coord.name()]['cpt_name']
    
    # Add optional CPT tags.
    if xy_slice.name() != 'unknown':
        cpt_field_name = cpt_field_name_dict.get(xy_slice.name())
        if cpt_field_name is None:
            cpt_field_name = xy_slice.name()    
        header_dict['field'] = cpt_field_name
    
    if xy_slice.units != 'unknown':
        header_dict['units'] = str(xy_slice.units)
        
    if hasattr(xy_slice.data, 'fill_value'):
        header_dict['missing'] = xy_slice.data.fill_value
        
    return header_dict

def write_main_header(outfile, main_header_dict):
    """
    Write CPT tags which belong at the top of the file.
    
    """
    outfile.write('xmlns:cpt=http://iri.columbia.edu/CPT/v10/\n') 
    for key, val in main_header_dict.items():
        outfile.write('cpt:{var}={val}\n'.format(var=key, val=val))

def write_slice_header(outfile, slice_header_dict):
    """
    Write CPT tags which belong at the top of each xy data slice.
    
    """
    header_line = []
    for key, val in slice_header_dict.items():
        header_line.append('cpt:{var}={val}'.format(var=key, val=val))
    header_line = ', '.join(header_line)
    outfile.write(header_line)
    outfile.write('\n')
    

def write_data(outfile, xy_slice, x_coord, y_coord, simple=False,
                timestep=None):
    """
    Write the data with xy coordinates wrapped around.
    
    """
    x_values = xy_slice.coord(x_coord).points
    y_values = xy_slice.coord(y_coord).points
    
    # x and y dimensions must be in the right order.
    xy_slice = transpose_data(xy_slice, x_coord, y_coord)
    # Place y coordinate values in a new first column within the data.
    data = numpy.insert(xy_slice.data, 0, y_values, axis=1)
    data = numpy.flipud(data)
    
    if simple:
        assert timestep is not None, 'Timestep has not been calculated.'
        # If simple just add the time data to the first column of the first 
        # row.
        time = convert_time(xy_slice.coord('time'))[0]
        time = sort_date(time, timestep)
        outfile.write('%s\t%s\n' % (time, 
                                    '\t'.join([str(x) for x in x_values])))
    else:
        # Write the x values above the data. To make sure each x value lines up 
        # with its correct column, add a tab space so the first column (where y 
        # values have just been added) is skipped.
        outfile.write('\t%s\n' % '\t'.join([str(x) for x in x_values]))
    # Write the data underneath.
    numpy.savetxt(outfile, data, delimiter='\t')

def cpt_converter(loadpaths, savepath, constraints=None, simple=False):
    """
    Main function. Load data, convert it and save.
    
    """
    check_filepath(savepath)
    
    # Load in data.
    data_cubelist = iris.load(loadpaths, constraints)
    
    if simple: 
        with open(savepath, 'a') as outfile:
            for cube in data_cubelist:
                timestep = detect_timestep(cube)
                x_coord, y_coord = get_xy_coords(cube)
                cube = make_coords_dimensions(cube, x_coord, y_coord)
                # Data must be cut into 2 dimensional xy slices, then each 
                # slice converted and added to the outfile.
                for xy_cube in cube.slices([x_coord.name(), y_coord.name()]):
                    write_data(outfile, xy_cube, x_coord, y_coord, simple=True,
                               timestep=timestep)
    else:
        # Gather metadata and write to file.    
        with open(savepath, 'a') as outfile:
            main_header_dict = get_main_header_metadata(data_cubelist)
            write_main_header(outfile, main_header_dict)
    
            for cube in data_cubelist:
                timestep = detect_timestep(cube)
                x_coord, y_coord = get_xy_coords(cube)
                cube = make_coords_dimensions(cube, x_coord, y_coord)
                # Data must be cut into 2 dimensional xy slices, then each 
                # slice converted and added to the outfile.
                for xy_cube in cube.slices([x_coord.name(), y_coord.name()]):
                    slice_header_dict = get_slice_header_dict(xy_cube, x_coord, 
                                                              y_coord)
                    slice_header_dict['T'] = sort_date(slice_header_dict['T'], 
                                                        timestep)
                    if slice_header_dict.get('S'):
                        slice_header_dict['S'] = sort_date(slice_header_dict['S'], 
                                                           'daily')
                    write_slice_header(outfile, slice_header_dict)
                    write_data(outfile, xy_cube, x_coord, y_coord)
    
    return savepath
    
cpt_convert_dict = {'time'                    : {'cpt_name' : 'T',
                                                 'callback' : convert_time},
                    'forecast_reference_time' : {'cpt_name' : 'S',
                                                 'callback' : convert_time},
                    'latitude'                : {'cpt_name' : 'Y'},
                    'longitude'               : {'cpt_name' : 'X'},
                    'grid_latitude'           : {'cpt_name' : 'Y'},
                    'grid_longitude'          : {'cpt_name' : 'X'},
                    'pressure'                : {'cpt_name' : 'P'},
                    'realization'             : {'cpt_name' : 'M'},
                    'model_level_number'      : {'cpt_name' : 'level'}}
