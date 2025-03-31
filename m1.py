import netCDF4 as nc
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from datetime import datetime

def read_era5_data(wind_filepath, wave_filepath, time_range, lat_range, lon_range):
    """
    Reads ERA5 wind and wave data for a specified time range and geographical range from separate netCDF files.
    Wave data is regridded to match the resolution of the wind data (0.25°).

    Parameters:
        wind_filepath (str): Path to the ERA5 wind netCDF file.
        wave_filepath (str): Path to the ERA5 wave netCDF file.
        time_range (tuple): Tuple of start and end times in ISO 8601 format (e.g., ('2024-01-01T00:00:00', '2024-01-02T00:00:00')).
        lat_range (tuple): Latitude range as (min_lat, max_lat).
        lon_range (tuple): Longitude range as (min_lon, max_lon).

    Returns:
        dict: A dictionary containing the extracted and regridded data arrays for specific variables (e.g., u10, v10, wave_height).
    """
    def extract_data(filepath, time_range, lat_range, lon_range, target_variables):
        # Open the netCDF file
        dataset = nc.Dataset(filepath, 'r')

        # Extract variables
        time_var = dataset.variables['valid_time']
        lat_var = dataset.variables['latitude']
        lon_var = dataset.variables['longitude']

        # Convert time_range strings to datetime objects
        start_time = datetime.fromisoformat(time_range[0])
        end_time = datetime.fromisoformat(time_range[1])

        # Convert time_range to indices
        time_units = time_var.units
        time_calendar = getattr(time_var, 'calendar', 'standard')
        start_time_index = nc.date2index(
            start_time, time_var, calendar=time_calendar, select='nearest'
        )
        end_time_index = nc.date2index(
            end_time, time_var, calendar=time_calendar, select='nearest'
        )

        # Find latitude and longitude indices within the specified range
        lat_indices = np.where((lat_var[:] >= lat_range[0]) & (lat_var[:] <= lat_range[1]))[0]
        lon_indices = np.where((lon_var[:] >= lon_range[0]) & (lon_var[:] <= lon_range[1]))[0]

        # Extract data for the specified time and region
        data = {}
        for var_name in target_variables:
            if var_name in dataset.variables:
                var = dataset.variables[var_name]
                data[var_name] = var[start_time_index:end_time_index+1, lat_indices.min():lat_indices.max()+1, lon_indices.min():lon_indices.max()+1]

        # Extract latitude and longitude values for regridding
        lat_values = lat_var[lat_indices]
        lon_values = lon_var[lon_indices]

        # Close the dataset
        dataset.close()
        return data, lat_values, lon_values

    def regrid_wave_data(wave_data, wave_lat, wave_lon, wind_lat, wind_lon):
        """
        Regrids wave data to match the resolution of the wind data.
        """
        grid_lat, grid_lon = np.meshgrid(wind_lat, wind_lon, indexing='ij')
        regridded_data = {}

        for var_name, values in wave_data.items():
            # Flatten the input wave data
            points = np.array([(lat, lon) for lat in wave_lat for lon in wave_lon])
            values_flat = values.reshape(-1, values.shape[-2] * values.shape[-1]).mean(axis=0)  # Time average

            # Perform grid interpolation
            regridded_data[var_name] = griddata(
                points, values_flat, (grid_lat, grid_lon), method='linear'
            )

        return regridded_data

    # Define target variables for wind and wave datasets
    wind_variables = ['u10', 'v10']  # Adjust based on the wind dataset
    wave_variables = ['shww']  # Example wave variable (significant wave height)

    # Extract wind data
    wind_data, wind_lat, wind_lon = extract_data(wind_filepath, time_range, lat_range, lon_range, wind_variables)

    # Extract wave data
    wave_data, wave_lat, wave_lon = extract_data(wave_filepath, time_range, lat_range, lon_range, wave_variables)

    # Regrid wave data to match wind data resolution
    regridded_wave_data = regrid_wave_data(wave_data, wave_lat, wave_lon, wind_lat, wind_lon)

    # Combine wind and regridded wave data into a single dictionary
    combined_data = {**wind_data, **regridded_wave_data}

    return combined_data, wind_lat, wind_lon

def visualize_data(data, wind_lat, wind_lon):
    """
    Visualizes wind and wave data using matplotlib.

    Parameters:
        data (dict): Dictionary containing wind and wave data arrays.
        wind_lat (array): Latitude values for wind data.
        wind_lon (array): Longitude values for wind data.
    """
    u10 = data.get('u10')
    v10 = data.get('v10')
    swh = data.get('shww')

    # Plot u10 and v10 (wind vectors)
    plt.figure(figsize=(12, 6))
    plt.quiver(wind_lon, wind_lat, u10.mean(axis=0), v10.mean(axis=0), scale=50, color='blue', label='Wind Vectors')
    plt.title('Average Wind Field (u10, v10)')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.legend()
    plt.grid()
    plt.show()

    # Plot significant wave height (swh)
    if swh is not None:
        plt.figure(figsize=(12, 6))
        plt.contourf(wind_lon, wind_lat, swh.mean(axis=0), levels=20, cmap='viridis')
        plt.colorbar(label='Significant Wave Height (m)')
        plt.title('Average Significant Wave Height')
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        plt.grid()
        plt.show()


# Example usage
wind_filepath = r"D:\data\wind_2013.nc"
wave_filepath = r"D:\data\wave_2013.nc"
time_range = ("2013-06-13T09:00:00", "2013-06-18T09:00:00")
lat_range = (30, 40)  # Example latitude range
lon_range = (120, 130)  # Example longitude range

data, wind_lat, wind_lon = read_era5_data(wind_filepath, wave_filepath, time_range, lat_range, lon_range)
visualize_data(data, wind_lat, wind_lon)

