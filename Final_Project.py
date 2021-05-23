# Name: Dustin McGrew
# Created Date: 10/29/2018
# Last Modified Date: 10/12/2018
# Description: Final Project. Automated mapping and analysis of real-time USGS earthquake data.

import arcpy, numpy, requests, json, os, datetime

# Output shapefile parameter.
outFC = arcpy.GetParameterAsText(0)
# Input map document parameter.
inMxd = arcpy.GetParameterAsText(1)
# Input time period for Earthquake query.
input_time_period = arcpy.GetParameterAsText(2)
# Input magnitude for Earthquake query.
input_magnitude = arcpy.GetParameterAsText(3)
# Input cell height for quadrat grid.
input_cell_height = float(arcpy.GetParameterAsText(4))
# Input cell width for quadrat grid.
input_cell_width = float(arcpy.GetParameterAsText(5))
# Boolean variable for generating earthquake text file
generate_earthquake_text = arcpy.GetParameter(6)
# Boolean variable for generating global earthquake map.
generate_map = arcpy.GetParameter(7)
# Boolean variable for generating program log.
generate_log = arcpy.GetParameter(8)

# Set the arcpy environment workspace to the output shapefile directory.
arcpy.env.workspace = os.path.dirname(outFC)
# Set feature class variable path.
featClass = os.path.basename(outFC)

# Delete the feature class variable if it exists.
if arcpy.Exists(outFC):
    arcpy.Delete_management(outFC)

# Create the point feature class for earthquake locations.
arcpy.CreateFeatureclass_management(arcpy.env.workspace, featClass, "POINT")

# Create the attribute fields for the earthquake feature class.
arcpy.AddField_management(featClass, "USGS_ID", "TEXT")
arcpy.AddField_management(featClass, "Title", "TEXT")
arcpy.AddField_management(featClass, "Magnitude", "FLOAT")
arcpy.AddField_management(featClass, "Magn_Type", "TEXT")
arcpy.AddField_management(featClass, "Location", "TEXT")
arcpy.AddField_management(featClass, "Longitude", "FLOAT")
arcpy.AddField_management(featClass, "Latitude", "FLOAT")
arcpy.AddField_management(featClass, "Depth_km", "FLOAT")
arcpy.AddField_management(featClass, "Type", "TEXT")
arcpy.AddField_management(featClass, "Alert", "TEXT")
arcpy.AddField_management(featClass, "Tsunami", "TEXT")
arcpy.AddField_management(featClass, "Felt", "TEXT")
arcpy.AddField_management(featClass, "Link", "TEXT")

# Set variable for the URL linking to USGS earthquake data.
base_url = "https://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/"


# Create function for querying the specified time period from URL.
# Create additional function for querying the magnitude type from URL.
# Create the desired url based on query parameters.
def select_time_period(time, magnitude_type):

    new_url = "https://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/"

    def select_magnitude_type(magnitude_type):
        if magnitude_type == "Significant Earthquakes":
            new_url = base_url + "significant.geojson"
            return new_url
        elif magnitude_type == "M4.5+":
            new_url = base_url + "4.5.geojson"
            return new_url
        elif magnitude_type == "M2.5+":
            new_url = base_url + "2.5.geojson"
            return new_url
        elif magnitude_type == "M1.0+":
            new_url = base_url + "1.0.geojson"
            return new_url
        elif magnitude_type == "All Earthquakes":
            new_url = base_url + "all.geojson"
            return new_url

    if time == "Hour":
        update_url = select_magnitude_type(magnitude_type)
        if len(magnitude_type) == 23:
            update_url = update_url[:69] + '_hour' + update_url[69:]
            return update_url
        else:
            update_url = update_url[:61] + '_hour' + update_url[61:]
            return update_url

    elif time == "Day":
        update_url = select_magnitude_type(magnitude_type)
        if len(magnitude_type) == 23:
            update_url = update_url[:69] + '_day' + update_url[69:]
            return update_url
        else:
            update_url = update_url[:61] + '_day' + update_url[61:]
            return update_url

    elif time == "Week":
        update_url = select_magnitude_type(magnitude_type)
        if len(magnitude_type) == 23:
            update_url = update_url[:69] + '_week' + update_url[69:]
            return update_url
        else:
            update_url = update_url[:61] + '_week' + update_url[61:]
            return update_url

    elif time == "Month":
        update_url = select_magnitude_type(magnitude_type)
        if len(magnitude_type) == 23:
            update_url = update_url[:69] + '_month' + update_url[69:]
            return update_url
        else:
            update_url = update_url[:61] + '_month' + update_url[61:]
            return update_url


arcpy.AddMessage("Attempting to connect to the USGS server...")

try:
    # Attempt to connect to the USGS website and grab the json text.
    url = select_time_period(input_time_period, input_magnitude)
    req = requests.get(url)
    decoded = json.loads(req.text)
    # Create field list for storing the JSON text.
    fieldList = ["SHAPE@", "USGS_ID", "Title", "Magnitude", "Magn_Type",
                "Location", "Longitude", "Latitude", "Depth_km", "Type",
                "Alert", "Tsunami", "Felt", "Link"]
    arcpy.AddMessage("The connection was successful.")
except:
    print "Error connecting to the USGS website. Please check if the online server is available."
    quit()


# Create function for handling null values in JSON text.
def handle_null_values(input):
    if input is magnitude:
        if input is None:
            input = 0
            return input
        else:
            return input
    else:
        if input is None:
            input = "None"
            return input
        else:
            return input


arcpy.AddMessage("Attempting to add attributes to the feature class...")

# Open the feature class and its attribute table using insert cursor variable.
with arcpy.da.InsertCursor(featClass, fieldList) as isCursor:
    # Generate earthquake data text file if boolean parameter is true.
    if generate_earthquake_text is True:
        # Set the text file path, remove if it exists.
        text_path = os.path.dirname(outFC) + os.sep + "EarthquakeData.txt"
        if os.path.exists(text_path):
            os.remove(text_path)
        # Open and write the column titles to the text file.
        text_file = open(text_path, "w")
        text_file.write("USGS_ID, Title, Magnitude, Magn_Type, Location," \
                        "Longitude, Latitude, Depth_km, Type, Alert, Tsunami, Felt, Link")
        text_file.write("\n")

    # Loop through the features stored in the JSON text variable.
    for feature in decoded['features']:
        # Set geometry variables for earthquake points using JSON text query.
        geometry_data = feature['geometry']['coordinates']
        xCoord = geometry_data[0]
        yCoord = geometry_data[1]
        zDepth = geometry_data[2]
        # Set variables for attribute table using JSON text query.
        USGS_ID = feature['properties']["ids"]
        USGS_ID_Fix = USGS_ID[1:10]
        title = feature['properties']['title']
        magnitude = feature['properties']['mag']
        magnitude_type = feature['properties']['magType']
        location = feature['properties']['place']
        etype = feature['properties']['type']
        alert = feature['properties']['alert']
        tsunami = feature['properties']['tsunami']
        felt = feature['properties']['felt']
        url = feature['properties']['url']
        # Remove commas in certain variables for easy viewing in excel.
        title_fix = title.replace(",", "")
        location_fix = location.replace(",", "")
        # Pass the original attribute variables through the null function.
        new_title = handle_null_values(title)
        new_magnitude = handle_null_values(magnitude)
        new_magnitude_type = handle_null_values(magnitude_type)
        new_location = handle_null_values(location)
        new_etype = handle_null_values(etype)
        new_alert = handle_null_values(alert)
        new_tsunami = handle_null_values(tsunami)
        new_felt = handle_null_values(felt)
        new_url = handle_null_values(url)

        collect = USGS_ID_Fix, title, magnitude
        # Write to the earthquake text data file if boolean parameter is true.
        if generate_earthquake_text is True:
            text_file.write("%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s"
                            % (USGS_ID_Fix, title_fix, magnitude, magnitude_type,
                                location_fix, xCoord, yCoord, zDepth, etype, alert,
                                tsunami, felt, url))
            text_file.write("\n")
        # Collect the attribute data into a single variable.
        newPoint = [arcpy.Point(xCoord, yCoord), USGS_ID_Fix, new_title, new_magnitude,
                    new_magnitude_type, new_location, xCoord, yCoord, zDepth, new_etype,
                    new_alert, new_tsunami, new_felt, new_url]
        # Add earthquake point features to ArcMap with the specified attributes.
        isCursor.insertRow(newPoint)


if generate_earthquake_text is True:
    text_file.close()

arcpy.AddMessage("The feature class was successfully updated.")

# Delete the insert cursor variable from memory.
del isCursor

arcpy.AddMessage("Attempting to create quadrat grid...")

# Set input variables for generating quadrat grid.
grid_path = arcpy.env.workspace + os.sep + "grid.shp"
origin_coordinate = "-180.0001388888889 -85.22193775799991"
y_axis_coord = "-179.9999999999999 -75.22193775799991"
opposite_corner = '179.9998611111108 83.999861111111'
cellSizeHeight = str(input_cell_height)
cellSizeWidth = str(input_cell_width)
numRows = ''
numColumns = ''
labels = 'NO_LABELS'
templateExtent = '#'
geometryType = 'POLYGON'
# Run the arcpy fishnet_management function to generate quadrat grid.
arcpy.CreateFishnet_management(grid_path, origin_coordinate, y_axis_coord, cellSizeWidth,
                                    cellSizeHeight, numRows, numColumns, opposite_corner, labels,
                                    templateExtent, geometryType)

arcpy.AddMessage("The quadrat grid was successfully created.")

arcpy.AddMessage("Attempting to calculate grid statistics...")

# Credit goes to Dr. Patricia Drews for this method for calculating VMR in ArcMap.
# Set variable paths for spatial join function.
join_path = arcpy.env.workspace + os.sep + "join.shp"
join_table = arcpy.env.workspace + os.sep + "join.dbf"
temp_table = arcpy.env.workspace + os.sep + "join_stat.txt"

# Perform spatial join on quadrat grid and earthquake points.
arcpy.SpatialJoin_analysis(grid_path, outFC, join_path)

# Calculate the statistics for the joined shapefile.
# This is necessary to get the standard deviation and variance of frequency count.
arcpy.Statistics_analysis(join_table, temp_table, [["Join_Count", "STD"]])
# Set variable for the joined statistics shapefile.
open_stat_file = open(arcpy.env.workspace+os.sep+"join_stat.txt", "r")
stat_file = open_stat_file.readlines()
stat_file_list = []
open_stat_file.close()
# Loop through the statistics file and append the data to list.
for x in stat_file:
    stat_file_list.append(x)
# Split the file data.
split_index = stat_file_list[1].split(",")

# Set variables for the number of earthquake points and individual grids.
# This is necessary to calculate the mean for the frequency count.
point_count = str(arcpy.GetCount_management(outFC))
grid_count = str(arcpy.GetCount_management(arcpy.env.workspace+os.sep+"grid.shp"))

# Calculate the variance, mean, and variance-to-mean ratio.
# for earthquake points in each grid.
arcpy.AddMessage("Calculating VAR.")
VAR = float(split_index[2])
arcpy.AddMessage("Calculating MEAN.")
MEAN = float(point_count) / float(grid_count)
arcpy.AddMessage("Calculating VMR.")

if MEAN == 0:
    arcpy.AddMessage("Error! No earthquakes were found.")
    arcpy.Delete_management(featClass)
    arcpy.Delete_management(grid_path)
    arcpy.Delete_management(join_path)
    arcpy.Delete_management(temp_table)
    arcpy.Delete_management(arcpy.env.workspace + os.sep + "schema.ini")
    quit()
else:
    VMR = VAR / MEAN

# Function for determining the data point pattern based on VMR.
# VMR > 1.0 = Clustered.
# VMR = 1.0 = Random.
# VMR < 1.0 = Dispersed.
def pattern_func(vmr):
    if vmr > 1.0:
        return "Clustered"
    elif vmr == 1.0:
        return "Random"
    elif vmr < 1.0:
        return "Dispersed"


if generate_map is True:
    # Generate a global earthquake map if the boolean parameter is true.
    arcpy.AddMessage("Attempting to generate earthquake map...")

    # Create feature layers for the earthquake points and quadrat grid.
    arcpy.MakeFeatureLayer_management(featClass, "Earthquake_points")
    arcpy.MakeFeatureLayer_management(grid_path, "Quadrat_grid")

    # Set variable for the ArcMap document path.
    curMxd = arcpy.mapping.MapDocument(inMxd)
    # Set variable for the data frame in the map document.
    data_frame = arcpy.mapping.ListDataFrames(curMxd)[0]
    # Set variable for the legend in the data frame.
    legend = arcpy.mapping.ListLayoutElements(curMxd, "LEGEND_ELEMENT", "Legend")[0]
    # Create arcmap.mapping layers for each input layer in the map.
    earthquake_layer = arcpy.mapping.Layer("Earthquake_points")
    quadrat_layer = arcpy.mapping.Layer("Quadrat_grid")
    hillshade_layer = arcpy.mapping.Layer(arcpy.env.workspace+os.sep+"HillshadeRaster.lyr")
    elevation_layer = arcpy.mapping.Layer(arcpy.env.workspace+os.sep+"ElevationRaster.lyr")
    ocean_layer = arcpy.mapping.Layer(arcpy.env.workspace+os.sep+"Ocean.lyr")
    plate_layer = arcpy.mapping.Layer(arcpy.env.workspace+os.sep+"PlateBoundaries.lyr")
    # Import the layer symbology for the earthquake and quadrat layers.
    arcpy.ApplySymbologyFromLayer_management(earthquake_layer, arcpy.env.workspace+os.sep+"PointSymbology.lyr")
    arcpy.ApplySymbologyFromLayer_management(quadrat_layer, arcpy.env.workspace+os.sep+"GridSymbology.lyr")
    # Disable automatically adding features to the legend.
    legend.autoAdd = False
    # Add the hillshade, elevation, and ocean layers to the map.
    arcpy.mapping.AddLayer(data_frame, hillshade_layer, "BOTTOM")
    arcpy.mapping.AddLayer(data_frame, elevation_layer, "TOP")
    arcpy.mapping.AddLayer(data_frame, ocean_layer, "TOP")
    # Enable automatically adding features to the legend.
    legend.autoAdd = True
    # Add the tectonic plate, earthquake, and quadrat layers to the map.
    arcpy.mapping.AddLayer(data_frame, plate_layer, "TOP")
    arcpy.mapping.AddLayer(data_frame, earthquake_layer, "TOP")
    arcpy.mapping.AddLayer(data_frame, quadrat_layer, "TOP")
    # Set the data frame view.
    newExtent = data_frame.extent
    newExtent.XMin, newExtent.YMin = -180.0, -80.0
    newExtent.XMax, newExtent.YMax = 180.0, 90.0
    data_frame.extent = newExtent
    # Set the data frame scale.
    data_frame.scale = 105000000
    arcpy.RefreshActiveView()
    arcpy.RefreshTOC()
    # Set variable for recording the current date.
    date = datetime.datetime.now()
    cur_date = str(date.month)+"/"+str(date.day)+"/"+str(date.year)
    # Set variable for the number of earthquake points.
    eq_count = arcpy.GetCount_management(featClass)


    # Set function for calculating the mean value of earthquake properties.
    def calculate_mean_value(table, field):
        na = arcpy.da.TableToNumPyArray(table, field)
        return numpy.mean(na[field])


    # Set function for calculating the max value of earthquake properties.
    def calculate_max_value(table, field):
        na = arcpy.da.TableToNumPyArray(table, field)
        return numpy.max(na[field])


    # Set function for calculating the min value of earthquake properties.
    def calculate_min_value(table, field):
        na = arcpy.da.TableToNumPyArray(table, field)
        return numpy.min(na[field])


    # Set variables for the average, max, and min magnitude.
    avg_magnitude = float(calculate_mean_value(featClass, "Magnitude"))
    max_magnitude = float(calculate_max_value(featClass, "Magnitude"))
    min_magnitude = float(calculate_min_value(featClass, "Magnitude"))

    # Set variables for the average, max, and min depth.
    avg_depth = float(calculate_mean_value(featClass, "Depth_km"))
    max_depth = float(calculate_max_value(featClass, "Depth_km"))
    min_depth = float(calculate_min_value(featClass, "Depth_km"))

    # Set variable for the map text.
    text_element = arcpy.mapping.ListLayoutElements(curMxd, "TEXT_ELEMENT")[0]

    # Add the following variables dynamically to the map text.
    text_element.text = "Current Date: %s" % cur_date + "\n" + \
                        "Selected Time Period: %s" % input_time_period + "\n" + \
                        "Selected Magnitude: %s" % input_magnitude + "\n" + \
                        "Quadrat Cell Size: %f, %f" % (input_cell_height, input_cell_width) + "\n" + \
                        "VMR: %f" % VMR + "\n" + \
                        "Pattern: %s" % pattern_func(VMR) + "\n" + \
                        "Number of Earthquakes: %s" % eq_count + "\n" + \
                        "Average Magnitude: %f" % avg_magnitude + "\n" + \
                        "Maximum Magnitude: %f" % max_magnitude + "\n" + \
                        "Minimum Magnitude: %f" % min_magnitude + "\n" + \
                        "Average Depth: %f" % avg_depth + "\n" + \
                        "Maximum Depth: %f" % max_depth + "\n" + \
                        "Minimum Depth: %f" % min_depth

    # Export the pdf map to the workspace environment.
    arcpy.mapping.ExportToPDF(curMxd, arcpy.env.workspace+os.sep+"Earthquake.pdf")

    arcpy.AddMessage("The earthquake map was successfully generated.")

if generate_log is True:
    # A program log is generated if the boolean is true.
    # The log contains information for the earthquakes gathered on that date.
    arcpy.AddMessage("Generating program log...")
    cur_date = str(date.month) + "/" + str(date.day) + "/" + str(date.year)
    eq_count = arcpy.GetCount_management(featClass)
    avg_magnitude = float(calculate_mean_value(featClass, "Magnitude"))
    max_magnitude = float(calculate_max_value(featClass, "Magnitude"))
    min_magnitude = float(calculate_min_value(featClass, "Magnitude"))
    avg_depth = float(calculate_mean_value(featClass, "Depth_km"))
    max_depth = float(calculate_max_value(featClass, "Depth_km"))
    min_depth = float(calculate_min_value(featClass, "Depth_km"))

    log_path = arcpy.env.workspace + os.sep + "Log.txt"

    text_file = open(log_path, "w")
    text_file.write("Current Date: %s" % cur_date + "\n" + \
                    "Selected Time Period: %s" % input_time_period + "\n" + \
                    "Selected Magnitude: %s" % input_magnitude + "\n" + \
                    "Quadrat Cell Size: %f, %f" % (input_cell_height, input_cell_width) + "\n" + \
                    "VMR: %f" % VMR + "\n" + \
                    "Pattern: %s" % pattern_func(VMR) + "\n" + \
                    "Number of Earthquakes: %s" % eq_count + "\n" + \
                    "Average Magnitude: %f" % avg_magnitude + "\n" + \
                    "Maximum Magnitude: %f" % max_magnitude + "\n" + \
                    "Minimum Magnitude: %f" % min_magnitude + "\n" + \
                    "Average Depth: %f" % avg_depth + "\n" + \
                    "Maximum Depth: %f" % max_depth + "\n" + \
                    "Minimum Depth: %f" % min_depth)
    text_file.close()
arcpy.AddMessage("Deleting unnecessary files...")

# File Cleanup.
arcpy.Delete_management(grid_path)
arcpy.Delete_management(join_path)
arcpy.Delete_management(temp_table)
arcpy.Delete_management(arcpy.env.workspace+os.sep+"schema.ini")
