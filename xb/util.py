import numpy as np
import pandas as pd
import tifffile
import shapely
from rasterio import features
from rasterio import Affine
import xml.etree.ElementTree as ET

def get_image_shape(image_path):
    """ Get image shape
    """
    
    with tifffile.TiffFile(image_path) as tif:
        img_shape = tif.pages[0].shape
        
    return img_shape

def get_ome_schema(image_path):
    with tifffile.TiffFile(image_path) as tif:
        ome_xml = tif.ome_metadata
    ome_xml_root = ET.fromstring(ome_xml)

    # Extract the namespace associated with the OME tag
    ome_namespace = ome_xml_root.tag.split('}')[0].strip('{')

    return ome_namespace
    
def extract_physical_sizes(image_path):
    """ Extract physical sizes from OME-XML assuming schema version 2016-06.
    """
    
    schema2016 = "http://www.openmicroscopy.org/Schemas/OME/2016-06"
    
    assert str(image_path).endswith('.ome.tif'), 'Image must be an OME-TIFF.'
    schema = get_ome_schema(image_path)
    assert schema == schema2016, f"Unexpected schema: {schema}"
    
    with tifffile.TiffFile(image_path) as tif:
        ome_xml = tif.ome_metadata

    ome_xml_root = ET.fromstring(ome_xml)

    image_data = {}

    for image in ome_xml_root.findall('.//ns0:Image', namespaces={'ns0': schema2016}):
        image_name = image.attrib.get('Name', None)

        pixels_element = image.find('ns0:Pixels', namespaces={'ns0': schema2016})

        if pixels_element is not None:
            physical_size_x = pixels_element.attrib.get('PhysicalSizeX', None)
            physical_size_y = pixels_element.attrib.get('PhysicalSizeY', None)
            physical_size_z = pixels_element.attrib.get('PhysicalSizeZ', None)

            image_data[image_name] = {
                'PhysicalSizeX': physical_size_x,
                'PhysicalSizeY': physical_size_y,
                'PhysicalSizeZ': physical_size_z
            }

    return image_data


def convert_polygons_to_label_image_xenium(
    df: pd.DataFrame, 
    img_shape: tuple, 
    x_col:str = "vertex_x", 
    y_col:str = "vertex_y", 
    label_col:str = "label_id",
    verbose: bool = False,
) -> np.array:
    ''' Create label image from a dataframe with polygons 
    
    Xenium files to load as `df`: cell_boundaries.parquet or nucleus_boundaries.parquet (pd.read_parquet(path_to_file)).
    Note that polygon coordinates need to be transformed into pixel coordinates of the according image 
    (morphology.ome.tif).
    
    Arguments
    ---------
    df: pd.Dataframe
        Dataframe containing polygon coordinates. Columns: "cell_id", "vertex_x", "vertex_y", "label_id"
    img_shape: tuple
        Shape of the image the polygons are drawn on.
    x_col: str
        Column name of the polygon vertices' x-coordinates in the dataframe.
    y_col: str
        Column name of the polygon vertices' y-coordinates in the dataframe.
    label_col: str
        Column name of the polygon/cell label in the dataframe.
    verbose: bool
        If True, print warnings for invalid polygons.

    Returns:
    ----------
    np.array 
        Label image with the same shape as the input image.
    '''    
    
    # Initialize label image
    labels = df[label_col].unique()
    max_label = np.max(labels)
    dtype = np.uint32 if max_label < np.iinfo(np.uint32).max else np.uint64
    
    assert max_label < np.iinfo(dtype).max, f"Label values exceed {dtype} range ({max_label})."
    
    label_image = np.zeros(img_shape, dtype=dtype)
    
    # Assert that min and max x and y are within the image shape
    x_min, x_max, y_min, y_max = df[x_col].min(), df[x_col].max(), df[y_col].min(), df[y_col].max()
    assert x_min >= 0 and x_max < img_shape[1], f"Polygon X coords ({x_min}, {x_max}) exceed image shape {img_shape}."
    assert y_min >= 0 and y_max < img_shape[0], f"Polygon Y coords ({y_min}, {y_max}) exceed image shape {img_shape}."
        
    # Iterate over each label id and map the corresponding polygon to the label image
    label_grouped_dfs = df.groupby(label_col)[[x_col, y_col]]
    
    for label_id, df_ in label_grouped_dfs:
        
        # Skip polygons with less than 3 vertices
        if len(df_) < 3:
            if verbose:
                print(f"Skipping invalid Polygon at cell_id = {label_id}")
            continue
        
        # Get polygon and crop dimensions
        polygon = shapely.geometry.Polygon(df_[[x_col, y_col]].values)
        
        minx, miny, maxx, maxy = polygon.bounds
        minx, miny, maxx, maxy = int(minx), int(miny), int(maxx), int(maxy)
        
        # Skip polygons with zero width or height
        if (int(maxx - minx)==0) or (int(maxy-miny)==0): 
            if verbose:
                print(f"Skipping invalid Polygon at cell_id = {label_id}")
            continue
        
        # Rasterize polygon on little crop of the image
        cell_image_crop = features.rasterize(
            [(polygon, label_id)],
            out_shape=(maxy-miny, maxx-minx),
            transform = Affine.translation(minx, miny),
            fill=0,
            dtype=dtype
        )
        
        # Update label image
        label_image[miny:maxy, minx:maxx] = np.where(
            label_image[miny:maxy, minx:maxx] == 0, cell_image_crop, label_image[miny:maxy, minx:maxx]
        )
       
    return label_image 