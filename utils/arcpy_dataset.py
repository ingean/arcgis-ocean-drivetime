import os
import arcpy
from .log import Log
from .filepaths import Path
log = Log(arcgis=True)

def create_featureclass(fc: str, geometry_type: str = 'POINT', wkid: int = 25833, template_fc: str = None, delete_existing: bool = False) -> str:
  """Creates a new feature class (will delete first if already exist)
  Args:
    fc (str): Full path to the feature class
    geometry_type (str): POINT, MULTIPOINT, POLYGON, POLYLINE, or MULTIPATCH  
    wkid (int): Well know ID of the spatial reference (e.g. 25833)
    template_fc (str): Optional; Full path to template feature class
    delete_existing (bool): Deletes existing feature class, default is False
  """
  if delete_existing:
    delete_featureclass(fc)

  log.info(f"Creating feature class: {fc}...")

  if (template_fc is not None):
    spatial_ref = arcpy.Describe(template_fc).spatialReference
  else:
    spatial_ref = arcpy.SpatialReference(wkid)
 
  path = os.path.split(fc)
  return arcpy.CreateFeatureclass_management(path[0], path[1], geometry_type, template_fc, 
                                    "DISABLED", "ENABLED", spatial_ref)

def delete_featureclass(dataset: str) -> None:
  """Deletes the feature class if it exists
  
  Args:
    dataset (str): Full path to the feature class or table
  """
  if arcpy.Exists(dataset):
    log.info("Deleting existing featureclass" + dataset + "...")
    arcpy.Delete_management(dataset)

def featureclass_to_dict(dataset: str, fields: list = "*") -> dict:
  """Returns a dictionary of all the records with the specified fields
  
  Args:
    dataset (str): Full path to the table or feature class
    fields (list): List of strings with the names of the fields to include in the dictionary. Default is all fields

  Returns:
    dict: A dictionary including all the records and specified fields of the feature class
  """

  return {row[0]:row[1:] for row in arcpy.da.SearchCursor(dataset, fields)} 

def featureclass_to_list(dataset: str, fields: list = "*", sort_field_index: int = -1) -> list:
  """Returns a list of all the records with the specified fields
  
  Args:
    dataset (str): Full path to the table or feature class
    fields (list): Optional; List of strings with the names of the fields to include in the dictionary. Default is all fields
    sort_field_index (int): Optional; Index of field to sort the list. Default is not to sort the list.

  Returns:
    dict: A list including all the records and specified fields of the feature class
  """
  data = [row for row in arcpy.da.SearchCursor(dataset, fields)]

  if sort_field_index > 0: 
   return data.sort(key=lambda tup: tup[sort_field_index])
  else:
    return data  

def add_fields(dataset: str, field_defs: list) -> None:
  log.info(f'Adding {len(field_defs)} fields to input fc {dataset}...')
  """Adds a new fields to a dataset if not already existing."""

  for field_def in field_defs:
    field_type = field_def[1] if len(field_def) == 2 else None
    add_field(dataset, field_def[0], field_type)

def add_field(dataset: str, field_name: str, field_type: str = "TEXT") -> int:
  """Adds a new field to a dataset if not already existing.

  Args:
    dataset (str): Full path to the table or feature class
    field_name (str): Name of new field.
    field_type (str): TEXT, FLOAT, DOUBLE, SHORT, LONG, DATE, BLOB, RASTER or GUID

  Returns:
    int: The index of the new field in the table or feature class
  """
  if not field_exists(dataset, field_name):
    log.info(f'Adding field {field_name} of type {field_type} to {dataset}')
    arcpy.management.AddField(dataset, field_name, field_type)
  else:
    log.info(f'Field {field_name} already exists in {dataset}')

  
  return find_field_index(dataset, field_name)

def field_exists(dataset: str, field_name: str) -> bool:
  """Return boolean indicating if field exists in the specified dataset."""
  return field_name in [field.name for field in arcpy.ListFields(dataset)]

def find_field_index(dataset: str, field_name: str) -> int:
  """Returns the index of a field in a feature class or table."""
  return [field.name for field in arcpy.ListFields(dataset)].index(field_name)

def export_dataset(input_dataset: str, output_dataset: str, query: str = None) -> arcpy.Result:
  """Copy the features from the input dataset meeting the query expression to the output dataset"""
  
  outpath = Path(output_dataset)
  
  return arcpy.conversion.FeatureClassToFeatureClass(
    input_dataset, outpath.path, 
    outpath.filename, query
  )

def count_rows(input_dataset: str) -> int:
  result = arcpy.GetCount_management(input_dataset)
  return int(result[0])