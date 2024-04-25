import os
import arcpy

arcpy.env.overwriteOutput = True

SOURCE_GDB = r'D:\Data\Redningsselskapet\GeoAnalyse\RS_GeoAnalyse.gdb'
TARGET_GDB = r'D:\Data\Redningsselskapet\GeoAnalyse\RS_Analyseresultater.gdb'

MAX_DISTANCE = 46300
STATIONS = os.path.join(SOURCE_GDB, 'Stasjoner_vestland')
COST_RASTER = os.path.join(SOURCE_GDB, 'KOST_Distance_Shore_Speed')

def create_drivetime_area_pr_feature(fc):
    
  with arcpy.da.SearchCursor(fc, ['OID@', 'SHAPE@XY', 'Navn']) as cursor:
    for row in cursor:   
      print(f'Creating drivetime area for feature: {row[2]}...')

      
      station_pt = arcpy.Point(row[1][0], row[1][1])
      station_geometry = arcpy.PointGeometry(station_pt, spatial_reference=25833) 
     
      result = arcpy.sa.DistanceAccumulation(station_geometry, in_cost_raster=COST_RASTER)
      result.save(os.path.join(TARGET_GDB, f'drivetime_{row[0]}'))
      
     

def remove_drivetime_overlap(fc):
  station_count = count(STATIONS)
  
  with arcpy.da.SearchCursor(fc, ['OID@']) as cursor:
    for row in cursor:
      print(f'Removing overlap for feature: {row[0]}')

      for i in range(station_count):
        id = i + 1 
        if row[0] != id:
          result = remove_overlap(row[0], id)

      result = arcpy.sa.Con(result, result, "", f"VALUE <= {MAX_DISTANCE}")
      result.save(os.path.join(TARGET_GDB, f'drivetime_{row[0]}_{MAX_DISTANCE}'))


def remove_overlap(idx1, idx2):  
  raster1 = os.path.join(TARGET_GDB, f'drivetime_{idx1}')
  raster2 = os.path.join(TARGET_GDB, f'drivetime_{idx2}')    
  overlap = arcpy.sa.Raster(raster1) < arcpy.sa.Raster(raster2)

  if idx1 < idx2:
    expr = "VALUE = 1"
  else: 
    expr = "VALUE IS NULL"

  result = arcpy.sa.Con(overlap, raster1, None, expr)    
  result.save(os.path.join(TARGET_GDB, f"drivetime_{idx1}"))  
  return result

def convert_to_polygons():
  station_count = count(STATIONS)
  
  for i in range(station_count):
    id = i + 1
    print(f"Converting raster: drivetime_{id}_{MAX_DISTANCE}")
    
    raster = os.path.join(TARGET_GDB, f'drivetime_{id}_{MAX_DISTANCE}')
    fc = os.path.join(TARGET_GDB, f'drivetime_{id}_poly')

    reclass_raster = reclass(raster)
    raster_to_poly(reclass_raster, fc)

def reclass(raster):
  return arcpy.sa.Reclassify(
      in_raster=raster,
      reclass_field="Value",
      remap="0 46300 1",
      missing_values="DATA" 
    )

def raster_to_poly(raster, fc):  
  arcpy.conversion.RasterToPolygon(
      in_raster=raster,
      out_polygon_features=fc,
      simplify="SIMPLIFY",
      raster_field="Value",
      create_multipart_features="MULTIPLE_OUTER_PART",
      max_vertices_per_feature=None
    )

def count(dataset):
  result = arcpy.management.GetCount(dataset)
  print(f'{dataset} has {result[0]} records')
  return int(result[0])

def main():
  arcpy.CheckOutExtension("Spatial")
  #create_drivetime_area_pr_feature(STATIONS)
  #remove_drivetime_overlap(STATIONS)
  convert_to_polygons()

###############################################################################
if __name__ == "__main__":
  main()

