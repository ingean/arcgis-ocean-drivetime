import os
import arcpy

arcpy.env.overwriteOutput = True

SOURCE_GDB = r'D:\Data\Redningsselskapet\GeoAnalyse\RS_GeoAnalyse.gdb'
TARGET_GDB = r'D:\Data\Redningsselskapet\GeoAnalyse\RS_Analyseresultater.gdb'

MAX_DISTANCE = 46300
STATION_COUNT = 6
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
  with arcpy.da.SearchCursor(fc, ['OID@']) as cursor:
    for row in cursor:
      print(f'Removing overlap for feature: {row[0]}')

      if row[0] > 1:
        result = remove_overlap(row[0], row[0] - 1)

      if row[0] < STATION_COUNT:
        result = remove_overlap(row[0], row[0] + 1)

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

def main():
  arcpy.CheckOutExtension("Spatial")
  create_drivetime_area_pr_feature(STATIONS)
  remove_drivetime_overlap(STATIONS)

###############################################################################
if __name__ == "__main__":
  main()

