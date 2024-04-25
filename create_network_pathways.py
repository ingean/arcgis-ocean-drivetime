import os
import numpy as np
import arcpy
from utils.utilities import unique_id

arcpy.env.overwriteOutput = True

SOURCE_GDB = r'D:\Data\Redningsselskapet\GeoAnalyse\RS_GeoAnalyse.gdb'
TARGET_GDB = r'D:\Data\Redningsselskapet\GeoAnalyse\RS_Nettverk.gdb'

AOI = os.path.join(SOURCE_GDB, 'AOI_50nm_hav')
LATTICE_LINES = os.path.join(TARGET_GDB, 'lattice_lines')

def create_square_diagonals(clipped_lattice):
  print(f"Creating diagonal lines for {clipped_lattice}...")        
  sr = arcpy.Describe(AOI).spatialReference
  
  internal_lines = None
  line_coords_array = []

  with arcpy.da.SearchCursor(clipped_lattice, ['SHAPE@']) as rows:
    for i, row in enumerate(rows):
      geom = row[0]
      center = geom.centroid
      for part in geom:
        line_coords_array.append([part[0].X, part[0].Y, part[2].X, part[2].Y])
        line_coords_array.append([part[1].X, part[1].Y, part[3].X, part[3].Y])
        #for pnt in part:
          #if pnt:
            #line_coords_array.append([center.X, center.Y, pnt.X, pnt.Y])
        del part
      del geom
      del center

  np_array = np.array(line_coords_array)
  struct_array = np.core.records.fromarrays(np_array.transpose(), np.dtype([('StartX', 'f8'), ('StartY', 'f8'),
                                                                            ('EndX', 'f8'), ('EndY', 'f8')]))
  xy_table = os.path.join(TARGET_GDB, f'xy_table_{unique_id()}')
  arcpy.da.NumPyArrayToTable(struct_array, xy_table)
  
  diagonal_lines = os.path.join(TARGET_GDB, f'diagonal_lines_{unique_id()}')
  if arcpy.Exists(xy_table):
    internal_lines = arcpy.management.XYToLine(xy_table, diagonal_lines,
                                              "StartX", "StartY", "EndX", "EndY", "GEODESIC", None, sr).getOutput(0)
  else:
    print("Failed to create diagonal lines!")

  del line_coords_array
  del np_array
  del struct_array

  return internal_lines
       
def get_extent(fc):
  desc = arcpy.Describe(fc)
  return desc.extent

def create_lattice(template_fc):
  print(f"Creating lattice for template: {template_fc}")
  ext = get_extent(template_fc)
  return arcpy.management.CreateFishnet(
          out_feature_class=r"memory\fishnet",
          origin_coord=f"{ext.XMin} {ext.YMin}",
          y_axis_coord=f"{ext.XMin} {ext.YMin + 10}",
          #origin_coord="-158646.971714386 6701394.84091875",
          #y_axis_coord="-158646.971714386 6701404.84091875",
          cell_width=1000,
          cell_height=1000,
          number_rows=None,
          number_columns=None,
          corner_coord=f"{ext.XMax} {ext.YMax}",
          #corner_coord="133897.261703125 7083731.82920096",
          labels="NO_LABELS",
          template='-158646.971714386 6701394.84091875 133897.261703125 7083731.82920096 PROJCS["ETRS_1989_UTM_Zone_33N",GEOGCS["GCS_ETRS_1989",DATUM["D_ETRS_1989",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",15.0],PARAMETER["Scale_Factor",0.9996],PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]',
          geometry_type="POLYGON"
        ).getOutput(0)

def clip_lattice(lattice, clip_fc):
  print(f"Clipping lattice for area: {clip_fc}")  
  return arcpy.analysis.Clip(lattice, clip_fc, r'memory\clipped_lattice').getOutput(0)

def lattice_to_lines(lattice):
  return arcpy.management.FeatureToLine(lattice, r'memory\fishnet_lines', "", "NO_ATTRIBUTES").getOutput(0)

def merge_lattice_with_lines(lattice_lines, internal_lines):
  print(f"Merging lattice and diagonals...")
  result = os.path.join(TARGET_GDB, f'lattice_lines_{unique_id()}')
  return arcpy.management.Merge([lattice_lines, internal_lines], result)

def main():
  lattice = create_lattice(AOI)
  clipped_lattice = clip_lattice(lattice, AOI)
  diagonal_lines = create_square_diagonals(clipped_lattice)
  lattice_lines = lattice_to_lines(clipped_lattice)
  merged_lines = merge_lattice_with_lines(lattice_lines, diagonal_lines)
 
if __name__ == "__main__":
  main()