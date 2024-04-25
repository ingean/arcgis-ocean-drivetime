import os
import requests
import arcpy
from utils.arcpy_dataset import create_featureclass, add_fields

TARGET_GDB = r'D:\Data\Redningsselskapet\GeoAnalyse\RS_GeoAnalyse.gdb'
stations = os.path.join(TARGET_GDB, 'rs_stasjoner')

URL = 'https://prod-rsfeed-xml2json-proxy.rs-marine-services.rs.no/prefetch/getstations'


def get_stations():
  res = requests.get(URL)
  return res.json()

def add_fields_from_response(fc):
  field_defs = [
    ["name", "TEXT", "Navn"],
    ["type", "TEXT", "Type"],
    ["region", "TEXT", "Region"],
    ["phone", "TEXT", "Telefonnummer"],
    ["address", "TEXT", "Adresse"],
    ["ziplocation", "TEXT", "Poststed"],
    ["zipcode", "TEXT", "Postnummer"],
    ["coordinate_type", "TEXT", "Koordinattype"],
    ["latitude", "DOUBLE", "Breddegrad"],
    ["longitude", "DOUBLE", "Lengdegrad"],
    ["descr_approach_sea", "TEXT", "Ankomst, sjø"],
    ["descr_approach_land", "TEXT", "Ankomst, land"],
    ["descr_depths", "TEXT", "Dybdeforhold"],
    ["descr_special_risks", "TEXT", "Risikofaktorer"],
    ["descr_bunkers", "TEXT", "Bunkre"],
    ["descr_dangerous_materials_delivery", "TEXT", "Leveranse av farlig avfall"],
    ["descr_other_information", "TEXT", "Annen informasjon"],
    ["descr_communication", "TEXT", "Kommunikasjon"],
    ["descr_mooring", "TEXT", "Fortøyning"],
    ["descr_tide", "TEXT", "Tidevann"],
    ["descr_warehouse", "TEXT", "Varehus"],
    ["descr_insurance", "TEXT", "Forsikring"],
    ["descr_television", "TEXT", "TV"],
    ["descr_power_supplier", "TEXT", "Strømselskap"],
    ["descr_shipwreck_area", "TEXT", "Område med skipsvrak"],
    ["descr_other_ports", "TEXT", "Andre havner"],
    ["descr_cooperation", "TEXT", "Samarbeidspartnere"],
    ["descr_district", "TEXT", "Distrikt"],
    ["descr_deals", "TEXT", "Avtaler"],
    ["descr_recreation", "TEXT", "Fritidsaktiviteter"],
    ["municipality_web", "TEXT", "Kommunale nettsider"]
  ]
  arcpy.management.AddFields(stations, field_defs)

def add_data_from_response(fc, response):
  fields = ['SHAPE@', 'name', 'type']
  with arcpy.da.InsertCursor(fc, fields) as cursor:
    for s in response:
      geom = arcpy.PointGeometry(arcpy.Point(s["longitude"],s["latitude"]),arcpy.SpatialReference(4326)).projectAs(arcpy.SpatialReference(25833))
      cursor.insertRow((geom, s["name"], s["type"]))

def main():
  response = get_stations()
  create_featureclass(fc=stations, delete_existing=True)
  add_fields_from_response(stations)
  add_data_from_response(stations, response["stations"])

if __name__ == "__main__":
  main()