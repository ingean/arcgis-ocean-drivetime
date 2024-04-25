from datetime import datetime

class Log:
  def __init__(self, **options) -> None:
    self.timeformat = options['timeformat'] if 'timeformat' in options  else '%H:%M:%S'
    self.debugging = options['debug'] if 'debug' in options else False
    self.arcgis = options['arcgis'] if 'arcgis' in options else False
    if self.arcgis:
      import arcpy
      self.arcpy = arcpy

  def msg(self, level:str , message: str):
    if level != 'Debug' or self.debugging == True:
      timestamp = datetime.strftime(datetime.now(), self.timeformat)
      if self.arcgis:
        self.arcpy.AddMessage(f'{timestamp} {level}: {message}')
      else:
        print(f'{timestamp} {level}: {message}')

  def debug(self, message:str):
    self.msg('Debug', message)

  def info(self, message:str):
    self.msg('Info', message)
  
  def warning(self, message:str):
    self.msg('Warning', message)

  def error(self, message:str):
    self.msg('Error', message)
