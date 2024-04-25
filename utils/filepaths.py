import os

class Path:
  """Utility class for working with file paths, names and extensions.

  Args:
      path (str): String representing a path, path with file name or just a file name.
      path (list): list of strings representing a path, filename and file extension.
      
  Attributes:
      path (str): Returns the path (but not the file name).
      fullpath (str): Returns path and file name with extension.
      fullfilename (str): Returns the file name with extension.
      filename (str): Returns the file name with no extension.
      extension (str): Returns the file extension if exists.
      exists (boolean): Returns True if it is an existing file or directory.
      isfile (boolean): Returns True if it is an existing file.
      isdir (boolean): Returns True if it is an existing directory

  """
  def __init__(self, path):
    if type(path) is list:
      fne = path[1]
      if len(path) > 2:
        e = path[2].replace('.', '')
        fne = f'{path[1]}.{e}'

      path = os.path.join(path[0], fne)

    fn = os.path.basename(path)
    p = os.path.dirname(path)
    fn_parts = fn.split('.')
    
    self.path = p
    self.fullpath = path
    self.fullfilename = fn
    self.filename = fn_parts[0]
    self.extension = fn_parts[1] if len(fn_parts) > 1 else None
    self.exists = os.path.exists(path)
    self.isfile = os.path.isfile(path)
    self.isdir = os.path.isdir(self.path)



""" # Tests
#p = Path(r'D:\Data\Imagery\Ortofoto\P10_200407-Orthomosaic.tiff')
#p = Path(r'D:\Data\Imagery\Ortofoto\P10_200407-Orthomosaic')
#p = Path(r'D:\Data\Imagery\'')
p = Path([r'D:\Data\Imagery', 'P10_200407-Orthomosaic', '.tiff'])

print(f'Path: {p.path}')
print(f'Full path: {p.fullpath}')
print(f'Full filename: {p.fullfilename}')
print(f'Filename: {p.filename}')
print(f'File extension: {p.extension}')
print(f'Full path exists: {p.exists}')
print(f'Input is file: {p.isfile}')
print(f'Input is dir: {p.isdir}') """
