""" module hfile.py

A set of routines to help read histogram files and convert them to
numpy arrays that can then be manipulated before plotting. This module
also includes some helper routines for manipulating the results
(e.g. rebinning) and searching through the file.

"""
from __future__ import division
from __future__ import print_function
from builtins import range
# past is included in the futures pip install
from past.builtins import basestring
from builtins import object

# This file is based on an earlier version written by Gavin P. Salam
# (2011)

# Comments may appear excessive in places. This a consequence of this
# file being one of the author's early attempts at python programming.



import numpy      # for arrays
import string    
import re         # for regular expressions
import io   
#import scipy.interpolate
#import sys

# a few variables that users can set(?)
default_lw = 3
default_border_lw = 1

default_encoding='utf-8'


#----------------------------------------------------------------------
class Error(Exception):
    """Base class for exceptions in this module."""
    def __init__(self, msg):
      self.msg = msg

#----------------------------------------------------------------------
def array(file, regexp=None, fortran=False):
  "Same as get_array"
  return get_array(file, regexp, fortran)

#----------------------------------------------------------------------
def open_any(filename, mode=None):
    """Return either a gzipped file or a normal file, based on the extension"""
    if (len(filename)>3 and filename[-3:] == '.gz'):
        import gzip
        return gzip.GzipFile(filename, mode)
    else:
        return open(filename, mode)

def decode(line):
    if (isinstance(line,bytes)): return line.decode(default_encoding)
    else: return line
    
#----------------------------------------------------------------------
def get_2darray(file,regexp=None):
    """
    For data in the form
    
     1  1   1
     1  2   5
    
     2  1   5
     2  2   8
    
    the resulting array will have the following contents
    
         array[0,0,:] = [1,1,1]
         array[0,1,:] = [1,2,5]
         ...
         array[1,1,:] = [2,2,8]
    
    The array is deemed to terminate when there are two blank lines
    """
    # if file is a string, assume it's a filename, otherwise a filehandle
    if (isinstance(file,basestring)) : file = open_any(file, 'r')
    if (regexp != None)              : search(file,regexp)
    
    # read in the first slice [0th index = 0]
    slices = []
    slices.append(get_array(file))
    # this will leave us in a situation where we've read in the first blank
    # line that follows a block of numbers
    
    # so we try to read in the next line, and if it's the end of the
    # file, or it's a second blank line, then we deem that we have
    # reached the end of the block, and try to go back to the original
    # place in the file.
    while True:
        # record where we are
        pos = file.tell()
        # read in a line
        line = file.readline()
        #print (line)
        if (not line or line.strip() == ""): break
        # go back to the start of the line
        file.seek(pos)
        slices.append(get_array(file))
    
    result = numpy.ndarray((len(slices), slices[0].shape[0], slices[0].shape[1]))
    for i in range(len(slices)):
        result[i,:,:] = slices[i]
    
    return result
      
#----------------------------------------------------------------------
def get_array(file, regexp=None, fortran=False):
  """
  Returns a 2d array that contains the next block of numbers in this
  file (which can be a filehandle or a filename) 
  
  - if a regexp is provided, then first that is searched for
  - Subsequently, any line starting with a non-numeric character is ignored
  - The block ends the first time a blank line is encountered
  
  For data in the form
  
   1  1
   2  4
   3  9
  
  the resulting array will have the following contents
  
       array[:,0] = [1, 2, 3]
       array[:,1] = [1, 4, 9] 
  
  """
  # NB: various potential improvements may be possible here, e.g.
  # - using # things like "isspace" to detect blank lines,
  # - using numpy.loadtxt # to do the string to array conversion,
  # - adding option of getting multiple blocks (e.g. as a python array/list)
  # - and adding option of reading gnuplot-style 3d data
  
  # if file is a string, assume it's a filename, otherwise a filehandle
  if (isinstance(file,basestring)) : file = open_any(file, 'r')
  if (regexp != None)              : search(file,regexp)

  # handle case where we numbers such as 0.4d3 (just replace d -> e)
  fortranRegex = re.compile(r'd', re.IGNORECASE)
  
  lines = []    # temporary store of lines, before conversion to array
  started = False
  while True:
    line = decode(file.readline())
    if (not line)               : break        # empty line = end-of-file
    line = line.rstrip()                       # strips trailing blanks, \n
    line = line.lstrip()                       # strips leading blanks
    if (not line) :
      if (started) : break                     # empty line = end-of-block
      else         : continue                  
    if (not re.match('[-0-9]', line)) : continue
    if (fortran): line = re.sub(fortranRegex,'e',line) # handle fortran double-prec
    lines.append(line)                         # collect the line
    started = True
  # do some basic error checking
  if (len(lines) < 1):
    raise Error("Block in get_array had 0 useful lines")
  # now we know the size, transfer the information to a numpy ndarray
  ncol = len(lines[0].split())                
  num_array = numpy.empty( (len(lines), ncol) )
  for i in range(len(lines)):
    num_array[i,:] = lines[i].split()
  return num_array


#----------------------------------------------------------------------
class XSection(object):
    """\
    Class that contains the result of an "xsc" entry in the file.
    (This is intended to work with files from the CSS framework)
    
    Relevant members are:
    obj.xsc : the cross section
    obj.err : its error
    obj.nentries: the number of entries of which it is built
    obj.name: the name
    obj.line: the full line
    """
    def __init__(self, line):
        self.line = line
        # got this from
        # http://stackoverflow.com/questions/6508043/regular-expression-to-find-any-number-in-a-string
        number_regex = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        # use the last 3 numbers 
        numbers = number_regex.findall(line)
        self.xsc = float(numbers[-3])
        self.err = float(numbers[-2])
        self.nentries = float(numbers[-1])
        # recall that "group" gets the first matching group here...
        self.units = re.search(r"([a-z]+) \(n entries",line).group((1))
        self.name = re.search(r"# (.*): xsc =",line).group((1))
        return None

    def __str__(self):
        return "{0} +- {1} {2} (n entries = {3})".format(self.xsc,self.err,self.units,self.nentries)

#----------------------------------------------------------------------
def get_xsection(file,regexp):
    if (isinstance(file,basestring)) : file = open_any(file, 'r')
    fullregexp = regexp+r'.*xsc ='
    line = search(file, fullregexp, return_line = True)
    #
    result = XSection(line)
    return result

#----------------------------------------------------------------------
def search(filehandle, regexp, return_line=False):
  """ looks through the file described by handle until it finds the regexp
  that's mentioned"""
  regexp_compiled = re.compile(regexp)
  while True:
    line = decode(filehandle.readline())
    if (line == "") : break        # empty line = end-of-file
    #if (re.search(regexp,line)) : 
    if (regexp_compiled.search(line)) : 
        if (return_line): return line
        else            : return filehandle
  return None


#----------------------------------------------------------------------
def reformat(*columns, **keyw):
  """returns a string containing each of the columns placed
      side by side. If instead of columns, 2d numpy arrays
      are supplied, then these are output sensibly too

      For now, it assumes that columns are numpy arrays,
      but this could in principle be relaxed.

      Among the keyword arguments, the only one currently supported is
      "format", which should be a string used to format the output

  """

  ncol = len(columns)
  shapes = []
  ncols  = []
  for i in range(ncol):
    shapes.append(columns[i].shape)
    if    (len(shapes[i]) == 1): ncols.append(0)
    elif  (len(shapes[i]) == 2): ncols.append(shapes[i][1])
    else: raise Error(" a 'column' appears not to be 1 or 2-dimensional" )

  # lazily assume that all lengths are the same
  nlines = shapes[0][0]

  output = io.BytesIO()
  if ("format" in keyw):
    frm=keyw["format"]
    
    for i in range(nlines) :
      for j in range(ncol):
        if (ncols[j] == 0): print(frm.format(columns[j][i]), end=' ', file=output)    # trailing comma kills newline
        else: 
          for k in range (ncols[j]):
            print(frm.format(columns[j][i,k]), end=' ', file=output)
      print(file=output)
  else:
    for i in range(nlines) :
      for j in range(ncol):
        if (ncols[j] == 0): print(columns[j][i], end=' ', file=output)    # trailing comma kills newline
        else: 
          for k in range (ncols[j]):
            print(columns[j][i,k], end=' ', file=output)
      print(file=output)

      
  return output.getvalue()

#----------------------------------------------------------------------
def cumul(array, nx=3, renorm=False, reverse = False):
  """
  Returns the accumulated version of the array. Currently assumes uniform
  spacing.

  \param nx      indicates how many columns of x-axis info there are. If
                 nx = 1, the spacing must be uniform and the x axis is adjusted
                 so as to be shifted by half a bin in the right direction

  \param renorm  if True, it assumes that we are dealing with a 
                 differential distribution (so sum of bins must be divided 
                 by rebin); NB: this is buggy if bin spacings not uniform

  \param reverse accumulates negatively from the high end downwards if
                 true. The high bin will always be at xmax+dx/2,
                 regardless of the value of reverse, which means that
                 if reverse is True, then the high-end y bins will be
                 zero. This facilitates addition of results with and
                 without reverse, since the x bins should always match
  """
  ny = len(array[:,0])
  if (ny < 2) : raise Error("number of bins must be >= 2")

  result = numpy.empty(array.shape)
  # copy x and y axes
  result[:,:nx] = array[:,:nx]
  if (reverse):
    # some care is needed to make sure that we align the bins properly
    result[ny-1,nx:] = 0
    result[:ny-1,nx:] = -array[ny-1:0:-1,nx:].cumsum(axis=0)[::-1]
  else:
    result[:,nx:] = array[:,nx:].cumsum(axis=0)

  # sort out shifting of bins and renormalisation
  if (nx == 1): 
    dx = array[1,0]-array[0,0]
    result[:,0] += 0.5*dx
    if (renorm): result[:,nx:] *= dx  
  else:
    if (renorm): result[:,nx:] *= (result[:,nx-1]-result[:,0]).reshape(ny,1)
  
  # finally handle the reverse accumulation, keeping identical binning
  # to what we would have had in the forward direction (this
  # effectively loses us one bin of information, but it adds
  # convenience for adding things produced in different directions)
  #if (reverse): result[:,nx:] = -result[-1,nx:] + result[:,nx:]
    
  return result

#----------------------------------------------------------------------
def rebin(array, rebin=2, nx=3, renorm=False, 
                 min=-1e300, max=1e300):
  """ Returns a rebinned version of an array. 
  
  \param rebin   indicates how many bins to combine
  \param nx      indicates how many columns of x-axis info there are
  \param renorm  if True, it assumes that we are dealing with a 
                 differential distribution (so sum of bins must be divided 
                 by rebin)
  \param min     starts only from bins >= min
  \param max     goes only up to bins <= max
  """
  ilo = 0
  ihi = len(array[:,0])-1

  if (ihi < 1) : raise Error("number of bins must be >= 2")
  
  if   (nx == 1) : minslice = array[:,0]; maxslice = array[:,0]
  elif (nx == 2) : minslice = array[:,0]; maxslice = array[:,1]
  elif (nx == 3) : minslice = array[:,0]; maxslice = array[:,2]
  else : raise Error("nx is not 1, 2 or 3")
  
  # get the bounds
  for i in range(len(minslice)):
    if (min <= minslice[i]): ilo = i; break;
  for i in range(len(minslice)-1,-1,-1):
    if (max >= maxslice[i]): ihi = i; break;

  # now figure out how many bins we're going to need
  nbins = int(ihi+1-ilo)//rebin
  #print ilo, ihi, minslice[ilo], maxslice[ihi]
  #print nbins
  rebinned = numpy.empty( (nbins, len(array[0,:])) )
  
  # and now do the rebinning
  ir = 0
  for i in range(ilo, ilo+nbins*rebin, rebin):
    # first handle the x axes
    if   (nx == 1): rebinned[ir, 0] = sum(array[i:i+rebin,0])/rebin
    elif (nx >= 2): 
      rebinned[ir, 0   ] = array[i        , 0    ]
      rebinned[ir, nx-1] = array[i+rebin-1, nx-1 ]
      if (nx == 3) : rebinned[ir, 1] = 0.5*(rebinned[ir,0]+rebinned[ir,2])
    # then the remaining ones
    rebinned[ir, nx:] = numpy.sum(array[i:i+rebin, nx:], axis=0)
    ir += 1

  if (renorm) : rebinned[:,nx:] /= rebin
  return rebinned

#----------------------------------------------------------------------
def extended_binary(dict_or_list_of_arrays, binary_fn, subset = None):
  """
  Return an array each of whose elements is given by the repeated
  pairwise application of binary_fn on the corresponding elements in the
  dict_or_list_of_arrays; if subset is not None, then only the keys in
  the subset are used
  """
  
  if (isinstance(dict_or_list_of_arrays, dict)):
    if (subset is None): subset = list(dict_or_list_of_arrays.keys())
    arrays = [dict_or_list_of_arrays[key] for key in subset]
  else:
    arrays = dict_or_list_of_arrays

  first = True
  for array in arrays:
    if (first): result = array
    else      : result = binary_fn(result,array)
    first = False

  return result

#----------------------------------------------------------------------
def maximum(dict_or_list_of_arrays, subset = None):
  """
  Return an array each of whose elements is given by the maximum of
  the corresponding elements in the dict_or_list_of_arrays; if subset is not
  None, then only the keys in the subset are used
  """
  return extended_binary(dict_or_list_of_arrays, numpy.maximum, subset = subset)

#----------------------------------------------------------------------
def minimum(dict_or_list_of_arrays, subset = None):
  """
  Return an array each of whose elements is given by the minimum of
  the corresponding elements in the dict_or_list_of_arrays; if subset is not
  None, then only the keys in the subset are used
  """
  return extended_binary(dict_or_list_of_arrays, numpy.minimum, subset = subset)


#----------------------------------------------------------------------
def minmax(dict_or_list_of_arrays, subset = None):
  """ Return a tuple such that [0][...] is the min
  and [1][...] is the max across the dictionary of arrays
  """
  return (minimum(dict_or_list_of_arrays,subset),
          maximum(dict_or_list_of_arrays,subset))

#------------------------------------------------------------
# should go into a base class?
def opt_or_default(keyw, key, default = None, delete=False):
  """
  If keyw[key] exists, returns it and deletes keyw[key], otherwise
  returns the default value
  """
  if key in keyw:
    result = keyw[key]
    if delete: del keyw[key]
    return result
  else: return default

