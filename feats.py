#!/usr/bin/env python

from __future__ import print_function

_usage="""
   ShowFeats [optional args] <filenames>
"""

import math
#set up the logger:
from rdkit import RDLogger as logging
logger = logging.logger()
logger.setLevel(logging.INFO)

from rdkit import Geometry
from rdkit.Chem.Features import FeatDirUtilsRD as FeatDirUtils

def _getVectNormal(v,tol=1e-4):
  if math.fabs(v.x)>tol:
    res = Geometry.Point3D(v.y,-v.x,0)
  elif math.fabs(v.y)>tol:
    res = Geometry.Point3D(-v.y,v.x,0)
  elif math.fabs(v.z)>tol:
    res = Geometry.Point3D(1,0,0)
  else:
    raise ValueError('cannot find normal to the null vector')
  res.Normalize()
  return res

_canonArrowhead=None
def _buildCanonArrowhead(headFrac,nSteps,aspect):
  global _canonArrowhead
  startP = RDGeometry.Point3D(0,0,headFrac)
  _canonArrowhead=[startP]

  scale = headFrac*aspect
  baseV = RDGeometry.Point3D(scale,0,0)
  _canonArrowhead.append(baseV)

  twopi = 2*math.pi
  for i in range(1,nSteps):
    v = RDGeometry.Point3D(scale*math.cos(i*twopi),scale*math.sin(i*twopi),0)
    _canonArrowhead.append(v)

def ShowMolFeats(mol,factory,radius=0.5,confId=-1,showOnly=True,
                 name='',transparency=0.0,colors=None,excludeTypes=[],
                 useFeatDirs=True,featLabel=None,dirLabel=None,includeArrowheads=True,
                 writeFeats=False,showMol=True,featMapFile=False):
  if not name:
    if mol.HasProp('_Name'):
      name =  mol.GetProp('_Name')
    else:
      name = 'molecule'

  molFeats=factory.GetFeaturesForMol(mol)
  if not featLabel:
    featLabel='%s-feats'%name
  if not dirLabel:
    dirLabel=featLabel+"-dirs"

  for i,feat in enumerate(molFeats):
    family=feat.GetFamily()
    if family in excludeTypes:
      continue
    pos = feat.GetPos(confId)
    nm = '%s(%d)'%(family,i+1)

    if writeFeats:
      aidText = ' '.join([str(x+1) for x in feat.GetAtomIds()])
      print('%s\t%.3f\t%.3f\t%.3f\t1.0\t# %s'%(family,pos.x,pos.y,pos.z,aidText))

    if featMapFile:
      print("  family=%s pos=(%.3f,%.3f,%.3f) weight=1.0"%(family,pos.x,pos.y,pos.z),end='',file=featMapFile)

    if featMapFile:
      aidText = ' '.join([str(x+1) for x in feat.GetAtomIds()])
      print('# %s'%(aidText),file=featMapFile)

# --- ----  --- ----  --- ----  --- ----  --- ----  --- ----
import sys,os,getopt
from rdkit import RDConfig
from optparse import OptionParser
parser=OptionParser(_usage)

parser.add_option('-x','--exclude',default='',
                  help='provide a list of feature names that should be excluded')
parser.add_option('-f','--fdef',default=os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef'),
                  help='provide the name of the feature definition (fdef) file.')
parser.add_option('--noDirs','--nodirs',dest='useDirs',default=True,action='store_false',
                  help='do not draw feature direction indicators')
parser.add_option('--noHeads',dest='includeArrowheads',default=True,action='store_false',
                  help='do not draw arrowheads on the feature direction indicators')
parser.add_option('--noClear','--noClear',dest='clearAll',default=False,action='store_true',
                  help='do not clear PyMol on startup')
parser.add_option('--noMols','--nomols',default=False,action='store_true',
                  help='do not draw the molecules')
parser.add_option('--writeFeats','--write',default=False,action='store_true',
                  help='print the feature information to the console')
parser.add_option('--featMapFile','--mapFile',default='',
                  help='save a feature map definition to the specified file')
parser.add_option('--verbose',default=False,action='store_true',
                  help='be verbose')

if __name__=='__main__':
  from rdkit import Chem
  from rdkit.Chem import AllChem
  from rdkit.Chem.PyMol import MolViewer

  options,args = parser.parse_args()
  if len(args)<1:
    parser.error('please provide either at least one sd or mol file')

  try:
    fdef = open(options.fdef,'r').read()
  except IOError:
    logger.error('ERROR: Could not open fdef file %s'%options.fdef)
    sys.exit(1)

  factory = AllChem.BuildFeatureFactoryFromString(fdef)

  if options.writeFeats:
     print('# Family   \tX    \tY   \tZ   \tRadius\t # Atom_ids')

  if options.featMapFile:
    if options.featMapFile=='-':
      options.featMapFile=sys.stdout
    else:
      options.featMapFile=file(options.featMapFile,'w+')
    print("BeginPoints", file=options.featMapFile)

  i = 1
  for midx,molN in enumerate(args):
    if molN!='-':
      featLabel='%s_Feats'%molN
    else:
      featLabel='Mol%d_Feats'%(midx+1)

    dirLabel=featLabel+"-dirs"

    if molN != '-':
      try:
        ms = Chem.SDMolSupplier(molN)
      except:
        logger.error('Problems reading input file: %s'%molN)
        ms = []
    else:
      ms = Chem.SDMolSupplier()
      ms.SetData(sys.stdin.read())

    for m in ms:
      nm = 'Mol_%d'%(i)
      if m.HasProp('_Name'):
        nm += '_'+m.GetProp('_Name')
      if options.verbose:
        if m.HasProp('_Name'):
          print("#Molecule: %s"%m.GetProp('_Name'))
        else:
          print("#Molecule: %s"%nm)
      ShowMolFeats(m,factory,transparency=0.25,excludeTypes=options.exclude,name=nm,
                   showOnly=False,
                   useFeatDirs=options.useDirs,
                   featLabel=featLabel,dirLabel=dirLabel,
                   includeArrowheads=options.includeArrowheads,
                   writeFeats=options.writeFeats,showMol = False,
                   featMapFile=options.featMapFile)
      i += 1
      if not i%100:
        logger.info("Done %d poses"%i)

  if options.featMapFile:
    print("EndPoints",file=options.featMapFile)


  sys.exit(0)
