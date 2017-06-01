# vim: set expandtab shiftwidth=2 softtabstop=2:

# util
# A set of handy functions used in both tool and command

from TEMPy.ProtRep_Biopy import BioPy_Structure, BioPyAtom
from TEMPy.EMMap import Map 

def chimera_to_tempy_atom(atom, serial):
  """ Convert one of the ChimeraX atoms to a Tempy Style one """

  ta = BioPyAtom([])

  ta.serial = serial
  ta.atom_name = atom.name
  ta.alt_loc = "" # TODO - Not sure what to put here 

  ta.res = atom.residue.name
  ta.chain = atom.chain_id
  ta.res_no = atom.residue.number
  ta.model = atom.chain_id # PDB number?
  ta.icode = "" # TODO - Not sure about this
  #if atom.is_disordered()==1:
  #    self.icode = "D"
       # 1 if the residue has disordered atoms
  #            self.icode = pdbString[26].strip()#code for insertion residues
  #             # Starting co-ordinates of atom.
  ta.init_x = atom.coord[0]
  ta.init_y = atom.coord[1]
  ta.init_z = atom.coord[2]

  ta.x = atom.coord[0]
  ta.y = atom.coord[1]
  ta.z = atom.coord[2]

  ta.occ = atom.occupancy
  ta.temp_fac = atom.bfactor
  ta.elem = atom.element.name
  ta.charge=""  
  ta.record_name = "HETATM"
  # TODO not sure about this
  if atom.in_chain:
    ta.record_name = "ATOM"

  #Mass of atom as given by atomicMasses global constant. Defaults to 1.
  ta.mass = atom.element.mass

  ta.vdw = 1.7
  ta.isTerm = False
  ta.grid_indices = []

  return ta

def chimera_to_tempy_map(cmap):
  """ Convert a Chimera Map into a tempy one """
  
  apix = 1.0 # Not sure what tis is
  origin = cmap.data.origin  
  
  # Grid_Data is the ChimeraX class we are after here
  # TODO - cmap.matrix *should* return a numpy array which is what we are after I *think*
  # Don't know if its actually correct
  actual_map = cmap.matrix()
  mt = Map(actual_map, origin, apix, "nofilename") 

  return mt


