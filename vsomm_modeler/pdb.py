from collections import Counter, OrderedDict
import numpy as np

__author__ = "Yerko Escalona"
__version__ = "1.0"

class PDB():
    """
    Simple PDB class

    =============  ============  ===========  =============================================
    COLUMNS        DATA  TYPE    FIELD        DEFINITION
    =============  ============  ===========  =============================================
    1 -  6         Record name   "CRYST1"
    7 - 15         Real(9.3)     a              a (Angstroms).
    16 - 24        Real(9.3)     b              b (Angstroms).
    25 - 33        Real(9.3)     c              c (Angstroms).
    34 - 40        Real(7.2)     alpha          alpha (degrees).
    41 - 47        Real(7.2)     beta           beta (degrees).
    48 - 54        Real(7.2)     gamma          gamma (degrees).

    1 -  6         Record name   "ATOM  "
    7 - 11         Integer       serial       Atom  serial number.
    13 - 16        Atom          name         Atom name.
    17             Character     altLoc       Alternate location indicator.
    18 - 21        Residue name  resName      Residue name.
    22             Character     chainID      Chain identifier.
    23 - 26        Integer       resSeq       Residue sequence number.
    27             AChar         iCode        Code for insertion of residues.
    31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    55 - 60        Real(6.2)     occupancy    Occupancy.
    61 - 66        Real(6.2)     tempFactor   Temperature  factor.
    67 - 76        String        segID        (unofficial CHARMM extension ?)
    77 - 78        LString(2)    element      Element symbol, right-justified.
    79 - 80        LString(2)    charge       charge  on the atom.
    =============  ============  ===========  ====================================
    """

    def __init__(self):
        self.__version__ = "1.0"
        # self.atoms = 0
        self.serials = []
        self.names = []
        self.altlocs = []
        self.resnames = []
        self.chainids = []
        self.resseqs = []  # for proteins
        # self.icodes = []  # for proteins
        self.coords = np.empty(shape=[0, 3])
        self.occupancies = []
        self.tempfactors = []
        # self.segids = []  # for VMD
        self.elements = []
        # self.charge = []  # for VMD

    def parse_list(self, pdblist):
        for line in pdblist:
            # self.atoms += 1
            # self.serials.append(self.atoms)
            self.serials.append(line[6:11])
            self.names.append(line[12:16].strip())
            self.altlocs.append(line[16:17].strip())
            self.resnames.append(line[17:21].strip())
            self.chainids.append(line[21:22].strip())
            self.resseqs.append(line[22:26])
            self.coords = np.vstack(
                [self.coords, [float(line[30:38]), float(line[38:46]), float(line[46:54])]])

            self.occupancies.append(line[54:60])
            self.tempfactors.append(line[60:66])  # AKA bfactor

            # self.segids.append(line[66:76].strip()) # VMD format
            self.elements.append(line[76:78].strip())
            # self.atomtypes.append(line[76:78].strip())

    def __add__(self, other):
        output = PDB()

        output.serials = self.serials + other.serials
        output.names = self.names + other.names
        output.altlocs = self.altlocs + other.altlocs
        output.resnames = self.resnames + other.resnames
        output.chainids = self.chainids + other.chainids
        output.resseqs = self.resseqs + other.resseqs
        # output.icodes = # self.icodes = []  # for proteins
        output.coords = self.coords = np.vstack((self.coords, other.coords))
        output.occupancies = self.occupancies + other.occupancies
        output.tempfactors = self.tempfactors + other.occupancies
        # output.segids = # self.segids = []  # for VMD
        output.elements = self.elements + other.elements

        return output

    def output(self):
        # pdb format
        # ATOM = "ATOM  {serial:5d} {name:<4s}{altLoc:<1s}{resName:<4s}"
        #       "{chainID:1s}{resSeq:4d}{iCode:1s}"
        #       "   {pos[0]:8.2f}{pos[1]:8.2f}{pos[2]:8.2f}{occupancy:6.2f}"
        #       "{tempFactor:6.2f}      {segID:<4s}{element:>2s}\n"

        # custom pdb format
        ATOM = ("ATOM  {serial:5} {name:<4s}{altLoc:<1s}{resName:<4s}"
                "{chainID:1s}{resSeq:4} "
                "   {coord[0]:8.2f}{coord[1]:8.2f}{coord[2]:8.2f}{occupancy:8s}"
                "{tempFactor:8s}          {element:>2s}\n")
        output = []

        for i in range(len(self.serials)):
            vals = {}

            # add truncate function in case of a serial number > 9999
            vals['serial'] = self.serials[i]
            vals['name'] = self.names[i]
            vals['altLoc'] = self.altlocs[i][:1]
            vals['resName'] = self.resnames[i][:4]
            vals['chainID'] = self.chainids[i][:1]
            vals['resSeq'] = self.resseqs[i]
            # vals['iCode'] = icodes[i][:1] # for proteins
            # don't take off atom so conversion works
            vals['coord'] = self.coords[i]
            vals['occupancy'] = self.occupancies[i]
            vals['tempFactor'] = self.tempfactors[i]
            # vals['segID'] = segids[i][:4] # VMD format
            vals['element'] = self.elements[i]

            output.append(ATOM.format(**vals))

        return output


def get_mol_coordinates(filename):
    "Obtain data from the PDB file"

    models = OrderedDict()
    modelflag = False
    resname = "resname"
    with open(filename, 'r') as pdbfile:
        for line in pdbfile:
            if line.startswith('MODEL'):
                model = []
                modelflag = True
            elif line.startswith('TER') or line.startswith("\n"):
                continue
            elif modelflag and not line.startswith('ENDMDL'):
                model.append(line)
                resname = line[17:23]
                resname = resname.replace(" ", "")
            elif line.startswith('ENDMDL'):
                modelflag = False
                pdb = PDB()
                pdb.parse_list(model)
                models[resname] = pdb

    return models
