from re import finditer
from collections import namedtuple

__author__ = "Yerko Escalona"
__version__ = "1.0"


class TOP(object):

    class _Solute(object):
        """Inner class to store information about a molecule type."""

        def __init__(self):
            self.atomNumbers = []
            self.resNumbers = []
            self.atomNames = []
            self.atomTypes = []
            self.atomMasses = []
            self.atomCharges = []
            self.chargeGroups = []
            self.excludedAtoms = []
            self.excluded14Atoms = []
            self.molecules = []

    class _Solvent(object):
        """Inner class to store information about a molecule type."""
        def __init__(self):
            self.atomNumbers = []
            self.atomNames = []
            self.atomTypes = []
            self.atomMasses = []
            self.atomCharges = []

    class _MoleculeType(object):
        """Inner class to store information about a molecule type."""
        def __init__(self):
            self.atoms = []
            self.bonds = []
            self.angles = []
            self.impropers = []
            self.dihedrals = []

    def _split_iter(self, string):
        """the same as s.split() but it makes a generator out of it"""
        for match in finditer("\S+", string):
            yield match.group(0)

    def _processBlocks(self, filename):
        with open(filename) as file:
            lines = ""
            for line in file:
                if line.startswith('#'):
                    continue
                elif '#' in line:
                    # TODO: check possible errors by using sharps
                    line = line[:line.index('#')]
                elif line.startswith('END'):
                    yield lines
                    lines = ""
                else:
                    lines += line

    def _processBlock(self, block):
        title = ""
        for b in block:
            title += b+" "

    def _processAtomTypeName(self, block):
        atomTypes = int(next(block))
        for atomType in range(atomTypes):
            self.atomTypeNames.append(next(block))

    def _processResidueName(self, block):
        residues = int(next(block))
        for residue in range(residues):
            self.residueNames.append(next(block))

    def _processSoluteAtom(self, block):
        atoms = int(next(block))
        for atom in range(atoms):
            self.solute.atomNumbers.append(int(next(block)))
            self.solute.resNumbers.append(int(next(block)))
            self.solute.atomNames.append(next(block))
            self.solute.atomTypes.append(int(next(block)))
            self.solute.atomMasses.append(float(next(block)))
            self.solute.atomCharges.append(float(next(block)))
            self.solute.chargeGroups.append(int(next(block)))

            excludedAtoms = []
            excluded_atoms = int(next(block))
            for exclude_atom in range(excluded_atoms):
                excludedAtoms.append(int(next(block)))
            self.solute.excludedAtoms.append(excludedAtoms)

            excluded14Atoms = []
            excluded_atoms14 = int(next(block))
            for exclude_atom14 in range(excluded_atoms14):
                excluded14Atoms.append(int(next(block)))
            self.solute.excluded14Atoms.append(excluded14Atoms)

    def _processBond(self, block):
        bonds = int(next(block))
        for bond in range(bonds):
            i = int(next(block))
            j = int(next(block))
            bondType = int(next(block))
            self.MoleculeType.bonds.append(self.Bond(i, j, bondType))

    def _processAngle(self, block):
        angles = int(next(block))
        for angle in range(angles):
            i = int(next(block))
            j = int(next(block))
            k = int(next(block))
            angleType = int(next(block))
            self.MoleculeType.angles.append(self.Angle(i, j, k, angleType))

    def _processImproper(self, block):
        impropers = int(next(block))
        for improper in range(impropers):
            i = int(next(block))
            j = int(next(block))
            k = int(next(block))
            l = int(next(block))
            improperType = int(next(block))
            self.MoleculeType.impropers.append(self.Improper(i, j, k, l, improperType))

    def _processDihedral(self, block):
        dihedrals = int(next(block))
        for dihedral in range(dihedrals):
            i = int(next(block))
            j = int(next(block))
            k = int(next(block))
            l = int(next(block))
            dihedralType = int(next(block))
            self.MoleculeType.dihedrals.append(self.Dihedral(i, j, k, l, dihedralType))

    def _processSoluteMolecules(self, block):
        molecules = int(next(block))
        for molecule in range(molecules):
            self.soluteMolecules.append(int(next(block)))

    def _processSolventAtom(self, block):
        atoms = int(next(block))
        for atom in range(atoms):
            self.solvent.atomNumbers.append(int(next(block)))
            self.solvent.atomNames.append(next(block))
            self.solvent.atomTypes.append(float(next(block)))
            self.solvent.atomMasses.append(float(next(block)))
            self.solvent.atomCharges.append(float(next(block)))

    def processFile(self, filename):
        """Process one line from a file."""
        blocks = self._processBlocks(filename)

        for block in blocks:
            _currentBlock = self._split_iter(block)
            _currentCategory = next(_currentBlock)

            if _currentCategory == 'ATOMTYPENAME':
                self._processAtomTypeName(_currentBlock)
            if _currentCategory == 'RESNAME':
                self._processResidueName(_currentBlock)
            elif _currentCategory == 'SOLUTEATOM':
                self._processSoluteAtom(_currentBlock)
            elif _currentCategory == 'BONDH':
                self._processBond(_currentBlock)
            elif _currentCategory == 'BOND':
                self._processBond(_currentBlock)
            elif _currentCategory == 'BONDANGLEH':
                self._processAngle(_currentBlock)
            elif _currentCategory == 'BONDANGLE':
                self._processAngle(_currentBlock)
            elif _currentCategory == 'IMPDIHEDRALH':
                self._processImproper(_currentBlock)
            elif _currentCategory == 'IMPDIHEDRAL':
                self._processImproper(_currentBlock)
            elif _currentCategory == 'DIHEDRALH':
                self._processDihedral(_currentBlock)
            elif _currentCategory == 'DIHEDRAL':
                self._processDihedral(_currentBlock)
            elif _currentCategory == 'SOLUTEMOLECULES':
                self._processSoluteMolecules(_currentBlock)
            elif _currentCategory == 'SOLVENTATOM':
                self._processSolventAtom(_currentBlock)

    def generate_subtopologies(self):
        i = 0
        for soluteMolecule in self.soluteMolecules:
            first_gi = i + 1
            last_gi = soluteMolecule + 1
            
            tmp_top = TOP()
            
            tmp_top.atomTypeNames = self.atomTypeNames
            
            tmp_top.solute.atomNumbers = [a-i for a in self.solute.atomNumbers[i:soluteMolecule]]
            tmp_top.solute.atomTypes = self.solute.atomTypes[i:soluteMolecule]
            
            new_resNumbers = []
            for resNumber in self.solute.resNumbers[i:soluteMolecule]:
                new_resNumbers.append(resNumber - self.solute.resNumbers[i] + 1)
            
            tmp_top.solute.resNumbers = new_resNumbers
            
            new_residueNames = []
            for resNumber in set(self.solute.resNumbers[i:soluteMolecule]):
                new_residueNames.append(self.residueNames[resNumber-1])
                
            tmp_top.residueNames = new_residueNames
            
            tmp_top.solute.atomNames = self.solute.atomNames[i:soluteMolecule]
            tmp_top.solute.chargeGroups = self.solute.chargeGroups[i:soluteMolecule]
            tmp_top.solute.atomCharges = self.solute.atomCharges[i:soluteMolecule]
            tmp_top.solute.atomMasses = self.solute.atomMasses[i:soluteMolecule]
            
            new_excludedAtoms = []
            for excludedAtoms in self.solute.excludedAtoms[i:soluteMolecule]:
                new_excludedAtoms.append([a-i for a in excludedAtoms])            
            tmp_top.solute.excludedAtoms = new_excludedAtoms
            
            new_excluded14Atoms = []
            for excluded14Atoms in self.solute.excluded14Atoms[i:soluteMolecule]:
                new_excluded14Atoms.append([a-i for a in excluded14Atoms])
            tmp_top.solute.excluded14Atoms = new_excluded14Atoms
            
            for bond in self.MoleculeType.bonds:
                for gi in range(first_gi, last_gi):
                    if gi == bond.i or gi == bond.j:
                        new_bond = self.Bond(bond.i-i, bond.j-i, bond.type)
                        tmp_top.MoleculeType.bonds.append(new_bond)
                        break
            
            for angle in self.MoleculeType.angles:         
                for gi in range(first_gi, last_gi):
                    if gi == angle.i or gi == angle.j or gi == angle.k:
                        new_angle = self.Angle(angle.i-i, angle.j-i, angle.k-i, angle.type)
                        tmp_top.MoleculeType.angles.append(new_angle)
                        break
                        
            for improper in self.MoleculeType.impropers:
                for gi in range(first_gi, last_gi):
                    if gi == improper.i or gi == improper.j or gi == improper.k or gi == improper.l:
                        new_improper = self.Improper(improper.i-i, improper.j-i, improper.k-i, improper.l-i, improper.type)
                        tmp_top.MoleculeType.impropers.append(new_improper)
                        break
            
            for dihedral in self.MoleculeType.dihedrals:
                for gi in range(first_gi, last_gi):
                    if gi == dihedral.i or gi == dihedral.j or gi == dihedral.k or gi == dihedral.l:
                        new_dihedral = self.Dihedral(dihedral.i-i, dihedral.j-i, dihedral.k-i, dihedral.l-i, dihedral.type)
                        tmp_top.MoleculeType.dihedrals.append(new_dihedral)
                        break
            
            tmp_top.solvent = self.solvent
        
            ###
            
            i = soluteMolecule
            yield tmp_top

    def generate_itp(self, moleculeName="mol"):

        lines = """[ moleculetype ]
; Name            nrexcl
{}     0
""".format(moleculeName)

        lines += """
[ atoms ]
;   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB
"""
    
        cg = 1
        for i, gi in enumerate(self.solute.atomNumbers):
            lines += "{:>6d}{:>11s}{:>7d}{:>7s}{:>7s}{:7d}{:>11s}{:>11s}\n".format(gi,  
                  self.atomTypeNames[self.solute.atomTypes[i] - 1],
                  self.solute.resNumbers[i],
                  self.residueNames[self.solute.resNumbers[i] - 1],
                  self.solute.atomNames[i],
                  cg,
                  str(self.solute.atomCharges[i]), 
                  str(self.solute.atomMasses[i]))

            if self.solute.chargeGroups[i] == 1:
                cg += 1

        lines += """
[ bonds ]
;  ai    aj funct            c0            c1            c2            c3
"""
        for bond in self.MoleculeType.bonds:
            lines += "{:>5d}{:>6d}{:>6s}    {:}\n".format(bond.i, bond.j, '2', 'gb_' + str(bond.type))
        
        lines += """
[ exclusions ]
;  ai    aj
"""
        for i, gi in enumerate(self.solute.atomNumbers):
            for gj in self.solute.excludedAtoms[i]:
                lines += "{:>5d}{:>6d}\n".format(gi, gj)
            
        for i, gi in enumerate(self.solute.atomNumbers):
            for gj in self.solute.excluded14Atoms[i]:
                lines += "{:>5d}{:>6d}\n".format(gi, gj)
                
        lines += """
[ pairs ]
;  ai    aj funct
"""
        for i, gi in enumerate(self.solute.atomNumbers):
            for gj in self.solute.excluded14Atoms[i]:
                lines += "{:>5d}{:>6d}     1\n".format(gi, gj)
                
        lines += """
[ angles ]
;  ai    aj    ak funct            c0            c1            c2            c3        
"""
        
        for angle in self.MoleculeType.angles:
            lines += '{:>5d}{:>6d}{:>6d}{:>6s}    {:}\n'.format(angle.i, angle.j, angle.k,
                                                                '2', 'ga_' + str(angle.type))
        
        lines += """
[ dihedrals ]
;  ai    aj    ak    al funct            c0            c1            c2            c3            c4            c5      
"""
        for dihedral in self.MoleculeType.dihedrals:
            lines += "{:>5d}{:>6d}{:>6d}{:>6d}{:>6s}    {:}\n".format(dihedral.i, dihedral.j, dihedral.k, dihedral.l,
                                                                     '1', 'gd_' + str(dihedral.type))
        for improper in self.MoleculeType.impropers:
            lines += "{:>5d}{:>6d}{:>6d}{:>6d}{:>6s}    {:}\n".format(improper.i, improper.j, improper.k, improper.l,
                                                                     '2', 'gi_' + str(improper.type))
        lines += "\n"
                
        return lines

    def __init__(self):
        self.solute = self._Solute()
        self.solvent = self._Solvent()
        self.MoleculeType = self._MoleculeType()
        self.atomTypeNames = []
        self.residueNames  = []
        self.soluteMolecules = []
        
        self.Bond = namedtuple('Bond', ['i', 'j', 'type'])
        self.Angle = namedtuple('Angle', ['i', 'j', 'k', 'type'])
        self.Improper = namedtuple('Improper', ['i', 'j', 'k', 'l', 'type'])
        self.Dihedral = namedtuple('Dihedral', ['i', 'j', 'k', 'l', 'type'])


