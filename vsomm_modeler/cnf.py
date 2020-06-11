import sys
import numpy as np
import getpass
from datetime import datetime

__author__ = "Yerko Escalona"


class CNF():
    """
    Simple CNF class

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

                   Block name    POSITION
    1 -  5         Integer       resSeq       Residue sequence number
    7 -  11        Residue name  resName      Residue name
    13 - 16        Atom          name         Atom name
    21?- 24        Integer       serial       Atom serial number
    28 - 39        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    43 - 54        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    58 - 69        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.

                   Block name    LATTICESHIFTS
    9 -  10        Integer       x            Orthogonal coordinates for X in Angstroms.
    19 - 20        Integer       y            Orthogonal coordinates for Y in Angstroms.
    29 - 30        Integer       z            Orthogonal coordinates for Z in Angstroms.

                   Block name    VELOCITY
    1 -  5         Integer       resSeq       Residue sequence number
    7 -  11        Residue name  resName      Residue name
    13 - 16        Atom          name         Atom name
    21?- 24        Integer       serial       Atom serial number
    28 - 39        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    43 - 54        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    58 - 69        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.

                   Block name    GENBOX
    1:8            Integer       box          box type
    2:1- 15        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    2:16-30        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    2:31-45        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    3:1- 15        Real(8.3)     x            degrees.
    3:16-30        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    3:31-45        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    4:1- 15        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    4:16-30        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    4:31-45        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    5:1- 15        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    5:16-30        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    5:31-45        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    =============  ============  ===========  =============================================
    """

    def __init__(self, file):
        self.__version__ = "1.0"
        self.blocks = []

        self.resseqs = []
        self.resnames = []
        self.names = []
        self.atoms = []
        self.positions = np.empty(shape=[0, 3], dtype=float)
        self.velocities = np.empty(shape=[0, 3], dtype=float)
        self.latticeshifts = np.empty(shape=[0, 3], dtype=int)
        self.boxsize = np.empty(shape=[0, 3], dtype=float)

        self.parse_file(file)

    def item_generator(self, src):
        line = src.readline()
        block = []
        while line:
            while line and not line.startswith("END"):
                if line.startswith("#"):
                    line = src.readline()
                    continue
                block.append(line)
                line = src.readline()
            yield block
            block = []
            line = src.readline()

    def get_blocks(self, src):
        return {item[0].strip(): item[1:] for item in self.item_generator(src)}

    def parse_file(self, cnf_filename):
        with open(cnf_filename) as file:
            self.blocks = self.get_blocks(file)
            self.parse_blocks(self.blocks)

    def parse_blocks(self, blocks):

        try:
            block_lines = blocks["TIMESTEP"]
        except:
            print("WARNING: No TIMESTEP block")
            blocks["TIMESTEP"] = "              0    0.0        \n"

        try:
            block_lines = blocks["POSITION"]
            for line in block_lines:
                self.resseqs.append(line[1:6].strip())
                self.resnames.append(line[6:12].strip())
                self.names.append(line[12:17].strip())
                self.atoms.append(line[20:25].strip())
                self.positions = np.vstack(
                    [self.positions, [float(line[24:39]),
                                      float(line[39:54]),
                                      float(line[54:69])]])
        except:
            print("No POSITION block")
            sys.exit(1)

        try:
            block_lines = blocks["LATTICESHIFTS"]
            for line in block_lines:
                self.latticeshifts = np.vstack(
                    [self.latticeshifts, [int(line[0:10]),
                                          int(line[10:20]),
                                          int(line[20:30])]])
        except:
            print("WARNING: No LATTICESHIFTS block")
            for i in range(len(self.atoms)):
                self.latticeshifts = np.vstack([self.latticeshifts, [0, 0, 0]])

        try:
            block_lines = blocks["VELOCITY"]
            for line in block_lines:
                self.velocities = np.vstack(
                    [self.velocities, [float(line[24:39]),
                                       float(line[39:54]),
                                       float(line[54:69])]])
        except:
            print("WARNING: No VELOCITY block")
            for i in range(len(self.atoms)):
                self.velocities = np.vstack([self.velocities, [0, 0, 0]])

        try:
            block_lines = blocks["GENBOX"]
            self.ntb = int(block_lines[0])
            self.boxsize = np.array([float(block_lines[1][0:15]),
                                     float(block_lines[1][16:30]),
                                     float(block_lines[1][31:45])])
            self.angle = np.array([float(block_lines[2][0:15]),
                                   float(block_lines[2][16:30]),
                                   float(block_lines[2][31:45])])
            self.euler = np.array([float(block_lines[3][0:15]),
                                   float(block_lines[3][16:30]),
                                   float(block_lines[3][31:45])])
            self.cartesian = np.array([float(block_lines[4][0:15]),
                                       float(block_lines[4][16:30]),
                                       float(block_lines[4][31:45])])
        except:
            print("No GENBOX block")
            sys.exit(1)

    def output(self, output_filename="out.cnf"):
        
        with open(output_filename, 'w') as file:

            file.write("TITLE\n")
            file.write("  CNF file automatically generated by using CNF class v" + self.__version__ + "\n")
            file.write("  by " + getpass.getuser() + " at " + str(datetime.now()) + "\n")
            file.write("END\n")

            file.write("TIMESTEP\n")
            output.extend(self.blocks["TIMESTEP"])
            file.write("END\n")

            file.write("POSITION\n")
            file.write("# first 24 chars ignored\n")
            position_format = ("{resSeq:>5} {resName:<5s} {name:<5s} {atom:>6}"
                               "{position[0]: 15.9f}{position[1]: 15.9f}{position[2]: 15.9f}\n")
            for i in range(len(self.atoms)):
                vals = {}
                vals['resSeq'] = self.resseqs[i]
                vals['resName'] = self.resnames[i]
                vals['name'] = self.names[i]
                vals['atom'] = self.atoms[i]
                vals['position'] = self.positions[i]
                file.write(position_format.format(**vals))
            file.write("END\n")

            file.write("LATTICESHIFTS\n")
            latticeshift_format = "{latticeshift[0]:10}{latticeshift[1]:10}{latticeshift[2]:10}\n"
            for i in range(len(self.atoms)):
                vals = {}
                vals['latticeshift'] = self.latticeshifts[i]
                file.write(latticeshift_format.format(**vals))
            file.write("END\n")

            file.write("VELOCITY\n")
            file.write("# first 24 chars ignored\n")
            velocity_format = ("{resSeq:>5} {resName:<5s} {name:<5s} {atom:>6}"
                               "{velocity[0]: 15.9f}{velocity[1]: 15.9f}{velocity[2]: 15.9f}\n")

            for i in range(len(self.atoms)):
                vals = {}
                vals['resSeq'] = self.resseqs[i]
                vals['resName'] = self.resnames[i]
                vals['name'] = self.names[i]
                vals['atom'] = self.atoms[i]
                vals['velocity'] = self.velocities[i]
                file.write(velocity_format.format(**vals))
            file.write("END\n")

            file.write("GENBOX\n")
            genbox_format = ("{ntb:>5}\n"
                             "{boxsize[0]: 15.9f}{boxsize[1]: 15.9f}{boxsize[2]: 15.9f}\n"
                             "{angle[0]: 15.9f}{angle[1]: 15.9f}{angle[2]: 15.9f}\n"
                             "{euler[0]: 15.9f}{euler[1]: 15.9f}{euler[2]: 15.9f}\n"
                             "{cartesian[0]: 15.9f}{cartesian[1]: 15.9f}{cartesian[2]: 15.9f}\n")
            vals = {}
            vals['ntb'] = self.ntb
            vals['boxsize'] = self.boxsize
            vals['angle'] = self.angle
            vals['euler'] = self.euler
            vals['cartesian'] = self.cartesian
            file.write(genbox_format.format(**vals))
            file.write("END\n")

    def output_gro(self, output_filename="out.gro"):
        
        with open(output_filename, 'w') as file:
            file.write("GRO file automatically generated by using CNF class v" + self.__version__ + "\n")
            file.write("%5d\n" % len(self.atoms))
        
            position_format = ("{resSeq:>5}{resName:<5s}{name:>5s}{atom:>5}"
                               "{position[0]: 8.3f}{position[1]: 8.3f}{position[2]: 8.3f}\n")
            for i in range(len(self.atoms)):
                vals = {}
                vals['resSeq'] = self.resseqs[i]
                vals['resName'] = self.resnames[i]
                vals['name'] = self.names[i]
                vals['atom'] = self.atoms[i]
                vals['position'] = self.positions[i]
                file.write(position_format.format(**vals))
            
            boxsize_format = "{boxsize[0]: 10.5f}{boxsize[1]: 10.5f}{boxsize[2]: 10.5f}"
            file.write(boxsize_format.format(boxsize=self.boxsize))
