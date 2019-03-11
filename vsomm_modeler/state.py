
import numpy as np

global CARBON_MASS
global HYDROGEN_MASS
global OXYGEN_MASS
global NITROGEN_MASS
global SULPHUR_MASS

CARBON_MASS = 12.011
HYDROGEN_MASS = 1.008
OXYGEN_MASS = 15.999
NITROGEN_MASS = 14.007
SULPHUR_MASS = 32.06

global PKA_COOH

# TODO: quitar esta variable y ocupar deprot_prot_fraction
PKA_COOH = 4.388



class State_SSNMR(object):
    """
    State class

    """

    def __init__(self):
        """ Initialization of the state of a system with multiple building
        blocks. """

        # properties obtained from the building blocks
        self.properties = ["carbon_count", "hydrogen_count", "oxygen_count",
                           "nitrogen_count", "sulphur_count",
                           "carbonyl_c_count", "carboxyl_c_count",
                           "o_aryl_count", "aryl_count", "di_o_alkyl_count",
                           "o_alkyl_count", "methoxyl_count", "alkyl_c_count",
                           "protonated", "deprotonated", "atom_count", "charge"]

        # initialize them with 0
        for prop in self.properties:
            setattr(self, prop, 0)

        self.mass = 0.0

        # properties to compare
        self.fractions = ["carbon_fraction", "nitrogen_fraction",
                          "hydrogen_fraction", "oxygen_fraction",
                          "sulphur_fraction", "carbonyl_c_fraction",
                          "carboxyl_c_fraction", "o_aryl_fraction",
                          "aryl_fraction", "di_o_alkyl_fraction",
                          "o_alkyl_fraction", "methoxyl_fraction",
                          "alkyl_c_fraction", "deprot_prot_fraction"]

        # initialize them with 0
        for frac in self.fractions:
            setattr(self, frac, 0.0)

        # reinitialize deprot_prot_fraction
        self.deprot_prot_fraction = -1.0

        # initialize logK1
        self.logK1 = 0.0

        # list of building blocks
        self.list_building_blocks = []

    def append(self, building_block):
        """Append the properties of a new building block to the molecule."""

        for prop in self.properties:
            value = getattr(self, prop)
            value += building_block[prop]
            setattr(self, prop, value)

        self.list_building_blocks.append(building_block['name'])

    def remove(self, building_block):

        for prop in self.properties:
            value = getattr(self, prop)
            value -= building_block[prop]
            setattr(self, prop, value)

    #def set_list(self, list_building_blocks):
    #    """Update the list of properties of the molecule."""

    #    for building_block in list_building_blocks:
    #        for prop in self.properties:
    #            value = getattr(self, prop)
    #            value += building_block[prop]
    #            setattr(self, prop, value)
    #        self.list_building_blocks.append(building_block['name'])

    #    # self.set_status()

    def modify(self, current_building_block, position, new_building_block):
        if self.list_building_blocks[position] is current_building_block['name']:
            for prop in self.properties:
                value = getattr(self, prop)
                # remove
                value -= current_building_block[prop]
                # add
                value += new_building_block[prop]
                setattr(self, prop, value)

            self.list_building_blocks[position] = new_building_block['name']
        else:
            print("[!] Error: No %s in position %s" % (current_building_block, position))


    def set_status(self):
        """Define the properties of the molecule."""

        # define mass
        self.mass = self.carbon_count * CARBON_MASS + \
            self.hydrogen_count * HYDROGEN_MASS + \
            self.oxygen_count * OXYGEN_MASS + \
            self.nitrogen_count * NITROGEN_MASS + \
            self.sulphur_count * SULPHUR_MASS

        # define carbon and nitrogen fractions
        self.carbon_fraction = self.carbon_count * CARBON_MASS / self.mass
        self.hydrogen_fraction = self.hydrogen_count * HYDROGEN_MASS / self.mass
        self.oxygen_fraction = self.oxygen_count * OXYGEN_MASS / self.mass
        self.nitrogen_fraction = self.nitrogen_count * NITROGEN_MASS / self.mass
        self.sulphur_fraction = self.sulphur_count * SULPHUR_MASS / self.mass

        # define carbohydrate composition
        prop_to_compare = ["carbonyl_c_count", "carboxyl_c_count", "o_aryl_count",
                           "aryl_count", "di_o_alkyl_count", "o_alkyl_count",
                           "methoxyl_count", "alkyl_c_count"]
        frac_to_compare = ["carbonyl_c_fraction", "carboxyl_c_fraction",
                           "o_aryl_fraction", "aryl_fraction",
                           "di_o_alkyl_fraction", "o_alkyl_fraction",
                           "methoxyl_fraction", "alkyl_c_fraction"]

        for prop, frac in zip(prop_to_compare, frac_to_compare):
            setattr(self, frac, getattr(self, prop) / float(self.carbon_count))

        if self.protonated != 0 and self.deprotonated != 0:
            self.deprot_prot_fraction = self.deprotonated / float(self.protonated)
        else:
            self.deprot_prot_fraction = -1

    def status(self):
        """Print status."""

        for frac in self.fractions:
            print("{:<20} {:8.2f}".format(frac, np.around(getattr(self, frac), decimals=3)))


class State_IHSS(object):
    """
    State class

    """

    #__slots__ = "properties", "fractions", "list_building_blocks", "carbon_count", "hydrogen_count", "oxygen_count", "nitrogen_count", "sulphur_count", "carbonyl_c_count", "carboxyl_c_count", "o_aryl_count", "aryl_count", "di_o_alkyl_count", "o_alkyl_count", "methoxyl_count", "alkyl_c_count", "protonated", "deprotonated", "atom_count", "charge", "carbon_fraction", "nitrogen_fraction", "hydrogen_fraction", "oxygen_fraction", "sulphur_fraction", "carbonyl_fraction", "carboxyl_fraction", "aromatic_fraction", "acetal_fraction", "heteroaliphatic_fraction", "aliphatic_fraction", "pH", "mass"

    def __init__(self):
        """ Initialization of the state of a system with multiple building blocks."""

        # properties obtained from the building blocks
        self.properties = ["carbon_count", "hydrogen_count", "oxygen_count",
                           "nitrogen_count", "sulphur_count",
                           "carbonyl_c_count", "carboxyl_c_count",
                           "o_aryl_count", "aryl_count", "di_o_alkyl_count",
                           "o_alkyl_count", "methoxyl_count",
                           "alkyl_c_count", "protonated", "deprotonated",
                           "atom_count", "charge"]

        # initialize them with 0
        for prop in self.properties:
            setattr(self, prop, 0)

        self.mass = 0.0

        # properties to compare
        self.fractions = ["carbon_fraction", "nitrogen_fraction",
                          "hydrogen_fraction", "oxygen_fraction",
                          "sulphur_fraction", "carbonyl_fraction",
                          "carboxyl_fraction", "aromatic_fraction",
                          "acetal_fraction", "heteroaliphatic_fraction",
                          "aliphatic_fraction", "deprot_prot_fraction"]

        # initialize them with 0
        for frac in self.fractions:
            setattr(self, frac, 0.0)

        # list of building blocks
        self.list_building_blocks = []

    def append(self, building_block):
        """Append the properties of a new building block to the molecule."""

        for prop in self.properties:
            value = getattr(self, prop)
            value += building_block[prop]
            setattr(self, prop, value)

        self.list_building_blocks.append(building_block['name'])

    def set_list(self, list_building_blocks):
        """Update the list of properties of the molecule."""

        for building_block in list_building_blocks:
            for prop in self.properties:
                value = getattr(self, prop)
                value += building_block[prop]
                setattr(self, prop, value)
            self.list_building_blocks.append(building_block['name'])

    def modify(self, current_building_block, position, new_building_block):
        if self.list_building_blocks[position] is current_building_block['name']:
            for prop in self.properties:
                value = getattr(self, prop)
                # remove
                value -= current_building_block[prop]
                # add
                value += new_building_block[prop]
                setattr(self, prop, value)

            self.list_building_blocks[position] = new_building_block['name']
        else:
            print("[!] Error: No %s in position %s" % (current_building_block, position))

    def set_status(self):
        """Define the properties of the molecule."""

        # define mass
        self.mass = self.carbon_count * CARBON_MASS + \
            self.hydrogen_count * HYDROGEN_MASS + \
            self.oxygen_count * OXYGEN_MASS + \
            self.nitrogen_count * NITROGEN_MASS + \
            self.sulphur_count * SULPHUR_MASS

        # define carbon and nitrogen fractions
        self.carbon_fraction = self.carbon_count * CARBON_MASS / self.mass
        self.hydrogen_fraction = self.hydrogen_count * HYDROGEN_MASS / self.mass
        self.oxygen_fraction = self.oxygen_count * OXYGEN_MASS / self.mass
        self.nitrogen_fraction = self.nitrogen_count * NITROGEN_MASS / self.mass
        self.sulphur_fraction = self.sulphur_count * SULPHUR_MASS / self.mass

        self.carbonyl_count = self.carbonyl_c_count
        self.carboxyl_count = self.carboxyl_c_count
        self.aromatic_count = (self.o_aryl_count + self.aryl_count)
        self.acetal_count = self.di_o_alkyl_count
        self.heteroaliphatic_count = self.o_alkyl_count
        self.aliphatic_count = self.methoxyl_count + self.alkyl_c_count

        # define carbohydrate composition
        self.carbonyl_fraction = self.carbonyl_c_count / self.carbon_count
        self.carboxyl_fraction = self.carboxyl_c_count / self.carbon_count
        self.aromatic_fraction = (self.o_aryl_count + self.aryl_count) / self.carbon_count
        self.acetal_fraction = self.di_o_alkyl_count / self.carbon_count
        self.heteroaliphatic_fraction = self.o_alkyl_count / self.carbon_count
        self.aliphatic_fraction = (self.methoxyl_count + self.alkyl_c_count) / self.carbon_count

        # define pH
        if self.protonated != 0 and self.deprotonated != 0:
            self.deprot_prot_fraction = self.deprotonated / float(self.protonated)
        else:
            self.deprot_prot_fraction = -1

    def status(self):
        """Print status."""

        for frac in self.fractions:
            print("{:<20} {:8.2f}".format(frac, np.around(getattr(self, frac), decimals=3)))
