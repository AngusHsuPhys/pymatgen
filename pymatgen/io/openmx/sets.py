import os
from pymatgen.io.openmx.inputs import System, Species, Scf, MD
from pymatgen.io.vasp.inputs import Structure
from monty.serialization import loadfn, dumpfn

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))



class ScfInputSet:
    def __init__(self, structure=None, **kwargs):
        self.structure = structure
        self.CONFIG = loadfn(os.path.join(MODULE_DIR, "ScfInputSet.yaml"))

        self.system()
        self.species()
        self.scf()
        self.md()

        self.set = self.set(**kwargs)


    def system(self):
        self.system = System(**self.CONFIG["system"])

    def species(self): 
        unique_elements = list(set(site.specie.symbol for site in self.structure))
        species_config = self.CONFIG["species"]["vpss_and_options"]
        species_list = [species_config[element] for element in unique_elements]
        merged_species = dict(item for species_dict in species_list for item in species_dict.items())
        self.species = Species.get_species_from_vps_and_option(merged_species)

    def scf(self):
        self.scf = Scf.get_scf_with_pmg_kgrid(self.structure, **self.CONFIG["scf"])

    def md(self):
        self.md = MD(**self.CONFIG["md"])

    def set(self, **kwargs):
        input = {}
        for obj in [self.system, self.species, self.scf, self.md]:
            input.update(obj.template)

        input.update(kwargs)
        return input

if __name__ == "__main__":
    structure = Structure.from_file("POSCAR")

    scf_input = ScfInputSet(structure, system_currentdirectory=".", level_of_stdout=2)
    print(scf_input.set)
    


