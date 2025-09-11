import os
from pymatgen.io.openmx.inputs import System, Species, Scf, MD
from pymatgen.io.vasp.inputs import Structure
from monty.serialization import loadfn, dumpfn

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))



class ScfInputSet:
    def __init__(self, structure=None, **kwargs):
        self.structure = structure
        self.CONFIG = loadfn(os.path.join(MODULE_DIR, "ScfInputSet.yaml"))
        self.input_params = kwargs if kwargs else {}
        self.system()
        self.species()
        self.scf()
        self.md()
        # self.optical_conductivity()

    def system(self):
        self.system = System(**self.CONFIG["system"])

    def species(self): 
        unique_elements = list(set(site.specie.symbol for site in self.structure))
        species_config = self.CONFIG["species"]["vpss_and_options"]
        species_list = [species_config[element] for element in unique_elements] # this looks like [{"Ag_PBE19", "Standard"}, {"Ga_PBE19", "Standard"}]
        
        if self.input_params.get("potcar_spec"): # this looks like {"Ag": "Quick", "Ga": "Quick"}
            new_species_list = []
            for element, option in self.input_params["potcar_spec"].items():
                vps = list(species_config[element].keys())[0] # this looks like "Ag_PBE19"
                for specie_dict in species_list: # this looks like {"Ag_PBE19": "Standard"}
                    if vps in specie_dict:
                        specie_dict[vps] = option
                        new_species_list.append(specie_dict)
                    else:
                        new_species_list.append(specie_dict)
            species_list = new_species_list

        merged_species = dict(item for species_dict in species_list for item in species_dict.items())
        self.species = Species.get_species_from_vps_and_option(merged_species)

    def scf(self):
        # if self.input_params.get("kppa"):
        #update the kppa in the scf config
        if self.input_params.get("kppa"):
            self.CONFIG["scf"]["kppa"] = self.input_params["kppa"]   
        self.scf = Scf.get_scf_with_pmg_kgrid(self.structure, **self.CONFIG["scf"])

    def md(self):
        self.md = MD(**self.CONFIG["md"])

    # def optical_conductivity(self):
    #     if self.input_params.get("CDDF") == "off":
    #         self.optical_conductivity = None
    #     elif self.input_params.get("CDDF") == "on":
    #         self.optical_conductivity = CDDF(**self.CONFIG["CDDF"])
    #         if self.input_params.get("kgrid_density"):
    #             # self.CONFIG["CDDF"]["kgrid_density"] = self.input_params["kgrid_density"]  
    #             self.optical_conductivity = CDDF.get_optcond_with_pmg_kgrid(self.structure, self.input_params["kgrid_density"] )


    def as_dict(self):
        input = {}

        for obj in [self.system, self.species, self.scf, self.md]:
            input.update(obj.template)

        # # Include OpticalConductivity settings if available
        # if self.optical_conductivity:
        #     input.update(self.optical_conductivity.template)
        # Remove kppa if it was in input_params
        if self.input_params.get("kppa"):
            del self.input_params["kppa"]
            
        input.update(self.input_params)
        return input


if __name__ == "__main__":
    structure = Structure.from_file("POSCAR")

    scf_input = ScfInputSet(structure, system_currentdirectory=".", level_of_stdout=2)
    print(scf_input.as_dict())
    


