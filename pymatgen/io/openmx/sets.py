import os
from pymatgen.io.openmx.inputs import System, Species, Scf, MD
from pymatgen.io.vasp.inputs import Structure
from monty.serialization import loadfn

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


class ScfInputSet:
    CONFIG = loadfn(os.path.join(MODULE_DIR, "ScfInputSet.yaml"))

    @classmethod
    def write_input(cls, system_name, structure):
        # Initialize system with configuration values
        system_config = cls.CONFIG["system"]
        system = System(
            system_name=system_name,
            system_current_dir=system_config["system_current_dir"],
            level_of_stdout=system_config["level_of_stdout"],
            level_of_fileout=system_config["level_of_fileout"],
        )

        # Extract unique elements from structure
        unique_elements = list(set(site.specie.symbol for site in structure))

        # Create species from configuration
        species_config = cls.CONFIG["species"]["vpss_and_options"]
        species_list = [species_config[element] for element in unique_elements]
        merged_species = dict(item for species_dict in species_list for item in species_dict.items())
        species = Species.get_species_from_vps_and_option(merged_species)


        # Create scf from configuration
        scf = Scf.get_scf_with_pmg_kgrid(structure=structure, **cls.CONFIG["scf"])

        # Create md from configuration
        md = MD(**cls.CONFIG["md"])

        # Concatenate input strings
        input = {}
        for obj in [system, species, scf, md]:
            input.update(obj.template)

        print(input)
        return input
    


if __name__ == "__main__":
    scf_input_set = ScfInputSet()
    structure = Structure.from_file("POSCAR")
    scf_input_set.write_input("GaAs", structure)


