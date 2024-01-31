import yaml, os
from pymatgen.io.openmx.inputs import System, Species, Atoms, Scf, MD
from pymatgen.io.vasp.inputs import Structure
from monty.serialization import loadfn

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


class ScfInputSet:
    CONFIG = loadfn(os.path.join(MODULE_DIR, "ScfInputSet.yaml"))

    @classmethod
    def write_input(cls, system_name, structure):
        system = System(
            system_name=system_name,
            system_current_dir=cls.CONFIG["system"]["system_current_dir"],
            level_of_stdout=cls.CONFIG["system"]["level_of_stdout"],
            level_of_fileout=cls.CONFIG["system"]["level_of_fileout"],
        )


        # get element list of structure
        element_list = []
        for site in structure:
            element_list.append(site.specie.symbol)
        
        #remove duplicate element
        element_list = list(set(element_list))

        # create species
        vpss_and_options = list()
        for element in element_list:
            vpss_and_options.append(cls.CONFIG["species"]["vpss_and_options"][element])
        
        vpss_and_options = dict(item for item in vpss_and_options for item in item.items())
        species = Species.get_species_from_vps_and_option(vpss_and_options)

        # create atoms
        atoms = Atoms.get_atoms_from_pmg_structure(
            structure=structure, 
            vpss=vpss_and_options.keys(),
            fractional_coordinates=cls.CONFIG["atoms"]["fractional_coordinates"],
            up_dn_diff=None if cls.CONFIG["atoms"]["up_dn_diff"] == "None" else cls.CONFIG["atoms"]["up_dn_diff"],
        )

        # create scf
        scf_config = cls.CONFIG["scf"]
        scf = Scf.get_scf_with_pmg_kgrid(
            structure=structure,
            **scf_config
        )

        # create md
        md_config = cls.CONFIG["md"]
        md = MD(**md_config)


        # write input
        input_str = ""
        input_str += system.get_string()
        input_str += species.get_string()
        input_str += atoms.get_string()
        input_str += scf.get_string()
        input_str += md.get_string()

        with open("input.dat", "w") as f:
            f.write(input_str)

        return input_str
    


if __name__ == "__main__":
    scf_input_set = ScfInputSet()
    structure = Structure.from_file("POSCAR")
    scf_input_set.write_input("GaAs", structure)


