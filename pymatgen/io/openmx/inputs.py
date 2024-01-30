import yaml


class BaseFileWriter:
    def __init__(self, template):
        self.template = template

    def get_string(self):
        # get file_string by filling in template with values from self object
        file_string = self.template.format(**self.__dict__)
        # make sure the lines are aligned
        file_string = "\n".join([line.strip() for line in file_string.split("\n")])
        # now remove all empty lines
        file_string = "\n".join([line for line in file_string.split("\n") if line != ""])
        return file_string

class System(BaseFileWriter):
    def __init__(self, system_name, system_current_dir=".", level_of_stdout=1, level_of_fileout=1):
        self.system_current_dir = system_current_dir
        self.system_name = system_name
        self.level_of_stdout = level_of_stdout
        self.level_of_fileout = level_of_fileout

        template = """\
        #
        # File Name
        #
        System.CurrrentDirectory         {system_current_dir} # default=./
        System.Name                      {system_name}
        level.of.stdout                   {level_of_stdout} # default=1 (1-3)
        level.of.fileout                  {level_of_fileout} # default=1 (0-2)
        """
        super().__init__(template)

class Species(BaseFileWriter):
    def __init__(self, species_number, species_definition):
        self.species_number = species_number
        self.species_definition = species_definition

        template = """\
        #
        # Definition of Atomic Species
        #
        Species.Number       {species_number}
        <Definition.of.Atomic.Species
        {species_definition}
        Definition.of.Atomic.Species>
        """
        super().__init__(template)

    @classmethod
    def get_species_from_vps_and_option(cls, vpss_and_options):
        with open('potential_table.yaml', 'r') as file:
            data = yaml.safe_load(file)

        # order the vpss_and_options by the alphabetical order of the element
        vpss_and_options = dict(sorted(vpss_and_options.items(), key=lambda item: item[0].split('_')[0]))
        
        output = ""
        for vps, option in vpss_and_options.items():
            # if option is not in the list of ["Quick", "Standard", "Precise"], raise ValueError
            if option not in ["Quick", "Standard", "Precise"]:
                raise ValueError(f"Option {option} not in list of options")
            
            # Find the dictionary with the matching VPS
            for d in data:
                if d['VPS'] == vps:
                    element = vps.split('_')[0]
                    pao_string = f"{element} {d[option]} {vps}"
                    output += pao_string + "\n"
                    break
            else:
                # If we get here, we didn't find a match
                raise ValueError(f"VPS {vps} not found in PAO-table.yaml")

        species_number = len(vpss_and_options)
        # call __init__ to get the template
        return Species(species_number=species_number, species_definition=output)




class Atoms(BaseFileWriter):
    def __init__(
        self, 
        atoms_number, 
        atoms_species_and_coordinates_unit, 
        atoms_species_and_coordinates, 
        atoms_unit_vectors_unit, 
        atoms_unit_vectors
    ):
        self.atoms_number = atoms_number
        self.atoms_species_and_coordinates_unit = atoms_species_and_coordinates_unit # FRAC # Ang|AU
        self.atoms_species_and_coordinates = atoms_species_and_coordinates
        self.atoms_unit_vectors_unit = atoms_unit_vectors_unit # Ang|AU
        self.atoms_unit_vectors = atoms_unit_vectors

        template = """\
        #
        # Atoms
        #
        Atoms.Number         {atoms_number}
        Atoms.SpeciesAndCoordinates.Unit   {atoms_species_and_coordinates_unit} # Ang|AU
        <Atoms.SpeciesAndCoordinates           
        {atoms_species_and_coordinates}
        Atoms.SpeciesAndCoordinates>
        Atoms.UnitVectors.Unit             {atoms_unit_vectors_unit}  # Ang|AU
        <Atoms.UnitVectors                     
        {atoms_unit_vectors}
        Atoms.UnitVectors>
        """
        super().__init__(template)

class Scf(BaseFileWriter):
    def __init__(
        self, 
        xc_type, 
        spin_polarization, 
        eigenvalue_solver,
        kgrid, 
        generation_kpoint, 
        mixing_type, 
        spin_orbit_coupling="off",
        electronic_temperature=300,
        energy_cutoff=150,
        max_iter=40,
        init_mixing_weight=0.30,
        min_mixing_weight=0.001,
        max_mixing_weight=0.40,
        mixing_history=5,
        start_pulay=6,
        every_pulay=6,
        criterion=1.0e-6,
        lapack_dste="dstevx"
    ):
        self.xc_type = xc_type
        self.spin_polarization = spin_polarization
        self.spin_orbit_coupling = spin_orbit_coupling
        self.electronic_temperature = electronic_temperature
        self.energy_cutoff = energy_cutoff
        self.max_iter = max_iter
        self.eigenvalue_solver = eigenvalue_solver
        self.kgrid = kgrid
        self.generation_kpoint = generation_kpoint
        self.mixing_type = mixing_type
        self.init_mixing_weight = init_mixing_weight
        self.min_mixing_weight = min_mixing_weight
        self.max_mixing_weight = max_mixing_weight
        self.mixing_history = mixing_history
        self.start_pulay = start_pulay
        self.every_pulay = every_pulay
        self.criterion = criterion
        self.lapack_dste = lapack_dste

        template = """\
        #
        # SCF or Electronic System
        #
        scf.XcType                 {xc_type} # LDA|LSDA-CA|LSDA-PW|GGA-PBE
        scf.SpinPolarization       {spin_polarization} # On|Off|NC
        scf.SpinOrbit.Coupling     {spin_orbit_coupling} # On|Off, default=off
        scf.ElectronicTemperature  {electronic_temperature} # default=300 (K)
        scf.energycutoff           {energy_cutoff}  # default=150 (Ry)
        scf.maxIter                {max_iter}       # default=40
        scf.EigenvalueSolver       {eigenvalue_solver} # DC|GDC|Cluster|Band
        scf.Kgrid                  {kgrid}         # means n1 x n2 x n3
        scf.Generation.Kpoint      {generation_kpoint} # regular|MP
        scf.Mixing.Type            {mixing_type}   # Simple|Rmm-Diis|Gr-Pulay|Kerker|Rmm-Diisk
        scf.Init.Mixing.Weight     {init_mixing_weight} # default=0.30 
        scf.Min.Mixing.Weight      {min_mixing_weight} # default=0.001 
        scf.Max.Mixing.Weight      {max_mixing_weight} # default=0.40 
        scf.Mixing.History         {mixing_history} # default=5
        scf.Mixing.StartPulay      {start_pulay} # default=6
        scf.Mixing.EveryPulay      {every_pulay} # default=6
        scf.criterion              {criterion}     # default=1.0e-6 (Hartree) 
        scf.lapack.dste            {lapack_dste} # dstevx|dstedc|dstegr,default=dstevx
        """
        super().__init__(template)


class MD(BaseFileWriter):
    def __init__(self, md_type, md_max_iter=1, md_time_step=0.5, md_opt_criterion=1.0e-4):
        self.md_type = md_type 
        self.md_max_iter = md_max_iter
        self.md_time_step = md_time_step
        self.md_opt_criterion = md_opt_criterion

        template = """\
        #
        # MD or Geometry Optimization
        #
        MD.Type                      {md_type}       # Nomd|Constant_Energy_MD|Opt
        MD.maxIter                    {md_max_iter}        # default=1
        MD.TimeStep                   {md_time_step}        # default=0.5 (fs)
        MD.Opt.criterion              {md_opt_criterion}        # default=1.0e-4 (Hartree/bohr)
        """
        super().__init__(template)





if __name__ == "__main__":
    filename = System(system_current_dir=".", system_name="test", level_of_stdout=1, level_of_fileout=1)
    print(filename.get_string())

    # Test get_species_def_from_vps_and_option
    inputs = {
        "Ni_PBE19H": "Quick",
        "Cu_PBE19H": "Quick",
        "Mn_PBE19": "Quick",
    }

    print(Species.get_species_from_vps_and_option(inputs).get_string())


    atoms_number = 2
    atoms_species_and_coordinates_unit = "FRAC"
    atoms_species_and_coordinates = """\
 1  Ga  0.0000  0.0000  0.0000   6.5 6.5  0.0 0.0 0.0 0.0 1
 2  As  0.2500  0.2500  0.2500   7.5 7.5  0.0 0.0 0.0 0.0 1
 """
    atoms_unit_vectors_unit = "Au"
    atoms_unit_vectors = """\
 5.367  0.000  5.367
 0.000  5.367  5.367
 5.367  5.367  0.000
 """
    filename = Atoms(atoms_number=atoms_number, atoms_species_and_coordinates_unit=atoms_species_and_coordinates_unit, atoms_species_and_coordinates=atoms_species_and_coordinates, atoms_unit_vectors_unit=atoms_unit_vectors_unit, atoms_unit_vectors=atoms_unit_vectors)
    print(filename.get_string())

    xc_type = "LSDA-CA"
    spin_polarization = "off"
    eigenvalue_solver = "band"
    kgrid = "7 7 7"
    generation_kpoint = "regular"
    mixing_type = "Rmm-Diis"
    
    filename = Scf(
        xc_type=xc_type, 
        spin_polarization=spin_polarization,
        eigenvalue_solver=eigenvalue_solver,
        kgrid=kgrid, 
        generation_kpoint=generation_kpoint, 
        mixing_type=mixing_type, 
    )
    print(filename.get_string())

    filename = MD(md_type="Opt", md_max_iter=1, md_time_step=0.5, md_opt_criterion=1.0e-4)
    print(filename.get_string())


