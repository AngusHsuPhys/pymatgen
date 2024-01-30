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
    def __init__(self, system_current_dir, system_name, level_of_stdout, level_of_fileout):
        self.system_current_dir = system_current_dir if system_current_dir else "."
        self.system_name = system_name
        self.level_of_stdout = level_of_stdout if level_of_stdout else 1
        self.level_of_fileout = level_of_fileout if level_of_fileout else 1

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
        spin_orbit_coupling, 
        electronic_temperature, 
        energy_cutoff, 
        max_iter, 
        eigenvalue_solver,
        kgrid, 
        generation_kpoint, 
        mixing_type, 
        init_mixing_weight, 
        min_mixing_weight, 
        max_mixing_weight, 
        mixing_history, 
        start_pulay, 
        every_pulay, 
        criterion, 
        lapack_dste
    ):
        self.xc_type = xc_type
        self.spin_polarization = spin_polarization
        self.spin_orbit_coupling = spin_orbit_coupling
        self.electronic_temperature = electronic_temperature if electronic_temperature else 300
        self.energy_cutoff = energy_cutoff if energy_cutoff else 150
        self.max_iter = max_iter if max_iter else 40
        self.eigenvalue_solver = eigenvalue_solver
        self.kgrid = kgrid
        self.generation_kpoint = generation_kpoint
        self.mixing_type = mixing_type
        self.init_mixing_weight = init_mixing_weight if init_mixing_weight else 0.30
        self.min_mixing_weight = min_mixing_weight if min_mixing_weight else 0.001
        self.max_mixing_weight = max_mixing_weight if max_mixing_weight else 0.40
        self.mixing_history = mixing_history if mixing_history else 5
        self.start_pulay = start_pulay if start_pulay else 6
        self.every_pulay = every_pulay if every_pulay else 6
        self.criterion = criterion if criterion else 1.0e-6
        self.lapack_dste = lapack_dste if lapack_dste else "dstevx"

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
    def __init__(self, md_type, md_max_iter, md_time_step, md_opt_criterion):
        self.md_type = md_type 
        self.md_max_iter = md_max_iter if md_max_iter else 1
        self.md_time_step = md_time_step if md_time_step else 0.5
        self.md_opt_criterion = md_opt_criterion if md_opt_criterion else 1.0e-4

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

    species_number = 3
    species_definition = """\
 Ga   Ga7.0-s2p2d1   Ga_CA19
 As   As7.0-s2p2d1   As_CA19
proj  As7.0-s1p1d1   As_CA19
        """
    filename = Species(species_number=species_number, species_definition=species_definition)
    print(filename.get_string())


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
