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
        self.system_current_dir = system_current_dir if system_current_dir != None else "./"
        self.system_name = system_name if system_name != None else "test"
        self.level_of_stdout = level_of_stdout if level_of_stdout != None else 1
        self.level_of_fileout = level_of_fileout if level_of_fileout != None else 1

        template = f"""\
        #
        # File Name
        #
        System.CurrrentDirectory         {self.system_current_dir} # default=./
        System.Name                      {self.system_name}
        level.of.stdout                   {self.level_of_stdout} # default=1 (0-2)
        level.of.fileout                  {self.level_of_fileout} # default=1 (0-2)
        """
        super().__init__(template)


class Atoms(BaseFileWriter):
    def __init__(self, atom_number, atom_coordinate_unit, atom_coordinate, lattice_unit, lattice_vector):
        self.atom_number = atom_number
        self.atom_coordinate_unit = atom_coordinate_unit # FRAC # Ang|AU
        self.atom_coordinate = atom_coordinate
        self.lattice_unit = lattice_unit # Ang|AU
        self.lattice_vector = lattice_vector

        template = """\
        #
        # Atoms
        #
        Species.Number       {atom_number}
        Atoms.SpeciesAndCoordinates.Unit       {atom_coordinate_unit}
        <Atoms.SpeciesAndCoordinates
        {atom_coordinate}
        Atoms.SpeciesAndCoordinates>
        Atoms.UnitVectors.Unit       {lattice_unit}
        <Atoms.UnitVectors
        {lattice_vector}
        Atoms.UnitVectors>
        """
        super().__init__(template)

    

if __name__ == "__main__":
    filename = System(system_current_dir=".", system_name="test", level_of_stdout=1, level_of_fileout=1)
    print(filename.get_string())

    atom_number = 2
    atom_coordinate_unit = "FRAC"
    atom_coordinate = """\
        1  Ga  0.0000  0.0000  0.0000   6.5 6.5  0.0 0.0 0.0 0.0 1
        2  As  0.2500  0.2500  0.2500   7.5 7.5  0.0 0.0 0.0 0.0 1
        """
    lattice_unit = "Ang"
    lattice_vector = """\
    5.367  0.000  5.367
    0.000  5.367  5.367
    5.367  5.367  0.000
        """
    filename = Atoms(atom_number=atom_number, atom_coordinate_unit=atom_coordinate_unit, atom_coordinate=atom_coordinate, lattice_unit=lattice_unit, lattice_vector=lattice_vector)
    print(filename.get_string())
