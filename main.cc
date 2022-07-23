#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <array>
#include <unordered_map>
#include <sstream>
#include <dftd4.h>

int number_of_element(std::string element) {
  std::unordered_map<std::string, size_t> symbols = {
    {"H", 1},
    {"C", 6},
    {"O", 8},
    {"P", 31},
    {"Pd", 46}
  };
  return symbols.at(element);
}

struct atom {

  atom() : number(0), positions() {}

  /**
   * atomic number
   */
  int number;

  /**
   * atomic symbol
   */
  std::string symbol;

  /**
   * atomic positions
   */
  std::array<double, 3> positions;

  // Extractor and inserter
  friend std::istream& operator>>(std::istream& is, atom& a) {
    is >> a.symbol >> a.positions[0] >> a.positions[1] >> a.positions[2]; 
    a.number = number_of_element(a.symbol);
    return is;
  }

  friend std::ostream& operator<<(std::ostream& os, const atom& a) { 
    os << "[ " << a.symbol << ", " << a.number << ", " << a.positions[0] << ", " << a.positions[1] << ", " << a.positions[2] << " ]"; 
    return os;
  }

};

struct molecule {

  /**
   * constructor
   */
  molecule() : n_atoms(0), numbers(), positions(), charge(), format('A') {}

  /**
   * convert to Bohr
   */
  void to_bohr() {
    if (format != 'B') {
      for (auto& e : positions) {
        e *= 1.88973;
      }
      format = 'B';
    }
  }

  void to_angstrom() {
    if (format != 'A') {
      for (auto& e : positions) {
        e *= 0.529177;
      }
      format = 'A';
    }
  }
  
  /**
   * number of atoms
   */
  int n_atoms;

  /**
   * atomic numbers
   */
  std::vector<int> numbers;

  /**
   * atomic coordinates in Bohr (one-dimensional array)
   */
  std::vector<double> positions;

  /**
   * charge
   */
  double charge;

  /**
   * Angstrom or Bohr
   */
  char format;

  void clear() {
    n_atoms = 0;
    charge = 0;
    numbers.clear();
    positions.clear();
    format = 'A';
  }

  friend std::istream& operator>>(std::istream& ifs, molecule& mol) {
    // clear old data
    mol.clear();

    // read all the lines from source stream
    size_t i = 0;
    for (std::string line{}; std::getline(ifs, line); ++i) {
      if (i > 1) {
        atom a;
        std::istringstream(line) >> a;
        mol.numbers.push_back(a.number);
        mol.positions.push_back(a.positions[0]);
        mol.positions.push_back(a.positions[1]);
        mol.positions.push_back(a.positions[2]);
      }
      mol.n_atoms = mol.numbers.size();
    }

    return ifs;
  }

  friend std::ostream& operator<<(std::ostream& os, const molecule& mol) {
    os << "[ " << std::endl;
    os << "# atoms: " << mol.n_atoms << ", charge: " << mol.charge << ", ";
    for (size_t i = 0; i < mol.n_atoms; ++i) {
      os << "element: " << mol.numbers[i] << ", positions: " << mol.positions[3 * i] << ", " << mol.positions[3 * i + 1] << ", " << mol.positions[3 * i + 2] << std::endl;
    }
    os << " ]";
    return os;
  }

};

int main(int argc, char** argv) {
  std::string usage = "dftb4-example filename.xyz charge";
  if (argc < 3) {
    std::cerr << usage << std::endl;
    return 1;
  }

  molecule mol;

  if (std::ifstream ifs{argv[1]}; ifs) {
    ifs >> mol;
    mol.charge = std::stod(argv[2]);
  }
  else {
    std::cerr << "Could not read in file." << std::endl;
    std::cerr << strerror(errno) << std::endl;
  }

  // DFT-D4 uses Bohr
  mol.to_bohr();

  dftd4_error error = dftd4_new_error();

  const int natoms = mol.n_atoms;
  const int* numbers = mol.numbers.data();
  const double* positions = mol.positions.data();
  const double* charge = &mol.charge;
  dftd4_structure structure = dftd4_new_structure(error, natoms, numbers, positions, charge, NULL, NULL);
  dftd4_model model = dftd4_new_d4_model(error, structure);

  std::vector<double> c_n;
  std::vector<double> charges;
  std::vector<double> c6;
  std::vector<double> alphas;
  c_n.resize(natoms);
  charges.resize(natoms);
  c6.resize(natoms * natoms);
  alphas.resize(natoms);

  dftd4_get_properties(error, structure, model, c_n.data(), charges.data(), c6.data(), alphas.data());

  for (size_t i = 0; i < natoms; ++i) {
    std::cout << "[ atom #" << i << ": CN " << c_n[i] << ", Q " << charges[i] << ", C6 " << c6[natoms + i] << ", Alpha " << alphas[i] << "]" << std::endl;
  }

  dftd4_delete_structure(&structure);
  dftd4_delete_model(&model);
  dftd4_delete_error(&error);
}
