#include <aspect/utilities.h>

#include <iostream>

void print_string_vector(std::vector<std::string> data)
{
  for (std::vector<std::string>::iterator it = data.begin(); it != data.end(); ++it)
    std::cout << "\t" << *it << std::endl;
}

int f()
{
  using namespace aspect::Utilities;

  std::vector<std::string> base_declarations;
  base_declarations.push_back("vector(  vec  )");
  base_declarations.push_back("tensor(  ten )");

  std::vector<std::string> var_decl_2 = expand_dimensional_variable_names<2> (base_declarations);
  std::vector<std::string> var_decl_3 = expand_dimensional_variable_names<3> (base_declarations);

  std::cout << "Dimension 2:" << std::endl;
  print_string_vector(var_decl_2);

  std::cout << "Dimension 3:" << std::endl;
  print_string_vector(var_decl_3);

  exit(0);

  // dummy return for compilation warnings

  return 0;
}

int ignored_variable = f();
