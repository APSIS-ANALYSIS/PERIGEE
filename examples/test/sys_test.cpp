#include "Vec_Tools.hpp"

class a
{
  public:
    a();
    
    ~a();

    void print_info() const;

    const std::vector<double>& get() const;

  private:
    std::vector<double> av;
};

int main( int argc, char * argv[] )
{
  a * example = new a();
  example -> print_info();
  
  std::cout<<example->get()[0]<<'\n';
  std::cout<<example->get()[1]<<'\n';
  std::cout<<example->get()[2]<<'\n';

  delete example;
  return EXIT_SUCCESS;
}

a::a()
{
  av.clear();
  av.push_back(3.14);
  av.push_back(3.14);
  av.push_back(0.14);
}

a::~a()
{
  VEC_T::clean(av);
}

void a::print_info() const
{
  for(unsigned int ii=0; ii<av.size(); ++ii)
    std::cout<<av[ii]<<'\n';
}

const std::vector<double>& a::get() const
{
  return av;
}

// EOF
