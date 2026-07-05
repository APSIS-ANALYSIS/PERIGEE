#ifndef HIMODEL_POWER_LAW_HPP
#define HIMODEL_POWER_LAW_HPP
// ============================================================================
// HIModel_PowerLaw.hpp
//
// Generic power-law hemolysis model:
// HI = C * t^alpha * tau_s^beta
//
// Date: Jul. 5 2026
// ============================================================================
#include "IHIModel.hpp"

class HIModel_PowerLaw : public IHIModel
{
  public:
    HIModel_PowerLaw(
        const std::string &in_model_name,
        const std::string &in_species_name,
        const std::string &in_reference_name,
        const double &in_C,
        const double &in_alpha,
        const double &in_beta )
    : model_name(in_model_name),
      species_name(in_species_name),
      reference_name(in_reference_name),
      C(in_C),
      alpha(in_alpha),
      beta(in_beta)
    {}

    virtual ~HIModel_PowerLaw() = default;

    void print_info() const override
    {
      SYS_T::commPrint("\t  HIModel_PowerLaw:: \n");
      SYS_T::commPrint("\t  Model     = %s \n", model_name.c_str());
      SYS_T::commPrint("\t  Species   = %s \n", species_name.c_str());
      SYS_T::commPrint("\t  Reference = %s \n", reference_name.c_str());
      SYS_T::commPrint("\t  C         = %.10e \n", C);
      SYS_T::commPrint("\t  alpha     = %.10e \n", alpha);
      SYS_T::commPrint("\t  beta      = %.10e \n", beta);
    }

    std::string get_model_name() const override { return model_name; }

    double get_C() const override { return C; }

    double get_alpha() const override { return alpha; }

    double get_beta() const override { return beta; }

    static std::unique_ptr<IHIModel> createPreset(const std::string &model_name)
    {
      if(model_name == "Giersiepen_1990")
        return SYS_T::make_unique<HIModel_PowerLaw>("Giersiepen_1990", "Human RBC",
            "Giersiepen et al., 1990", 3.6200e-7, 0.7850, 2.4160);

      if(model_name == "Heuser_Opitz_1980")
        return SYS_T::make_unique<HIModel_PowerLaw>("Heuser_Opitz_1980", "Porcine",
            "Heuser & Opitz, 1980", 1.8000e-6, 0.7650, 1.9910);

      if(model_name == "Zhang_2011")
        return SYS_T::make_unique<HIModel_PowerLaw>("Zhang_2011", "Ovine",
            "Zhang et al., 2011", 1.2280e-7, 0.6606, 1.9918);

      if(model_name == "Ding_2015_Porcine")
        return SYS_T::make_unique<HIModel_PowerLaw>("Ding_2015_Porcine", "Porcine",
            "Ding et al., 2015", 6.7010e-6, 0.2778, 1.0981);

      if(model_name == "Ding_2015_Human")
        return SYS_T::make_unique<HIModel_PowerLaw>("Ding_2015_Human", "Human",
            "Ding et al., 2015", 3.4580e-8, 0.2777, 2.0639);

      if(model_name == "Ding_2015_Bovine")
        return SYS_T::make_unique<HIModel_PowerLaw>("Ding_2015_Bovine", "Bovine",
            "Ding et al., 2015", 9.7720e-7, 0.2076, 1.4445);

      if(model_name == "Gesenhues_2016")
        return SYS_T::make_unique<HIModel_PowerLaw>("Gesenhues_2016", "Ovine",
            "Gesenhues et al., 2016", 2.3212e-6, 0.3300, 1.4924);

      return nullptr;
    }

  protected:
    const std::string model_name;
    const std::string species_name;
    const std::string reference_name;
    const double C;
    const double alpha;
    const double beta;

  private:
    HIModel_PowerLaw() = delete;
};

#endif
