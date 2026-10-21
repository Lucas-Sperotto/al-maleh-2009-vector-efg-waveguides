// SPDX-License-Identifier: MIT
#include "al_maleh/tm_validation.hpp"
#include <Eigen/Core>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
using namespace al_maleh;
namespace {
void check(bool yes,const char* m) {
 if(!yes)throw std::runtime_error(m);
}
void test_overlap() {
 const Rectangle rect{0,2,0,1};
 const auto nodes=make_regular_cloud(rect,5,5,1.65);
 const MLSEvaluator mls(nodes,1);
 const BackgroundQuadrature q{rect,4,4,3};
 const auto a=assemble_cartesian(mls,q);
 const auto spectrum=solve_generalized_symmetric(a.stiffness,a.mass);
 const auto result=validate_tm_modes(mls,q,a,spectrum,2,2);
 check(result.size()==4,"Expected TM11, TM21, TM12 and TM22");
 for(const auto& m:result) {
  check(m.nearest_index>=0 && m.strongest_index>=0,"Missing modal match");
  check(m.nearest_cluster_size>=1,"Empty eigenvalue cluster");
  check(std::isfinite(m.nearest_relative_error) &&
        m.nearest_relative_error>=0,"Invalid lambda error");
  check(m.nearest_overlap>=0 && m.nearest_overlap<=1,"Overlap range");
  check(m.nearest_cluster_overlap>=0 && m.nearest_cluster_overlap<=1,"Subspace overlap range");
  check(m.strongest_overlap>=0 && m.strongest_overlap<=1,"Best overlap range");
  check(std::isfinite(m.divergence_ratio)&&m.divergence_ratio>=0,"Divergence range");
  check(std::isfinite(m.normal_ratio)&&m.normal_ratio>=0,"Normal flux range");
 }
 MatchOptions invalid;invalid.minimum_shape_overlap=2;
 bool rejected=false;
 try{(void)validate_tm_modes(mls,q,a,spectrum,2,2,invalid);}
 catch(const std::invalid_argument&){rejected=true;}
 check(rejected,"Invalid overlap threshold not rejected");
 // A different integration order must still be well defined.
 const BackgroundQuadrature higher{rect,4,4,4};
 const auto ah=assemble_cartesian(mls,higher);
 const auto sh=solve_generalized_symmetric(ah.stiffness,ah.mass);
 const auto mh=validate_tm_modes(mls,higher,ah,sh,1,1);
 check(mh.size()==1 && std::isfinite(mh[0].nearest_overlap),
       "Higher-order quadrature failed");
}
}
int main() {
 try{
  test_overlap();
  std::cout<<"PASS: TM analytic curl overlaps, subspace and diagnostics\n";
 }catch(const std::exception& e){
  std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;
 }
}
